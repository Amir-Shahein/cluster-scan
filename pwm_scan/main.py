import re
import numpy as np
import pandas as pd


class PWMScan(object):
    def __init__(self):
        """
        Object attributes:

            self.pwm: the position weight matrix
                numpy 2D array, dtype=np.float

            self.psm: the position score matrix
                numpy 2D array, dtype=np.float

            self.sequence: the (genome) sequence to be scanned
                numpy 1D array, dtype=np.int

            self.annot: the genome annotation
                pandas DataFrame

            self.hits: the scanning result hits
                pandas DataFrame
        """
        self.pwm = None
        self.psm = None
        self.sequence = None
        self.annot = None
        self.hits = None

    def load_pwm(self, filename):
        """
        Args:
            filename: The input file could be comma, tab, space delimited.
        """
        sites = []
        with open(filename, 'r') as fh:
            S = fh.read().split() # split by space, tab or line break
            for s in S:
                sites = sites + s.split(',') # further split by ','

        for i in range(1, len(sites)):
            if len(sites[i]) != len(sites[0]):
                print('Input sites not identical in length!')
                return

        self.pwm = self.__gen_pwm(sites)

    def load_sequence(self, seq):
        """
        Args:
            seq: str, could be two things:
                (1) the fasta file name
                (2) the DNA sequence, containing only (A, C, G, T)
                    no space or other characters allowed
        """
        # If the unique characters in seq only contains the four DNA letters
        if set(seq.lower()) == set(['a', 'c', 'g', 't']):
            self.sequence = self.__str_to_np_seq(seq)
            return

        # If seq is not a DNA sequence, then seq should be the fasta file name
        self.sequence = self.__parse_fasta(seq)
        if self.sequence is None:
            print('Not valid fasta format.')
            return
        self.sequence = self.__str_to_np_seq(self.sequence)

    def load_annotation(self, filename):
        """
        Args:
            filename:
                str, the name of the annotation file in GTF format
                The GTF file by default has start, end, strand and attribute
                The attribute field should have 'gene_id' and 'name'
        """

        # Column names for the pandas data frame
        columns = ['Start', 'End', 'Strand', 'Gene ID', 'Name']

        # Create an empty dictionary for the pandas data frame
        D = {key:[] for key in columns}

        with open(filename, 'r') as fh:
            for line in fh:
                entry = line.strip().split('\t')

                # GTF format
                # 0        1       2        3      4    5      6       7      8
                # seqname, source, feature, start, end, score, strand, frame, attribute

                D['Start'].append(int(entry[3]))
                D['End'].append(int(entry[4]))
                D['Strand'].append(entry[6])

                attribute = entry[8]

                # Use regexp to get the gene_id and name
                gene_id = re.search('gene_id\s".*?";', attribute).group(0)[9:-2]
                name = re.search('name\s".*?";', attribute).group(0)[6:-2]

                D['Gene ID'].append(gene_id)
                D['Name'].append(name)

        # Create a data frame
        self.annot = pd.DataFrame(D, columns=columns)

    def load_count_matrix(self, matrix):
        """
        Args:
            matrix: 2-D numpy array, dtype = np.int
                the matrix should have 4 rows with
                row 0, 1, 2, 3 corresponding to A, C, G ,T
        """
        assert(matrix.shape[0] == 4)

        # Create a new float array and +1 for pseudocount
        pwm = np.array(matrix, dtype=np.float) + 1

        # Sum 'downwards' across rows to get position-specific sum
        # Usually every position should have the same number of sum
        pwm = pwm / np.sum(pwm, axis=0)

        self.pwm = pwm

    def launch_scan(self, filename=None, threshold=10., report_adjacent_genes=True,
                    promoter_length=500, use_genomic_GC=False):
        """
        Args:
            filename:
                str, the output excel filename

            threshold:
                float or int, threshold of the score above which the sequence motif is retained

            report_adjacent_genes:
                boolean

            promoter_length:
                int
                If reporting adjacent genes, the promoter range within which
                the hit is defined as adjacent

            use_genomic_GC:
                boolean, whether to use the genomic GC content as
                the background GC frequency to calculate score matrix

        Returns:
            self.hits:
                a pandas data frome with the following fields:
                ['Score', 'Sequence', 'Start', 'End', 'Strand']

                Also there's an option to write the result self.hits
                in the output file in csv format.
        """

        if self.pwm is None or self.sequence is None:
            return

        if use_genomic_GC:
            GC_content = self.__calc_GC(self.sequence)
            self.psm = self.__pwm_to_psm(self.pwm, GC_content)
        else:
            self.psm = self.__pwm_to_psm(self.pwm)

        self.hits = self.__psm_scan(self.psm, self.sequence, threshold)

        if report_adjacent_genes and not self.annot is None:
            self.__find_adjacent_genes(distance_range=promoter_length)

        if not filename is None:
            self.hits.to_csv(filename)

        return self.hits

    def __str_to_np_seq(self, str_seq):
        """
        A custom DNA base coding system with numbers.
        (A, C, G, T, N) = (0, 1, 2, 3, 0)

        A DNA string is converted to a numpy integer array (np.unit8) with the same length.
        """
        np_seq = np.zeros((len(str_seq), ), np.uint8)

        ref = {'A':0, 'a':0,
               'C':1, 'c':1,
               'G':2, 'g':2,
               'T':3, 't':3,
               'N':0, 'n':0} # N should be very rare in a valid genome sequence so just assign it as A

        for i, base in enumerate(str_seq):
            np_seq[i] = ref[base]

        return np_seq

    def __np_to_str_seq(self, np_seq):
        """
        Convert (0, 1, 2, 3, 4) base coding back to (A, C, G, T, N)
        """
        str_seq = ['A' for i in range(len(np_seq))]

        ref = {0:'A',
               1:'C',
               2:'G',
               3:'T'}

        for i, num in enumerate(np_seq):
            str_seq[i] = ref[num]

        return ''.join(str_seq)

    def __parse_fasta(self, filename):
        with open(filename, 'r') as fh:
            lines = fh.read().splitlines()
        first = lines.pop(0)
        if first.startswith('>'):
            return ''.join(lines)
        else:
            return None

    def __rev_comp(self, seq):
        """
        Reverse complementary. Input could be either a numpy array or a DNA string.
        """
        if isinstance(seq, np.ndarray):
            return 3 - seq[::-1]
        elif isinstance(seq, str):
            seq = 3 - __str_to_np_seq(seq)[::-1]
            return __np_to_str_seq(seq)

    def __calc_GC(self, seq):
        """
        Calculate GC content. Input could be either a numpy array or a DNA string.
        """
        if isinstance(seq, str):
            seq = self.__str_to_np_seq(seq)

        if isinstance(seq, np.ndarray):
            GC_count = np.sum(seq == 1) + np.sum(seq == 2)
            GC_content = GC_count / float(len(seq))
            return GC_content

    def __gen_pwm(self, sites):
        """
        Takes a list of sites (str with identical length).
        Returns the position weight matrix.

        Args:
            sites: list of str

        Returns:
            pwm: 2-D numpy arraym, dtype = np.float
        """

        sites = [self.__str_to_np_seq(s) for s in sites]

        # Start to build position weight matrix = pwm
        # Rows 0, 1, 2, 3 correspond to A, C, G, T, respectively, which is identical to the integer coding.
        n_mer = len(sites[0])
        pwm = np.ones((4, n_mer), np.float) # The position count matrix begin from 1 as pseudocount

        for s in sites:
            for i, base in enumerate(s): # For each base, add to the count matrix.
                pwm[base, i] += 1        # The integer coding (0, 1, 2, 3) = (A, C, G, T)
                                         # corresponds to the order of rows 0, 1, 2, 3 = A, C, G, T

        pwm = pwm / np.sum(pwm[:, 1]) # Compute the position weight matrix, i.e. probability matrix

        return pwm

    def __pwm_to_psm(self, pwm, GC_content=0.5):
        """
        Converts position weight matrix to position score matrix.
        The background GC content is taken into account to comput the likelihood.
        Score is the log2 likelihood.

        The default background GC_content = 0.5, which is strongly recommended.
        """
        psm = np.zeros(pwm.shape, np.float)

        psm[[0,3], :] = pwm[[0,3], :] / ((1 - GC_content) / 2) # Divide by the background frequency -> likelihood matrix
        psm[[1,2], :] = pwm[[1,2], :] / ((    GC_content) / 2)

        # log2 likelihood -> score matrix
        return np.log2(psm)

    def __psm_scan(self, psm, seq, thres):
        """
        The core function that performs PWM (PSM) scan
        through the (genomic) sequence.
        """
        if isinstance(seq, str):
            seq = __str_to_np_seq(seq)

        n_mer = psm.shape[1]    # length (num of cols) of the weight matrix
        cols = np.arange(n_mer) # column indices for psm from 0 to (n_mer-1)

        psm_rc = psm[::-1, ::-1] # Reverse complementary psm

        # Create an empty data frame
        colnames = ['Score', 'Sequence', 'Start', 'End', 'Strand']

        hits = pd.DataFrame({name:[] for name in colnames},
                            columns=colnames)

        # The main loop that scans through the (genome) sequence
        for i in range(len(seq) - n_mer + 1):

            window = seq[i:(i+n_mer)]

            # --- The most important line of code ---
            #     Use integer coding to index the correct score from column 0 to (n_mer-1)
            score = np.sum( psm[window, cols] )

            if score > thres:
                hits.loc[len(hits)] = [score                       , # Score
                                       self.__np_to_str_seq(window), # Sequence
                                       i + 1                       , # Start
                                       i + n_mer                   , # End
                                       '+'                         ] # Strand

            # --- The most important line of code ---
            #     Use integer coding to index the correct score from column 0 to (n_mer-1)
            score = np.sum( psm_rc[window, cols] )

            if score > thres:
                hits.loc[len(hits)] = [score                       , # Score
                                       self.__np_to_str_seq(window), # Sequence
                                       i + 1                       , # Start
                                       i + n_mer                   , # End
                                       '-'                         ] # Strand

        return hits

    def __find_adjacent_genes(self, distance_range):
        """
        Args:
            distance_range: distance of promoter range in bp

        Returns:
            None. This method modifies self.hits
        """

        # self.hits is a pandas data frame that contains the PWM hit sites

        # Adding three new fields (columns)
        for h in ('Gene ID', 'Name', 'Distance'):
            self.hits[h] = ''

        # Convert self.annot (data frame) into a list of dictionaries for faster search
        # list of dictionaries [{...}, {...}, ...]
        # Each dictionary has 'Start', 'End', 'Strand'
        intervals = []
        for j in range(len(self.annot)):
            entry = self.annot.iloc[j, :]
            intervals.append({'Start' :entry['Start'] ,
                              'End'   :entry['End']   ,
                              'Strand':entry['Strand']})

        # Outer for loop --- for each PWM hit
        for i in range(len(self.hits)):

            # ith row -> pandas series
            hit = self.hits.iloc[i, :]

            # The hit location is the mid-point of Start and End
            hit_loc = int( (hit['Start'] + hit['End']) / 2 )

            # Search through all annotated genes to see if
            # the hit location lies in the UPSTREAM promoter of any one of the genes

            # Create empty lists for each hit
            gene_id = []
            gene_name = []
            distance = []
            # Inner for loop --- for each annotated gene
            for j, intv in enumerate(intervals):

                if hit_loc < intv['Start'] - 500:
                    continue

                if hit_loc > intv['End'] + 500:
                    continue

                if intv['Strand'] == '+':
                    promoter_from = intv['Start'] - distance_range
                    promoter_to   = intv['Start'] - 1

                    if hit_loc >= promoter_from and hit_loc <= promoter_to:
                        entry = self.annot.iloc[j, :]
                        gene_id.append(entry['Gene ID'])
                        gene_name.append(entry['Name'])
                        dist = intv['Start'] - hit_loc
                        distance.append(dist)

                elif intv['Strand'] == '-':
                    promoter_from = intv['End'] + 1
                    promoter_to   = intv['End'] + distance_range

                    if hit_loc >= promoter_from and hit_loc <= promoter_to:
                        entry = self.annot.iloc[j, :]
                        gene_id.append(entry['Gene ID'])
                        gene_name.append(entry['Name'])
                        dist = hit_loc - intv['End']
                        distance.append(dist)

            # If the gene_id != []
            if gene_id:
                self.hits.loc[i, 'Gene ID'] = ' ; '.join(map(str, gene_id))
                self.hits.loc[i, 'Name'] = ' ; '.join(map(str, gene_name))
                self.hits.loc[i, 'Distance']  = ' ; '.join(map(str, distance))


