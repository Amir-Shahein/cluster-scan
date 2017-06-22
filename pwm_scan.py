import numpy as np
import timeit, csv
from Bio import SeqIO
import pandas



class PositionWeightMatrix(object):
    def __init__(self):
        '''
        List all data entities for performing a pwm scan:
            > The position weight matrix (and score matrix)
            > The (genome) sequence to be scanned
            > The genome annotation as a data frame
            > The scanning result hits as a data frame
        '''
        self.pwm = None
        self.psm = None
        self.sequence = None
        self.annot = None
        self.hits = None

    # Higher-level methods

    def load_pwm(self, csv_filename):
        if not csv_filename[-4:] == '.csv':
            return
        self.pwm = self.gen_pwm(csv_filename)

    def load_sequence(self, fasta_filename):
        if not fasta_filename[-4:] == '.fna' and not fasta_filename[-6:] == '.fasta':
            return
        fasta_reader = SeqIO.parse(fasta_filename, 'fasta')
        self.sequence = str(next(fasta_reader).seq)
        self.sequence = self.str_to_np_seq(self.sequence)

    def load_annotation(self, csv_filename):
        if not csv_filename[-4:] == '.csv':
            return
        self.annot = self.read_csv_as_data_frame(csv_filename)

    def launch_scan(self, csv_filename, thres, report_adjacent_genes=True, use_genomic_GC=False):
        # use_genomic_GC: decides whether to use the genomic GC content as the background GC frequency to calculate score matrix
        if self.pwm is None or self.sequence is None:
            return

        if use_genomic_GC:
            self.psm = self.pwm_to_psm(self.pwm, GC_content=self.calc_GC_content(self.sequence))
        else:
            self.psm = self.pwm_to_psm(self.pwm)

        self.hits = self.psm_scan(self.psm, self.sequence, thres)

        if report_adjacent_genes:
            self.find_adjacent_genes(distance_range=500)

        self.write_data_frame_to_csv(data_frame=self.hits, filename=csv_filename)

    # Lower-level methods

    def str_to_np_seq(self, str_seq):
        '''
        A custom DNA base coding system with numbers.
        (A, C, G, T, N) = (0, 1, 2, 3, 4)

        A DNA string is converted to a numpy integer array (np.unit8) with the same length.
        '''
        np_seq = np.zeros((len(str_seq), ), np.uint8)

        ref = {'A':0, 'a':0,
               'C':1, 'c':1,
               'G':2, 'g':2,
               'T':3, 't':3,
               'N':0, 'n':0} # N should be very rare in a valid genome sequence so just assign it as A

        for i, base in enumerate(str_seq):
            np_seq[i] = ref[base]

        return np_seq

    def np_to_str_seq(self, np_seq):
        '''
        Convert (0, 1, 2, 3, 4) base coding back to (A, C, G, T, N)
        '''
        str_seq = ['A' for i in xrange(len(np_seq))]

        ref = {0:'A',
               1:'C',
               2:'G',
               3:'T'}

        for i, num in enumerate(np_seq):
            str_seq[i] = ref[num]

        return ''.join(str_seq)

    def rev_comp(self, seq):
        '''
        Reverse complementary. Input could be either a numpy array or a DNA string.
        '''
        if isinstance(seq, np.ndarray):
            return 3 - seq[::-1]
        elif isinstance(seq, str):
            seq = 3 - str_to_np_seq(seq)[::-1]
            return np_to_str_seq(seq)

    def calc_GC_content(self, seq):
        '''
        Calculate GC content. Input could be either a numpy array or a DNA string.
        '''
        if isinstance(seq, str):
            seq = self.str_to_np_seq(seq)

        if isinstance(seq, np.ndarray):
            GC_count = np.sum(seq == 1) + np.sum(seq == 2)
            GC_content = GC_count / float(len(seq))
            return GC_content

    def gen_pwm(self, csv_filename):
        '''
        Takes an csv file containing all the binding sites (DNA string ACGT) with identical length.
        Returns the position weight matrix.
        '''
        binding_sites = []
        with open(csv_filename, 'rb') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')

            for line in csv_reader: # Each line is parsed as a list
                for word in line: # Each word is parsed as a str
                    word = self.str_to_np_seq(word) # Convert to integer coding: (A, C, G, T) -> (0, 1, 2, 3)
                    binding_sites.append(word)

        # Remove any binding sites in the list not equal to the length of the first binding site
        valid_sites = []
        for site in binding_sites:
            if len(site) == len(binding_sites[0]):
                valid_sites.append(site)
        binding_sites = valid_sites

        # Start to build position weight matrix = pwm
        # Rows 0, 1, 2, 3 correspond to A, C, G, T, respectively, which is identical to the integer coding.
        n_mer = len(binding_sites[0])
        pwm = np.ones((4, n_mer), np.float) # The position count matrix begin from 1 as pseudocount

        for each_site in binding_sites:
            for i, base in enumerate(each_site): # For each base, add to the count matrix.
                pwm[base, i] += 1                # The integer coding (0, 1, 2, 3) = (A, C, G, T)
                                                 # corresponds the order of rows 0, 1, 2, 3 = A, C, G, T

        pwm = pwm / np.sum(pwm[:, 1]) # Compute the position weight matrix, i.e. probability matrix

        return pwm

    def pwm_to_psm(self, pwm, GC_content=0.5):
        '''
        Converts position weight matrix to position score matrix.
        The background GC content is taken into account to comput the likelihood.
        Score is the log2 likelihood.

        The default background GC_content = 0.5, which is strongly recommended.
        '''
        psm = np.zeros(pwm.shape, np.float)

        psm[[0,3], :] = pwm[[0,3], :] / ((1 - GC_content) / 2) # Divide by the background frequency -> likelihood matrix
        psm[[1,2], :] = pwm[[1,2], :] / ((    GC_content) / 2)

        # log2 likelihood -> score matrix
        return np.log2(psm)

    def psm_scan(self, psm, seq, thres):
        '''
        The core function that performs position-weight-matrix scan through the (genomic) sequence.
        '''
        if isinstance(seq, str):
            seq = str_to_np_seq(seq)

        n_mer = psm.shape[1]    # length (num of cols) of the weight matrix
        cols = np.arange(n_mer) # column indices for psm from 0 to (n_mer-1)

        psm_rc = psm[::-1, ::-1] # Reverse complementary psm

        # The list 'hits' has the data-frame structure explained in the method 'self.read_csv_as_data_frame()'
        # hits is a list of dictionaries,
        # with each dictionary representing a sample
        hits = []

        # The main loop that scans through the (genome) sequence
        for i in xrange(len(seq) - n_mer + 1):

            window = seq[i:(i+n_mer)]

            # If there is any "N" in the current window sequence, then skip
            if np.sum(window == 4) > 0:
                continue

            # The most important line of code, use integer coding to index the correct score from column 0 to (n_mer-1)
            score = np.sum( psm[window, cols] )

            if score > thres:
                hits.append({'Score'   : score                     ,
                             'Sequence': self.np_to_str_seq(window),
                             'Start'   : i + 1                     ,
                             'End'     : i + n_mer                 ,
                             'Strand'  : '+'                       })

            # The most important line of code, again, with the reverse complementary psm
            score = np.sum( psm_rc[window, cols] )

            if score > thres:
                hits.append({'Score'   : score                     ,
                             'Sequence': self.np_to_str_seq(window),
                             'Start'   : i + 1                     ,
                             'End'     : i + n_mer                 ,
                             'Strand'  : '-'                       })

            # Show the current progress
            if i % 100000 == 0:
                print 'Currently at sequence location ' + str(i) + ' with ' + str(len(hits)) + ' hits'

        return hits

    def read_csv_as_data_frame(self, filename, header=True):
        '''
        Here I want to define a data structure to mimic 'data frame' in R.

        In R each row of a data frame is a sample. The order of the sample could matter.
        That's why sometimes we sort the order of samples (rows) based on some kind of ranking.
        In R each column of a data frame is an attribute of samples.
        Usually the order of columns does not matter, as attributes are usually not ranked
        and attributes are not more or less impportant than one another.

        Here in python I use a list containing dictionaries to mimic the data frame structure in R.
        Each item in the list is a sample. The content of each item contains data associated
        with that particular sample. Therefore each item in the list corresponds to a row in an R data frame.
        The content of each item in the list is organized as a dictionary,
        with keys as attribute names, and values as attribute values.

        With the above structure - a list containing dictionaries, which are samples -
        samples are rankable (because they are items of a list),
        but attributes are not and need not be rankable.
        Attribute values can be easily accessed because of dictionary.

        In R it is: data_frame$attribute_name[sample_id] or
                    data_frame[sample_id, 'attribute_name']

        Here it is: list[sample_id]['attribute_name']
        '''
        data_frame = []

        with open(filename, 'r') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')

            if header == True:
                headers = next(csv_reader) # csv_reader is a generator, use next() to take the first line as a list of headers
            else:
                headers = None

            for line in csv_reader: # Each line is parsed as a list
                if headers == None:
                    headers = range(len(line)) # Without headers simply use a integer list as headers

                sample = {} # Each line = a sample = a dictionary

                for i, word in enumerate(line): # Each word in a line is parsed as a str
                    sample[ headers[i] ] = word # Append an attribute to the dictionary, with attribute name = headers[i] and attribate value = word.

                data_frame.append(sample)

        return data_frame

    def write_data_frame_to_csv(self, data_frame, filename):
        '''
        The data structure of data frame is defined in read_csv_as_data_frame()
        '''
        attr_names = data_frame[0].keys()

        with open(filename, mode='w') as csv_file:
            # Write the attribute names (headers) in the first row
            for name in attr_names:
                csv_file.write(name + ',')
            csv_file.write('\n')

            for sample in data_frame:
                for name in attr_names:
                    if isinstance(sample[name], list):
                        str_list = map(str, sample[name]) # Convert every element in the list to str type
                        word = ' | '.join(str_list) + ','
                    else:
                        word = str(sample[name]) + ','
                    csv_file.write(word)
                csv_file.write('\n') # Change line for each sample

    def sort_hits(self, attribute, dataframe=None):
        '''
        This sorting function does not alter the data frame self.hits.
        It returns the sorted data frame.
        '''

        if dataframe == None:
            # The default data frame is self.hits
            dataframe = self.hits

        # Pull out the values of the particular attribute to form a list
        L = [ int(data[attribute]) for data in dataframe ]

        # Get the sorted indices
        indices = np.argsort(L)

        sorted_df = [None] * len(dataframe)
        for i, id in enumerate(indices):
            # Put in the values according to the sorted order
            sorted_df[i] = dataframe[id]

        return sorted_df

    def write_np_array_to_csv(self, np_array, filename):
        '''
        Writes a 2D numpy array (usually the position-weight matrix) to a csv file
        '''

        with open(filename, mode='w') as csv_file:
            rowN, colN = np_array.shape
            for r in xrange(rowN):
                for c in xrange(colN):
                    word = '%s%s' % (str(np_array[r, c]), ',')
                    csv_file.write(word)
                csv_file.write('\n') # Change line after completing a row of np_array

    def find_adjacent_genes(self, distance_range):

        # self.hits is the data frame structure defined in the method self.read_csv_as_data_frame()
        # which is a list containing dictionaries
        for i, hit in enumerate(self.hits):

            # Show current progress
            if i % 100 == 0:
                print 'Finding adjacent genes for the ' + str(i) + 'th hit '

            # For each hit, i.e. a dictionary, add three new attributes
            for new_attr in ['Locus Tag', 'Gene Name', 'Distance']:
                hit[new_attr] = []

            hit_location = (hit['Start'] + hit['End']) / 2

            # Search through all annotated genes to see if
            # the hit location lies in the UPSTREAM promoter of any one of the genes
            for gene in self.annot:
                if gene['Strand'] == '+':
                    promoter_from = int(gene['Start']) - distance_range
                    promoter_to   = int(gene['Start']) - 1

                    if hit_location >= promoter_from and hit_location <= promoter_to:
                        hit['Locus Tag'].append(gene['Locus Tag'])
                        hit['Gene Name'].append(gene['Name'])
                        distance = int(gene['Start']) - hit_location
                        hit['Distance'].append(distance)

                elif gene['Strand'] == '-':
                    promoter_from = int(gene['End']) + 1
                    promoter_to   = int(gene['End']) + distance_range

                    if hit_location >= promoter_from and hit_location <= promoter_to:
                        hit['Locus Tag'].append(gene['Locus Tag'])
                        hit['Gene Name'].append(gene['Name'])
                        distance = hit_location - int(gene['End'])
                        hit['Distance'].append(distance)



if __name__ == '__main__':

    PWM = PositionWeightMatrix()

    # Load genome sequence (fasta format)
    PWM.load_sequence('NC_008463_PA14_genome.fna')

    # Load annotation (csv format)
    PWM.load_annotation('PA14_annotation.csv')



    ''' '''

    # Load -10 sites stored in a csv file
    PWM.load_pwm('16_0900 RpoS PWM_02_minus ten.csv')

    # Scan with the -10 position weight matrix
    PWM.launch_scan('scan_results_minus10.csv',
                    thres=5, report_adjacent_genes=True, use_genomic_GC=False)

    # Store the scan results in a new variable (to be used later)
    minus_10_hits = [hit for hit in PWM.hits]



    ''' '''

    # Load -35 sites stored in a csv file
    PWM.load_pwm('16_0900 RpoS PWM_02_minus thirtyfive.csv')

    # Scan with the -35 position weight matrix
    PWM.launch_scan('scan_results_minus35.csv',
                    thres=6, report_adjacent_genes=True, use_genomic_GC=False)

    # Store the scan results in a new variable (to be used later)
    minus_35_hits = [hit for hit in PWM.hits]



    ''' '''

    # Combine the two data frames (row concatenation)
    minus_10_35_hits = minus_10_hits + minus_35_hits

    # Sort by the 'Start' position of each hit
    sorted_df = PWM.sort_hits(attribute='Start', dataframe=minus_10_35_hits)

    # Save as csv
    PWM.write_data_frame_to_csv(sorted_df, 'scan_results_minus10and35.csv')





    # all_hits = PWM.read_csv_as_data_frame('scan_results_minus10and35.csv')
    #
    # fasta_reader = SeqIO.parse('NC_008463_PA14_genome.fna', 'fasta')
    # PA14genome = str(next(fasta_reader).seq)
    #
    # RpoS_targets = []
    #
    # for i in xrange(1000):
    #     hit_1 = all_hits[i]
    #     hit_2 = all_hits[i+1]
    #
    #     if hit_1['Strand'] == '+' and hit_2['Strand'] == '+':
    #
    #         if len(hit_1['Sequence']) == 10 and len(hit_2['Sequence']) == 13:
    #
    #             d = int(hit_2['Start']) - int(hit_1['End'])
    #             if d >= 8 and d <= 50:
    #
    #                 start_1 = int(hit_1['Start'])
    #                 end_1 = int(hit_1['End'])
    #                 start_2 = int(hit_2['Start'])
    #                 end_2 = int(hit_2['End'])
    #
    #                 new = {'Strand'                  : '+',
    #                        'Location'                : int(hit_1['Start']),
    #                        'Sequence before -35'     : PA14genome[(start_1-10):start_1],
    #                        'Sequence -35'            : hit_1['Sequence'],
    #                        'Sequence between -10 -35': PA14genome[end_1:(start_2+1)],
    #                        'Sequence -10'            : hit_2['Sequence'],
    #                        'Locus Tag'               : hit_2['Locus Tag'],
    #                        'Gene Name'               : hit_2['Gene Name'],
    #                        'Distance'                : hit_2['Distance'],
    #                        'Score -35'               : hit_1['Score'],
    #                        'Score -10'               : hit_2['Score']}
    #
    #                 RpoS_targets.append(new)
    #
    # PWM.write_data_frame_to_csv(RpoS_targets, 'RpoS_targets.csv')









