#-------------------------------------------------------------------------------
# Version: 1.1

# Author: Yu-Cheng Lin

# Modifications:
#   Introduce pandas Data Frame
#-------------------------------------------------------------------------------

import numpy as np
import pandas as pd
from Bio import SeqIO



class PWMScan(object):
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

    # Publich methods

    def load_pwm(self, filename):
        '''
        :parameter
            filename: The input file could be comma, tab, space delimited.
        '''
        sites = []
        with open(filename, 'r') as fh:
            S = fh.read().split()
            for s in S:
                sites = sites + s.split(',')

        for i in xrange(1, len(sites)):
            if len(sites[i]) != len(sites[0]):
                print 'Input sites not identical in length!'
                return

        self.pwm = self.__gen_pwm(sites)

    def load_sequence(self, filename):
        '''
        :parameter
            filename: str, the name of the fasta file
        '''
        parser = SeqIO.parse(filename, 'fasta')
        self.sequence = str(next(parser).seq)
        self.sequence = self.__str_to_np_seq(self.sequence)
        # self.sequence = self.sequence[1:100000]

    def load_annotation(self, filename):
        '''
        :parameter
            filename: str, the name of the annotation file in csv format
        '''
        self.annot = pd.read_csv(filename)

    def launch_scan(self, filename, thres, report_adjacent_genes=True, use_genomic_GC=False):
        '''
        :parameter
            filename:
                the output filename

            thres:
                threshold of the score above which the sequence is retained

            report_adjacent_genes:
                True/False

            use_genomic_GC:
                True/False. Decides whether to use the genomic GC content as
                the background GC frequency to calculate score matrix

        :return:
            None. But write the results in the output file in csv format.
        '''

        if self.pwm is None or self.sequence is None:
            return

        if use_genomic_GC:
            GC_content = self.__calc_GC(self.sequence)
            self.psm = self.__pwm_to_psm(self.pwm, GC_content)
        else:
            self.psm = self.__pwm_to_psm(self.pwm)

        self.hits = self.__psm_scan(self.psm, self.sequence, thres)

        if report_adjacent_genes and not self.annot is None:
            self.__find_adjacent_genes(distance_range=500)

        self.hits.to_csv('output.csv')

    # Private methods

    def __str_to_np_seq(self, str_seq):
        '''
        A custom DNA base coding system with numbers.
        (A, C, G, T, N) = (0, 1, 2, 3, 0)

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

    def __np_to_str_seq(self, np_seq):
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

    def __rev_comp(self, seq):
        '''
        Reverse complementary. Input could be either a numpy array or a DNA string.
        '''
        if isinstance(seq, np.ndarray):
            return 3 - seq[::-1]
        elif isinstance(seq, str):
            seq = 3 - __str_to_np_seq(seq)[::-1]
            return __np_to_str_seq(seq)

    def __calc_GC(self, seq):
        '''
        Calculate GC content. Input could be either a numpy array or a DNA string.
        '''
        if isinstance(seq, str):
            seq = self.__str_to_np_seq(seq)

        if isinstance(seq, np.ndarray):
            GC_count = np.sum(seq == 1) + np.sum(seq == 2)
            GC_content = GC_count / float(len(seq))
            return GC_content

    def __gen_pwm(self, sites):
        '''
        Takes a list of sites (str with identical length).
        Returns the position weight matrix.
        '''

        sites = [self.__str_to_np_seq(s) for s in sites]

        # Start to build position weight matrix = pwm
        # Rows 0, 1, 2, 3 correspond to A, C, G, T, respectively, which is identical to the integer coding.
        n_mer = len(sites[0])
        pwm = np.ones((4, n_mer), np.float) # The position count matrix begin from 1 as pseudocount

        for s in sites:
            for i, base in enumerate(s): # For each base, add to the count matrix.
                pwm[base, i] += 1        # The integer coding (0, 1, 2, 3) = (A, C, G, T)
                                         # corresponds the order of rows 0, 1, 2, 3 = A, C, G, T

        pwm = pwm / np.sum(pwm[:, 1]) # Compute the position weight matrix, i.e. probability matrix

        return pwm

    def __pwm_to_psm(self, pwm, GC_content=0.5):
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

    def __psm_scan(self, psm, seq, thres):
        '''
        The core function that performs PWM (PSM) scan
        through the (genomic) sequence.
        '''
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
        for i in xrange(len(seq) - n_mer + 1):

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

        # self.hits is a pandas data frame that contains the PWM derected sites

        # Adding new fields (columns)
        for h in ('Locus Tag', 'Gene Name', 'Distance'):
            self.hits[h] = ''

        # list of dictionaries [{...}, {...}, ...]
        # Each dictionary has 'Start', 'End', 'Strand'
        intervals = []
        for j in xrange(len(self.annot)):
            gene = self.annot.iloc[j, :]
            intervals.append({'Start' :gene['Start'] ,
                              'End'   :gene['End']   ,
                              'Strand':gene['Strand']})

        for i in xrange(len(self.hits)):
            print i

            # ith row -> pandas series
            hit = self.hits.iloc[i, :]

            # The hit location is the mid-point of Start and End
            hit_loc = int( (hit['Start'] + hit['End']) / 2 )

            # Search through all annotated genes to see if
            # the hit location lies in the UPSTREAM promoter of any one of the genes

            locus_tag = []
            gene_name = []
            distance = []
            for j, intv in enumerate(intervals):

                if intv['Strand'] == '+':
                    promoter_from = intv['Start'] - distance_range
                    promoter_to   = intv['Start'] - 1

                    if hit_loc >= promoter_from and hit_loc <= promoter_to:
                        gene = self.annot.iloc[j, :]
                        dist = intv['Start'] - hit_loc
                        locus_tag.append(gene['Locus Tag'])
                        gene_name.append(gene['Name'])
                        distance.append(dist)

                elif intv['Strand'] == '-':
                    promoter_from = intv['End'] + 1
                    promoter_to   = intv['End'] + distance_range

                    if hit_loc >= promoter_from and hit_loc <= promoter_to:
                        gene = self.annot.iloc[j, :]
                        dist = hit_loc - intv['End']
                        locus_tag.append(gene['Locus Tag'])
                        gene_name.append(gene['Name'])
                        distance.append(dist)

            if locus_tag:
                self.hits.loc[i, 'Locus Tag'] = ' | '.join(map(str, locus_tag))
                self.hits.loc[i, 'Gene Name'] = ' | '.join(map(str, gene_name))
                self.hits.loc[i, 'Distance']  = ' | '.join(map(str, distance ))



if __name__ == '__main__':

    pwmscan = PWMScan()

    # Load binding sites to generate a position weight matrix
    pwmscan.load_pwm('Anr_sites.txt')

    # Load genome sequence (fasta)
    pwmscan.load_sequence('NC_008463_PA14_genome.fna')

    # Load annotation (csv)
    pwmscan.load_annotation('PA14_annotation.csv')

    pwmscan.launch_scan('output.csv', 10)


