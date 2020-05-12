# amir.shahein@epfl.ch
# This script uses some material from Yu-Cheng Lin's Github project: https://github.com/linyc74/pwm_scan (project accessed Feb. 28th 2020)

import re
import numpy as np
import pandas as pd
from tqdm import tqdm
from dataclasses import dataclass, field
from itertools import dropwhile

class PWMScan(object):
    def __init__(self):
        """
        Object attributes:
        
            self.sequence: the (genome) sequence to be scanned
                numpy 1D array, dtype=np.int
        
            self.annot: the genome annotation
                pandas DataFrame
        
            self.hits: the scanning result hits
                pandas DataFrame
                
            self.reg_hits: the hits result of scanning regulatory sequence
                pandas DataFrame
                
            self.PWM_Kdref:  PWM in the format of Kd/Kdref, so that to obtain the Kd 
                of a given sequence, you can multiply all the Kd/Kdref factors based on 
                positions in the PWM, and then multiply this by the Kd of the consensus
            
            self.n_mer: number of bases in a motif
            
            self.unpClusterList: unprocessed cluster list, containing ALL clusters
            with potential sites, even if binding motif size is 9 and overlap is (-8)
        """
        self.sequence = None
        self.annot = None
        self.hits = None
        self.reg_hits = None
        self.PWM_Kdref = None
        self.n_mer = None
        self.unpClusterList = None
        self.ovlpClusterList = None
        
        
    def load_Kd_Kdref_pwm(self, filename, n_mer):
        """
        Args: 
            filename: corresponds to a file that contains a Kd/Kdref PWM (ex: from Swank 2019)
            n_mer: number of bases in the motif
        """
        PWM_Kdref = np.ones((4, n_mer), np.float) #Start with 1s, because Kdref/Kdref = 1
        
        with open(filename, 'r') as fh:
            PWM_Kdref = [[float(B) for B in line.split()] for line in fh] #process the text file into a 2D list
        
        PWM_Kdref = np.array(PWM_Kdref, dtype=float) #convert the 2D list into a 2D numpy array
        
        PWM_Kdref[PWM_Kdref < 1] = 1 #minimum value is 1
                
        self.PWM_Kdref = PWM_Kdref
        self.n_mer = n_mer
        
        return
    
    def load_sequence(self, seq, seqtype, seqPercent=100):
        """
        Args:
            seq: str, could be two things:
                (1) the fasta file name of raw sequence
                (2) the DNA sequence, containing only (A, C, G, T)
                    no space or other characters allowed
            seqtype: str, 'FASTA' if (1), 'RAW' if (2)
            seqPercent: percentage of the sequence to load... (standard genome seqs can be huge)
        """
        # If the unique characters in seq only contains the four DNA letters
        if seqtype == "RAW":
            self.sequence = self.__str_to_np_seq(seq[:round(len(seq)*seqPercent/100)]) #just taking sequence to seqPercent
            return

        elif seqtype == "FASTA":
            self.sequence = self.__parse_fasta(seq)
            if seqPercent != 100:
                self.sequence = self.sequence[:round(len(self.sequence)*seqPercent/100)] #just taking sequence to seqPercent (note this process could be done in a more efficient way...)
            self.sequence = self.__str_to_np_seq(self.sequence)
            
        else:
            print("You didn't input the seqtype; input a string saying 'FASTA' or 'RAW'")
            
        if self.sequence is None:
            print('The sequence has not been input properly, it is None...')
            return        
        
    def load_annotation(self, filename):
        """
        NOTE: Based on my current workflow, I don't use this for the annotations'
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
            for line in dropwhile(self.is_comment, fh): #https://cmdlinetips.com/2018/01/3-ways-to-read-a-file-and-skip-initial-comments-in-python/
                entry = line.strip().split('\t')

                # GTF format
                # 0        1       2        3      4    5      6       7      8
                # seqname, source, feature, start, end, score, strand, frame, attribute

                D['Start'].append(int(entry[3]))
                D['End'].append(int(entry[4]))
                D['Strand'].append(entry[6])

                attribute = entry[8]
                print(attribute)
                # Use regexp to get the gene_id and name
                # attribute is the string to search within, in green is what to search for, .group() retrieves (is) the match, 
                gene_id = re.search('gene_id\s".*?";', attribute).group(0)[9:-2]
                print(gene_id)
                name = re.search('product\s".*?";', attribute).group(0)[9:-2] #Note that this needs to be fixed, 
                                                                            #because only exons and CDS regions have 'product' 
                                                                            #attribute, but this led me to question which regions 
                                                                            #(e.g. gene, exon) I should be using anyway, which led to a new approach

                D['Gene ID'].append(gene_id)
                D['Name'].append(name)

        # Create a data frame
        self.annot = pd.DataFrame(D, columns=columns)
        
    def is_comment(self, s):
        """function that checks if a line
            starts with some character, here
            # to identify that the line is a comment
        """
        #return true if a line starts with #
        return s.startswith('#')

    def launch_scan(self, output_file=None, threshold=100, report_adjacent_genes=False,
                    promoter_length=500, use_genomic_GC=False):
        """
        This function controls the flow of the script when running a scan to generate the single binding site dataframe
        Args:
            output_file:
                str, the output excel filename

            threshold:
                float or int, threshold of the score (Kd/Kdref ratio) below which the sequence motif is retained

            report_adjacent_genes:
                boolean

            promoter_length:
                int
                If reporting adjacent genes, the promoter range within which
                the hit is defined as adjacent

            use_genomic_GC: UNUSED CURRENTLY
                boolean, whether to use the genomic GC content as
                the background GC frequency... 

        Returns:
            self.hits:
                a pandas data frame with the following fields:
                ['Score', 'Sequence', 'Start', 'End', 'Strand']

                Also there's an option to write the result self.hits
                in the output file in csv format.
        """

        if self.PWM_Kdref is None or self.sequence is None:
            print('Either the PWM or the SEQUENCE has not been input properly...')
            return

        self.hits = self.__pwm_scan(self.PWM_Kdref, self.sequence, threshold) # RUN THE SCAN!

        if report_adjacent_genes and not self.annot is None:
            self.__find_adjacent_genes(distance_range=promoter_length)

        if not output_file is None:
            self.hits.to_csv(output_file)

        return self.hits
    
    def launch_scan_multifasta(self, input_file, output_file=None, threshold=50,
                    length5Prime=1000, length3Prime=100):
        """
        This function is an alternative scan for multifasta files, designed for the use-case where
        the user extracts a set of regulatory sequences (e.g. promoters from EPDnew) in a .fa file,
        and wants to cluster scan them, while retaining the labels (>the first line comment) for
        each of the sequences in a column in the output dataframe.
        Args:
            input_file:
                file name of the .fa file with multiple FASTA formatted sequences
                
            output_file:
                str, the output excel filename    
                
            threshold:
                float or int, threshold of the score (Kd/Kdref ratio) below which the sequence motif is retained
                
            length5Prime:
                This assumes the user is using promoter regions starting at the TSS, with an upstream and downstream range
                length upstream of TSS 
            
            length3Prime:
                This assumes the user is using promoter regions starting at the TSS, with an upstream and downstream range
                length downstream of TSS (note that length5Prime 100 and length3Prime 100 would correspond to 201bp)
        """
        
        if self.PWM_Kdref is None:
            print('First input the Kd/Kdref PWM...')
            return
        
        with open(input_file, 'r') as fh:
            lines = fh.read().splitlines()
        
        ###################### Declare things to feed into pwm_scan
        n_mer = self.PWM_Kdref.shape[1]    # length (num of cols) of the weight matrix
        cols = np.arange(n_mer) # array of column indices for PWM from 0 to (n_mer-1)
        PWM_rc = self.PWM_Kdref[::-1, ::-1] # Reverse complementary PWM
        # Create an empty data frame
        colnames = ['Score', 'Sequence', 'Start', 'End', 'Strand', 'EPDnew_ID']
        self.reg_hits = pd.DataFrame({name:[] for name in colnames},
                            columns=colnames)
        #######################
        
        TSS_EPDnew_ID = lines.pop(0).split()[1] #remove the first line from lines, and access its EPDnew ID
        iStart = 0
        
        for i in tqdm(range(len(lines))): #iterate through the lines in the file
            
            if lines[i].startswith('>'):
                
                seq = ''.join(lines[iStart:i]) #pull out the promoter sequence corresponding to TSS_EPDnew_ID
                
                seq = self.__str_to_np_seq(seq)
                
                self.__pwm_scan_multifasta(self.PWM_Kdref, PWM_rc, seq, threshold, n_mer, cols, TSS_EPDnew_ID, length5Prime)
                
                TSS_EPDnew_ID = lines[i].split()[1]
                
                iStart = i+1
                
                

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
               'N':0, 'n':0,
               'M':0, 'm':0,
               'R':0, 'r':0} # M, R, and N should be very rare in a valid genome sequence so just arbitrarily assign them to A

        for i, base in enumerate(str_seq): #i is the index, base is the letter  in str_seq
            np_seq[i] = ref.get(base, 0) #ref is a dictionary, using get allows to use base as a key to get associated numeric value, with
                                         #a default value of 0.

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
        
        """Convert from FASTA file to a more RAW sequence compatible with the __str_to_np_seq function"""
        
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

    def __pwm_scan(self, PWM, seq, thres):
        """
        The core function that performs the PWM scan
        through the (genomic) sequence.
        """
        if isinstance(seq, str): #if somehow the genome sequence is still a string, encode it as 0123
            seq = __str_to_np_seq(seq)

        n_mer = PWM.shape[1]    # length (num of cols) of the weight matrix
        cols = np.arange(n_mer) # array of column indices for PWM from 0 to (n_mer-1)

        PWM_rc = PWM[::-1, ::-1] # Reverse complementary PWM

        # Create an empty data frame
        colnames = ['Score', 'Sequence', 'Start', 'End', 'Strand']

        hits = pd.DataFrame({name:[] for name in colnames},
                            columns=colnames)

        # The main loop that scans through the (genome) sequence
        for i in tqdm(range(len(seq) - n_mer + 1)):

            window = seq[i:(i+n_mer)] #pull out the sequence we're comparing the PWM against.

            # --- The most important line of code ---
            #     Use integer coding to index the correct score from column 0 to (n_mer-1)
            score = np.prod( PWM[window, cols] ) #indexing with two np arrays subscripts members as coordinates

            if score < thres: #append a new row in the dataframe, with details
                hits.loc[len(hits)] = [score                       , # Score
                                       self.__np_to_str_seq(window), # Sequence
                                       i + 1                       , # Start
                                       i + n_mer                   , # End
                                       '+'                         ] # Strand

            # --- The most important line of code ---
            #     Use integer coding to index the correct score from column 0 to (n_mer-1)
            score = np.prod( PWM_rc[window, cols] )

            if score < thres:
                hits.loc[len(hits)] = [score                       , # Score
                                       self.__np_to_str_seq(window), # Sequence
                                       i + 1                       , # Start
                                       i + n_mer                   , # End
                                       '-'                         ] # Strand

        return hits
    
    def __pwm_scan_multifasta(self, PWM, PWM_rc, seq, thresh, n_mer, cols, TSS_EPDnew_ID, length5Prime):
        """
        The core function that performs the PWM scan
        through the regulatory sequences.
        """
        
        # The main loop that scans through the (genome) sequence
        
        for i in range(len(seq) - n_mer + 1):
        
            window = seq[i:(i+n_mer)] #pull out the sequence we're comparing the PWM against.
            
            # --- The most important line of code ---
            #     Use integer coding to index the correct score from column 0 to (n_mer-1)
            score = np.prod( PWM[window, cols] ) #indexing with two np arrays subscripts members as coordinates
            
            if score < thresh: #append a new row in the dataframe, with details
                self.reg_hits.loc[len(self.reg_hits)] = [score     , # Score
                                       self.__np_to_str_seq(window), # Sequence
                                       i + -length5Prime           , # Distance of binding site's first base relative to TSS
                                       i -length5Prime + n_mer -1  , # End
                                       'C'                         , # Coding Strand
                                       TSS_EPDnew_ID               ] # EPDnew ID
                
            # --- The most important line of code ---
            #     Use integer coding to index the correct score from column 0 to (n_mer-1)
            score = np.prod( PWM_rc[window, cols] )
            
            if score < thresh:
                self.reg_hits.loc[len(self.reg_hits)] = [score     , # Score
                                       self.__np_to_str_seq(window), # Sequence
                                       i -length5Prime             , # Distance of binding site's last base relative to TSS
                                       i -length5Prime + n_mer -1  , # End
                                       'N'                         , # Non-coding Strand
                                       TSS_EPDnew_ID               ] # EPDnew ID

        
        
    def __find_adjacent_genes(self, distance_range):
        
        """
        NOTE: Based on my current workflow, I don't use this for the annotations
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

        # Outer for loop --- for each PWM hit: Can probably just change this to FOR EACH CLUSTERS HIT!
        for i in tqdm(range(len(self.hits))):

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
                    continue #move to the next loop iteration

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
                self.hits.loc[i, 'Gene ID'] = ' ; '.join(map(str, gene_id)) #Turns each element of the iterable 'gene_id' into a string, and then concatenates the elements with a ; in between them
                self.hits.loc[i, 'Name'] = ' ; '.join(map(str, gene_name))
                self.hits.loc[i, 'Distance']  = ' ; '.join(map(str, distance))


    def generate_clusters(self, maxGap):
        
        """
        maxGap: the maximum gap distance between binding sites that is allowed, for
        binding sites in the same cluster (before truncating a cluster). 
        """
        
        unpClusterList = []
        
        sitesTemp = [] #Declare temporary variables used to build clusters and transfer to unpClusterList 
        affinitiesTemp = []
        startPosTemp = []
        strandTemp = []
        
        for i in tqdm(range(len(self.hits.index)-1)): #For each index in self.hits
            #gapDist = self.hits.loc[i+1].Start - self.hits.loc[i].End -1   
         
            #if the distance to the next site is less than the maximum gap distance input, but greater than the minimum gap distance (biggest overlap allowed)
            #if #((gapDist <= maxGap) and (gapDist >= minGap)) -- note, it will make some sense to incorporate this at some point.
            if ((self.hits.loc[i+1].Start - self.hits.loc[i].End -1) <= maxGap):
                sitesTemp.append(self.hits.loc[i].Sequence)
                affinitiesTemp.append(self.hits.loc[i].Score)
                startPosTemp.append(self.hits.loc[i].Start)
                strandTemp.append(self.hits.loc[i].Strand)
            
            elif (sitesTemp) : 
                #if sitesTemp is not empty, and therefore we're on the last member of a cluster
                sitesTemp.append(self.hits.loc[i].Sequence) #append the last member of the cluster to the temporary variables
                affinitiesTemp.append(self.hits.loc[i].Score)
                startPosTemp.append(self.hits.loc[i].Start)
                strandTemp.append(self.hits.loc[i].Strand)     
                
                #Transform temporary variables into a Clusters object and store it as an element of the list unpClusterList
                unpClusterList.append(Clusters(sitesTemp, affinitiesTemp, startPosTemp, strandTemp))
                
                sitesTemp = [] #Clear temporary variables
                affinitiesTemp = []
                startPosTemp = []
                strandTemp = []
        
        if (maxGap <= -1):
            self.ovlpClusterList = unpClusterList
        else:
            self.unpClusterList = unpClusterList
            
    def generate_clusters_multifasta(self, regHitsDF, maxGap):
        
        """
        This is the alternative method that compiles single binding site hits into 
        clusters when you're using extracted regulatory sequences passed as input
        through a file with multiple FASTA seqs.
        maxGap: the maximum gap distance between binding sites that is allowed, for
        binding sites in the same cluster (before truncating a cluster). 
        """
        
        unpClusterList = []
        sitesTemp = [] #Declare temporary variables used to build clusters and transfer to unpClusterList 
        affinitiesTemp = []
        startPosTemp = []
        strandTemp = []
        
        
        for i in tqdm(range(len(regHitsDF.index)-1)): #iterate over self.hits
            #gapDist = self.hits.loc[i+1].Start - self.hits.loc[i].End -1   
         
            #if the distance to the next site is less than the maximum gap distance input, but greater than the minimum gap distance (biggest overlap allowed)
            #if #((gapDist <= maxGap) and (gapDist >= minGap)) -- note, it will make some sense to incorporate this at some point.
            if (((regHitsDF.loc[i+1].Start - regHitsDF.loc[i].End -1) <= maxGap) and (regHitsDF.loc[i+1].EPDnew_ID == regHitsDF.loc[i].EPDnew_ID)): #also restrict cluster hits to same promoter
                sitesTemp.append(regHitsDF.loc[i].Sequence)
                affinitiesTemp.append(regHitsDF.loc[i].Score)
                startPosTemp.append(regHitsDF.loc[i].Start)
                strandTemp.append(regHitsDF.loc[i].Strand)
            
            elif (sitesTemp) : 
                #else if sitesTemp is not empty, or EPDnew_ID is not ==, and therefore we're on the last member of a cluster:
                sitesTemp.append(regHitsDF.loc[i].Sequence) #append the last member of the cluster to the temporary variables
                affinitiesTemp.append(regHitsDF.loc[i].Score)
                startPosTemp.append(regHitsDF.loc[i].Start)
                strandTemp.append(regHitsDF.loc[i].Strand)
                EPDnew_ID = regHitsDF.loc[i].EPDnew_ID #also append the cluster's promoter region
                
                #Transform temporary variables into a Clusters object and store it as an element of the list unpClusterList
                unpClusterList.append(AnnotatedClusters(sitesTemp, affinitiesTemp, startPosTemp, strandTemp, 9, EPDnew_ID)) #Note: change hard-encoded 9 to self.n_mer
                
                sitesTemp = [] #Clear temporary variables
                affinitiesTemp = []
                startPosTemp = []
                strandTemp = []
        
        if (maxGap <= -1):
            self.ovlpClusterList = unpClusterList
        else:
            self.unpClusterList = unpClusterList
    
@dataclass    
class Clusters:
    bindingSiteSeqs: list
    siteAffinities: list
    startPos: list
    endPos: list = field(init=False)
    strand: list
    spacing: list = field(init=False)
    siteLength: int = 9
    
    def __post_init__(self):
        self.endPos = [(x+(self.siteLength-1)) for x in self.startPos]
        self.spacing = [(self.startPos[y+1] - self.endPos[y] -1) for y in range(len(self.startPos)-1)]

@dataclass
class AnnotatedClusters(Clusters): # AnnotatedClusters is a child class of Clusters, with the additional attribute geneID
    geneID: str = 'UNKNOWN' #This default value should actually never get triggered


















