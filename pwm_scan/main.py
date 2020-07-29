# amir.shahein@epfl.ch
# This script uses a bit of material from Yu-Cheng Lin's Github project: https://github.com/linyc74/pwm_scan (project accessed Feb. 28th 2020)

import re
import numpy as np
import pandas as pd
from tqdm import tqdm
from dataclasses import dataclass, field
from itertools import dropwhile
import matplotlib.pyplot as plt 
from statMechModel import general_cluster_statmech_model
from statMechModel import single_site_statmech_model
from mpl_toolkits.mplot3d.axes3d import Axes3D

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
                
            self.restricted_reg_hits: reg_hits restsricted to an annotation (
            designed for annotation: ChIP-seq hit <=1kB from TSS --> i.e. bound promoter)
                
            self.PWM_Kdref:  PWM in the format of Kd/Kdref, so that to obtain the Kd 
                of a given sequence, you can multiply all the Kd/Kdref factors based on 
                positions in the PWM, and then multiply this by the Kd of the consensus
            
            self.n_mer: number of bases in a motif
                int
            
            self.unpClusterList: unprocessed cluster list, containing ALL clusters
            with potential sites, even if binding motif size is 9 and overlap is (-8)
                list of objects
            
            self.regElementList: list of RegulatoryElement objects
                list of objects
            
            self.flat_df_count: dataframe corresponding to spacings between 
            binding sites and the number of occurrences of these spacings in the
            input dataframe
                
            self.non_clus_reg_hits: dataframe of individual binding sites, i.e.
            binding sites that are not part of a homotypic cluster
        """
        self.sequence = None
        self.annot = None
        self.hits = None
        self.reg_hits = None
        self.restricted_reg_hits = None
        self.EPDnewIDchipValuesDict = None
        self.PWM_Kdref = None
        self.n_mer = None
        self.EPDnewIDconversions = None
        self.regElementList = None
        self.unpClusterList = None
        self.regElementDict = None
        self.flat_df_count = None
        self.non_clus_reg_hits = None
                
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
        
        self.regElementDict = {}
        
        TSS_EPDnew_ID = lines.pop(0).split()[1] #remove the first line from lines, and access its EPDnew ID
        iStart = 0

        for i in tqdm(range(len(lines))): #iterate through the lines in the file
        
            if lines[i].startswith('>'):
                
                seq = ''.join(lines[iStart:i]) #pull out the promoter sequence corresponding to TSS_EPDnew_ID
                
                self.regElementDict[TSS_EPDnew_ID] = seq #Store the EPDnewID: promoter sequence key:value pair in a dictionary for later
                
                seq = self.__str_to_np_seq(seq)
                
                self.__pwm_scan_multifasta(self.PWM_Kdref, PWM_rc, seq, threshold, n_mer, cols, TSS_EPDnew_ID, length5Prime)
                
                TSS_EPDnew_ID = lines[i].split()[1]
                
                iStart = i+1

        seq = ''.join(lines[iStart:]) #pull out the final promoter sequence corresponding to TSS_EPDnew_ID (note: this could have problems if only want to read part of a file)
        
        self.regElementDict[TSS_EPDnew_ID] = seq #Store the EPDnewID: promoter sequence key:value pair in a dictionary for later
        
        seq = self.__str_to_np_seq(seq)
        
        self.__pwm_scan_multifasta(self.PWM_Kdref, PWM_rc, seq, threshold, n_mer, cols, TSS_EPDnew_ID, length5Prime)
    
    def create_ID_conversion_dict(self, filename):
        """
        Load in the conversion file into a dict, with key = ENSEMBL ID, value = EPDnewID

        Parameters
        ----------
        filename : txt
            file from EPDnew database with EPDnewID - ENESMBL equivalence.

        Returns
        -------
        None.
        
        Modifies
        --------
        self.EPDnewIDconversions

        """
        
        self.EPDnewIDconversions = {}
        
        with open(filename) as fh:
            
            for line in fh:
                
                EPDnewID, ENSEMBL = line.strip().split()
                
                self.EPDnewIDconversions[ENSEMBL.strip()] = EPDnewID.strip()
         
    def restrict_hits_annot(self, dfOfHits, annotFilename):
        """
        Takes the dataframe (reg_hits) with EPDnewIDs and creates a new dataframe that
        is restricted to promoters with ChIP-seq signal 1kB from the TSS. Based on
        the ENCODE dataset, currently using data from
        https://www.frontiersin.org/articles/10.3389/fnbeh.2017.00035/full

        Parameters
        ----------
        dfOfHits : dataframe
            input dataframe with EPDnewID-labeled TFBSs or clusters as entries (typically reg_hits).
        annotFilename : csv file
            File containing annotations e.g. Promoter(<=1kb), ENSEMBL gene ID, ChIP score.

        Returns
        ----------
        restricted_reg_hits
        
        """
        
        self.restricted_reg_hits = pd.DataFrame(columns=dfOfHits.columns)
        
        if self.EPDnewIDconversions == None:
            print("load the ID conversion file")
            return
        
        self.EPDnewIDchipValuesDict = {}
        
        df = pd.read_csv('/Users/transcend/Python_Stuff/python_scripts/cluster_scan/input_data/EGR1_ENCODE_chip_peaks_gene_annotations.csv', header=1)
        
        df = df[df['annotation']=="Promoter (<=1kb)"] #Restrict the df to promoters with a ChIP peak 1kb from the TSS
        
        df = df[df['EnsemblID'].notnull()] #Some of the EnsemblIDs were not there, so skip them for now
                        
        for i in tqdm(df.index): #Iterate through refined annotation df
            
            if df.loc[i]["EnsemblID"] in self.EPDnewIDconversions: #If promoter with ChIP-peak is in EPDnew database
                
                self.EPDnewIDchipValuesDict[self.EPDnewIDconversions[df.loc[i]["EnsemblID"]]] = df.loc[i]["Score"] #Save Key-EPDnewID, Value-ChIP Peak
                
        for i in tqdm(dfOfHits.index):
                 
            if dfOfHits.loc[i]["EPDnew_ID"] in self.EPDnewIDchipValuesDict:
                
                self.restricted_reg_hits = self.restricted_reg_hits.append(dfOfHits.loc[i])
                
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
        
    def misc_process_df(self, df, threshold=None, resetIndex=False):
        """
        This method just groups some of the dataframe processing code that one might
        want to use during a scanning process.
        
        Operations: 1) Reset index 
                    2) Restrict regulatory hits to a threshold
        
        Attributes:
        reset: True or False
        threshold: Float, the Kd/Kdref threshold (preserve Scores lower than) to limit the dataframe to
        
        Returns
        -------
        None.

        """
        if threshold is not None:
            df = df[df.Score<threshold]   
            
        if resetIndex:
            df = df.reset_index(drop=True)
        
        return df
        
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
                unpClusterList.append(RegEle(sitesTemp, affinitiesTemp, startPosTemp, strandTemp))
                
                sitesTemp = [] #Clear temporary variables
                affinitiesTemp = []
                startPosTemp = []
                strandTemp = []
        
        if (maxGap <= -1):
            self.ovlpClusterList = unpClusterList
        else:
            self.unpClusterList = unpClusterList
            
    def generate_reg_elements_clusters(self, regHitsDF, maxGap, siteLength=None):
        
        """
        This method compiles single binding site hits dataframe into 
        clusters when you're using extracted regulatory sequences passed as input
        through a file with multiple FASTA seqs (typical case is promoters from EPDnew db).
        
        regHitsDF: dataframe corresponding to regHits, with single binding site hits annotated
        with the EPDnew ID of the associated TSS 
        
        maxGap: the maximum gap distance between binding sites that is allowed, for
        binding sites in the same cluster (before truncating a cluster). 
        
        Note: Might need to reset index of the regHitsDF that's being fed in here'
        """
        if siteLength is None:
            siteLength = self.n_mer
        if siteLength is None:
            print('ERROR: siteLength is Nonetype, either load the PWM or input a value directly')
            return
        
        self.non_clus_reg_hits = pd.DataFrame(columns=regHitsDF.columns)
        
        regElementList = []
        rSitesTemp = [] #variables starting with 'r' are temporary variables for the AnnotRegEle object corresponding to a Promoter (or other r'egulatory element)
        rAffinitiesTemp = []
        rStartPosTemp = []
        rStrandTemp = []
        
        unpClusterList = []
        cSitesTemp = [] #variables starting with 'c' are temporary variables for the AnnotRegEle object corresponding to a Cluster
        cAffinitiesTemp = []
        cStartPosTemp = []
        cStrandTemp = []
        
        for i in tqdm(range(len(regHitsDF.index)-1)): #iterate over the dataframe of hits
            
            samePromoter = (regHitsDF.loc[i+1].EPDnew_ID == regHitsDF.loc[i].EPDnew_ID) #determine if the next hit is in the same promoter
            
            #if the 'i+1' site is from the same promoter, then add the 'i' site to the promoter object
            if (samePromoter):
                rSitesTemp.append(regHitsDF.loc[i].Sequence)
                rAffinitiesTemp.append(regHitsDF.loc[i].Score)
                rStartPosTemp.append(regHitsDF.loc[i].Start)
                rStrandTemp.append(regHitsDF.loc[i].Strand)  
                
            else:  #if we're on the last site for the promoter (incl. if it's the only site in a promoter)
                rSitesTemp.append(regHitsDF.loc[i].Sequence)
                rAffinitiesTemp.append(regHitsDF.loc[i].Score)
                rStartPosTemp.append(regHitsDF.loc[i].Start)
                rStrandTemp.append(regHitsDF.loc[i].Strand)
                rEPDnew_ID = regHitsDF.loc[i].EPDnew_ID #append the promoter's ID

                #Transform temporary variables into a AnnotRegEle object and store it as an element of the list regElementList
                regElementList.append(AnnotRegEle(rSitesTemp, rAffinitiesTemp, rStartPosTemp, rStrandTemp, siteLength, rEPDnew_ID))
                rSitesTemp = [] #Clear temporary variables
                rAffinitiesTemp = []
                rStartPosTemp = []
                rStrandTemp = []                

            #if the distance to the next site is less than the maximum gap distance allowed, and the i+1 site is also from the same promoter
            if (((regHitsDF.loc[i+1].Start - regHitsDF.loc[i].End -1) <= maxGap) and (samePromoter)):
                cSitesTemp.append(regHitsDF.loc[i].Sequence)
                cAffinitiesTemp.append(regHitsDF.loc[i].Score)
                cStartPosTemp.append(regHitsDF.loc[i].Start)
                cStrandTemp.append(regHitsDF.loc[i].Strand)
            
            elif (cSitesTemp) : #if we're on the last site of the cluster
                cSitesTemp.append(regHitsDF.loc[i].Sequence)
                cAffinitiesTemp.append(regHitsDF.loc[i].Score)
                cStartPosTemp.append(regHitsDF.loc[i].Start)
                cStrandTemp.append(regHitsDF.loc[i].Strand)
                cEPDnew_ID = regHitsDF.loc[i].EPDnew_ID #also append the cluster's promoter region
                
                #Transform temporary variables into a AnnotRegEle object and store it as an element of the list unpClusterList
                unpClusterList.append(AnnotRegEle(cSitesTemp, cAffinitiesTemp, cStartPosTemp, cStrandTemp, siteLength, cEPDnew_ID))
                
                cSitesTemp = [] #Clear temporary variables
                cAffinitiesTemp = []
                cStartPosTemp = []
                cStrandTemp = []
                
            else: #if not from same promoter OR not within maxGap, and not the last site of the cluster, then it's a single binding site
                self.non_clus_reg_hits = self.non_clus_reg_hits.append(regHitsDF.loc[i], ignore_index=True)
        
        #Need to deal with the last iteration
        i=i+1
        
        rSitesTemp.append(regHitsDF.loc[i].Sequence)  #this makes sense because already checked this index above, and cleared if not the same promoter
        rAffinitiesTemp.append(regHitsDF.loc[i].Score)
        rStartPosTemp.append(regHitsDF.loc[i].Start)
        rStrandTemp.append(regHitsDF.loc[i].Strand)
        rEPDnew_ID = regHitsDF.loc[i].EPDnew_ID #append the promoter's ID
        regElementList.append(AnnotRegEle(rSitesTemp, rAffinitiesTemp, rStartPosTemp, rStrandTemp, siteLength, rEPDnew_ID))
        
        if (cSitesTemp): #if cSitesTemp is not empty, that means the final site is part of the cluster (otherwise it would have truncated and cleared already). If empty, then no cluster.
            cSitesTemp.append(regHitsDF.loc[i].Sequence)
            cAffinitiesTemp.append(regHitsDF.loc[i].Score)
            cStartPosTemp.append(regHitsDF.loc[i].Start)
            cStrandTemp.append(regHitsDF.loc[i].Strand)
            cEPDnew_ID = regHitsDF.loc[i].EPDnew_ID
            unpClusterList.append(AnnotRegEle(cSitesTemp, cAffinitiesTemp, cStartPosTemp, cStrandTemp, siteLength, cEPDnew_ID))
        
        else:  #if cSitesTemp is empty it's because the last site is not part of the cluster
            self.non_clus_reg_hits = self.non_clus_reg_hits.append(regHitsDF.loc[i], ignore_index=True)
                
        self.regElementList = regElementList
        
        self.unpClusterList = unpClusterList
        
    def plot_site_spacing(self, clusterDfList):
        
        """
        This method plots spacing between two binding sites on the X-axis and 
        count on the Y-axis.
        
        clusterDfList: the input list of cluster objects
        """

        df = pd.DataFrame([vars(s) for s in clusterDfList]) #Generate a dataframe from a list of Class objects

        #Make a series of lists into a one-dimensional list
        flat_df = []
        for sublist in df['spacing']:
            for item in sublist:
                flat_df.append(item)     

        #Turn the flat list into a flat DF
        flat_df = pd.DataFrame(flat_df)
        
        #Count the occurences of values in a DataFrame
        self.flat_df_count = pd.DataFrame(flat_df[0].value_counts())
        
        #Sort the DF by value
        order = np.arange(-8,36).tolist() #this needs to be changed, to account for loss of overlaps (-8 first, just add lists) 
        
        for i in order:
            
            if i not in self.flat_df_count.index: #If a spacing is not present, set its count to zero
                
                self.flat_df_count.loc[i] = 0
        
        self.flat_df_count = self.flat_df_count.loc[order]
        
        #Generating a bar plot from a flat dataframe
        fig, ax = plt.subplots()
        self.flat_df_count[0].plot(ax=ax, kind='bar')
    
    def ss_reglist_calc_meanOcc(self, df, highAffOccThresh=None, conc=16, consensusKd=16, concO=0.6):
        """
        This method takes as input a dataframe of binding site hits, and modifies the dataframe in order to have a meanOccupancy column.

        Parameters
        ----------
        df : dataframe
            the dataframe of hits and properties.
            
        highAffOccThresh: float
            if this variable is specified, it means that we only want to look at high affinity sites

        Returns
        -------
        The modified dataframe.

        """
        
        if highAffOccThresh:
            df = self.misc_process_df(df, threshold=highAffOccThresh, resetIndex=True)
        
        df["meanOcc"]=""
        
        for p in tqdm(range(len(df))):
            
            mOcc = single_site_statmech_model(df.loc[p].Score, conc, consensusKd, concO)
            
            df.at[p,'meanOcc'] = mOcc
                        
        return df
        
    def cluster_reglist_calc_meanOcc(self, regObjList, lowAffinityOccThresh=None, conc=16, consensusKd=16, concO=0.6):
        """
        Method calculates the mean occupancy, and converts the input list of objects into a dataframe, with a mean occupancy column.

        Parameters
        ----------
        regObjList : List of regulatory objects
            Object is from class RegEle.
            
        lowAffinityOccThresh : float
            If specified, this means we only care about the mean occupancy coming from low-affinity sites,
            which have Kd values above the threshold ratio stored in this variable. In order to implement
            this, for independent binding sites, we skip the occupancy contribution if the site's Kd is 
            below this threshold. For overlapping clusters, we multiply each state's weighted rel. mult.
            by the number of low-affinity binding sites instead of the total number of binding sites, before
            dividing by the partition function.

        Returns
        -------
        Dataframe corresponding to the input list of objects, with an additional column for the mean occupancy.

        """
        df = pd.DataFrame(columns=pd.DataFrame([vars(regObjList[0])]).columns) #Make an empty dataframe with columns correesponding to properties of regObjList
        df["meanOcc"]="" #add a 'meanOcc' column for mean occupancy
        
        for p in tqdm(range(len(regObjList))):
            
            meanOcc = general_cluster_statmech_model(regObjList[p], lowAffinityOccThresh, conc, consensusKd, concO)  #calculate the occupancy of the i'th object
            df = df.append(pd.DataFrame([vars(regObjList[p])]),ignore_index=True) #append the i'th object to df
            df.loc[p].meanOcc = meanOcc
        return df
                    
    def plot_meanOcc(self, dfCluster, dfSingle, binz=100, spacing=None, lowAffThresh=None):
        """
        Method plots the meanOcc column of the input dataframes histograms on the same plot.

        """
        binz = np.linspace(0,4.5,num=100)
        #binz = np.linspace(0,(dfCluster['meanOcc'].max()+0.5),num=100)
        
        plt.hist(dfCluster['meanOcc'], label='Clusters', alpha=0.4, bins=binz)
        plt.hist(dfSingle['meanOcc'], label='Single sites', alpha=0.4, bins=binz)
        plt.yscale('log', nonposy='clip')
        plt.legend(loc='best')
        
        if spacing and not lowAffThresh:
            plt.title('Spacing threshold: '+str(spacing)+'bp max between sites')
        
        if spacing and lowAffThresh:
            plt.title('Spacing threshold: '+str(spacing)+'bp max between sites '+'| Affinity threshold: '+str(lowAffThresh))
        
        plt.ylabel('Frequency')
        plt.xlabel('Mean occupancy')
        
        if spacing and not lowAffThresh:
            fname='/Users/transcend/SynologyDrive/Python_Stuff/python_scripts/cluster_scan/plots/mean_occupancy_distributions/clusters_vs_singlesites_histogram_meanOcc_threshold30_maxGap'+str(spacing)+'.jpeg'
            plt.savefig(fname, dpi=1200, quality=95)
        if spacing and lowAffThresh:
            fname='/Users/transcend/SynologyDrive/Python_Stuff/python_scripts/cluster_scan/plots/mean_occupancy_distributions/clusters_vs_singlesites_histogram_meanOcc_threshold30_maxGap'+str(spacing)+'_affThresh'+str(lowAffThresh)+'.jpeg'
            plt.savefig(fname, dpi=1200, quality=95)
        plt.show()
        
    def parameter_scan_spacing_affinity(self, regHitsDf, spacingList, lowAffinityOccThreshList, siteLen=9, plot=False):
        """
        This method takes as input a dataframe of regulatory hits (binding sites identified
        in promoter sequences), and then turns them into clusters vs. non-clusters through
        the generate_reg_elements_clusters method. Following this, this method takes the
        output (resets index where appropriate), and builds a dataframe with an extra meanOcc column,
        using ss_reglist_calc_meanOcc and cluster_regList_calc_meanOcc  methods. The sum of this column
        is calculated, and stored in ss_totalMeanOcc or clus_totalMeanOcc dataframes. This is done
        for each spacing value in spacingList. If plot==True, then a frequency vs. occupancy histogram 
        is plotted for each spacing, using plot_meanOcc method. This method finishes by returning 
        both ss_totalMeanOcc and clus_totalMeanOcc dataframes, and generating 
        a plot of spacing vs. total occupancy score.

        Parameters
        ----------
        plot : Boolean
            If true then plot the histograms corresponding to the frequency vs ooccupancy.
            
        regHitsDF: dataframe of regulatory hits
        
        spacingList: List of spacing values to calculate total mean occupancy sscoress for,
                    and for which to plot
                    
        lowAffinityOccThreshList: list
                If specified, this means we only care about the mean occupancy coming from low-affinity sites,
                which have Kd values above the threshold ratio stored in this variable. In order to implement
                this, for independent binding sites, we skip the occupancy contribution if the site's Kd is 
                below this threshold (done for each threshold, and each spacing). For overlapping clusters, 
                we multiply each state's weighted rel. mult. by the number of low-affinity binding sites 
                instead of the total number of binding sites, before dividing by the partition function.
        
        Returns
        -------
        Dataframes of spacing and/or threshold vs. totalMeanOcc, for single sites and clusters.

        """
        
        ss_totalMeanOcc = pd.DataFrame(columns=['Spacing', 'KdThreshold', 'Total_mean_occ'])
        clus_totalMeanOcc = pd.DataFrame(columns=['Spacing', 'KdThreshold', 'Total_mean_occ'])
        
        for i in range(len(spacingList)): # need to make sure that we're not actually impacting unpClusterList or non_clus_reg_hits
            
            self.generate_reg_elements_clusters(regHitsDf, spacingList[i], siteLength=siteLen)

            for j in tqdm(range(len(lowAffinityOccThreshList))):
                                
                clusDF = self.cluster_reglist_calc_meanOcc(self.unpClusterList, lowAffinityOccThresh=lowAffinityOccThreshList[j], conc=16, consensusKd=16, concO=0.6)
                
                clusSumMeanOcc = clusDF['meanOcc'].sum()
                                
                ssDF = self.ss_reglist_calc_meanOcc(self.non_clus_reg_hits, highAffOccThresh=lowAffinityOccThreshList[j], conc=16, consensusKd=16, concO=0.6)
                
                ssSumMeanOcc = ssDF['meanOcc'].sum()
                
                ss_totalMeanOcc = ss_totalMeanOcc.append({'Spacing': spacingList[i], 'KdThreshold': lowAffinityOccThreshList[j], 'Total_mean_occ': ssSumMeanOcc}, ignore_index=True)
                clus_totalMeanOcc = clus_totalMeanOcc.append({'Spacing': spacingList[i], 'KdThreshold': lowAffinityOccThreshList[j], 'Total_mean_occ': clusSumMeanOcc}, ignore_index=True)
                
                if plot==True:
                    self.plot_meanOcc(clusDF, ssDF, binz=100, spacing=spacingList[i], lowAffThresh=lowAffinityOccThreshList[j])
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        x = ss_totalMeanOcc['Spacing']
        y = ss_totalMeanOcc['KdThreshold']
        z = ss_totalMeanOcc['Total_mean_occ']
        ax.plot_trisurf(x,y,z, linewidth=0.2,label='Single Sites', color='blue')
        x = clus_totalMeanOcc['Spacing']
        y = clus_totalMeanOcc['KdThreshold']
        z = clus_totalMeanOcc['Total_mean_occ']
        ax.plot_trisurf(x,y,z, linewidth=0.2, label='Clusters', color='red')
        ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1, 1, 1.5, 1]))
        plt.xlabel('Spacing Threshold (bp)')
        plt.ylabel('Affinity Threshold (Kd/Kd$_(cons)$)')
        plt.show()
        
        return ss_totalMeanOcc, clus_totalMeanOcc
    
    def parameter_scan_spacing(self, regHitsDf, spacingList, siteLen=9, plot=False):
        """
        This method takes as input a dataframe of regulatory hits (binding sites identified
        in promoter sequences), and then turns them into clusters vs. non-clusters through
        the generate_reg_elements_clusters method. Following this, this method takes the
        output (resets index where appropriate), and builds a dataframe with an extra meanOcc column,
        using ss_reglist_calc_meanOcc and cluster_regList_calc_meanOcc  methods. The sum of this column
        is calculated, and stored in ss_totalMeanOcc or clus_totalMeanOcc dataframes. This is done
        for each spacing value in spacingList. If plot==True, then a frequency vs. occupancy histogram 
        is plotted for each spacing, using plot_meanOcc method. This method finishes by returning 
        both ss_totalMeanOcc and clus_totalMeanOcc dataframes, and generating 
        a plot of spacing vs. total occupancy score.

        Parameters
        ----------
        plot : Boolean
            If true then plot the histograms corresponding to the frequency vs ooccupancy.
            
        regHitsDF: dataframe of regulatory hits
        
        spacingList: List of spacing values to calculate total mean occupancy sscoress for,
                    and for which to plot
                    
        lowAffinityOccThreshList: list
                If specified, this means we only care about the mean occupancy coming from low-affinity sites,
                which have Kd values above the threshold ratio stored in this variable. In order to implement
                this, for independent binding sites, we skip the occupancy contribution if the site's Kd is 
                below this threshold (done for each threshold, and each spacing). For overlapping clusters, 
                we multiply each state's weighted rel. mult. by the number of low-affinity binding sites 
                instead of the total number of binding sites, before dividing by the partition function.
        
        Returns
        -------
        Dataframes of spacing and/or threshold vs. totalMeanOcc, for single sites and clusters.

        """
        
        ss_totalMeanOcc = pd.DataFrame(columns=['Spacing', 'Total_mean_occ'])
        clus_totalMeanOcc = pd.DataFrame(columns=['Spacing', 'Total_mean_occ'])
        
        for i in range(len(spacingList)): # need to make sure that we're not actually impacting unpClusterList or non_clus_reg_hits
            
            self.generate_reg_elements_clusters(regHitsDf, spacingList[i], siteLength=siteLen)
                                
            clusDF = self.cluster_reglist_calc_meanOcc(self.unpClusterList, conc=16, consensusKd=16, concO=0.6)
            
            clusSumMeanOcc = clusDF['meanOcc'].sum()
                            
            ssDF = self.ss_reglist_calc_meanOcc(self.non_clus_reg_hits, conc=16, consensusKd=16, concO=0.6)
            
            ssSumMeanOcc = ssDF['meanOcc'].sum()
            
            ss_totalMeanOcc = ss_totalMeanOcc.append({'Spacing': spacingList[i], 'Total_mean_occ': ssSumMeanOcc}, ignore_index=True)
            clus_totalMeanOcc = clus_totalMeanOcc.append({'Spacing': spacingList[i], 'Total_mean_occ': clusSumMeanOcc}, ignore_index=True)
            
            if plot==True:
                self.plot_meanOcc(clusDF, ssDF, binz=100, spacing=spacingList[i])
        
        plt.plot(ss_totalMeanOcc['Spacing'], ss_totalMeanOcc['Total_mean_occ'], label='Single sites', color='black', marker='o', linestyle='dashed',linewidth=2, markersize=5, alpha=0.7)
        plt.plot(clus_totalMeanOcc['Spacing'], clus_totalMeanOcc['Total_mean_occ'], label='Clusters', color='red', marker='o', linestyle='dashed',linewidth=2, markersize=5, alpha=0.7)
        plt.legend(loc='best')
        plt.ylabel('Total Mean Occ.  ')
        plt.xlabel('Between-site Spacing Threshold (bp)')
        fname='/Users/transcend/SynologyDrive/Python_Stuff/python_scripts/cluster_scan/plots/mean_occupancy_distributions/total_mean_occ_for_different_spacing_thresholds.jpeg'
        plt.savefig(fname, dpi=1200, quality=95)
        plt.show()   
        
        return ss_totalMeanOcc, clus_totalMeanOcc
    
    def ovlpClus_totalMeanOcc_affinityScan(self, regHitsDf, affThreshList, siteLen=9):
        """
        This method shows occupancy to overlapping binding sites vs. all other sites, as we vary the affinity threshold.
        How it works: takes regHitsDf, runs it through misc-process-df for each value
        in an affinityThreshList, and then calls generate_reg_elements_clusters 
        (spacing threshold=-1), then it  generates a meanOcc column for both 
        non_clus_reg_hits and unpClusterList, by calling ss_reglist_calc_meanOcc 
        and clus_reglist_calc_meanOcc. Then it sums the column and store the totals 
        in  a  meanOcc data frame beside  an affinityThreshold value. Then it  
        plots this data frame. This is valuable because it shows sensitivity of OVLP 
        clusters to affinity thresholds, to bypasss the argument that this is  just a product 
        of the primary binding site structure. Other approaches like restricting to occupancy
        contribution from  sites  below an affinity threshold in ovlpClusterss would not make sense  here. 

        Parameters
        ----------
        regHitsDf : TYPE
            DESCRIPTION.
        affThreshList : TYPE
            List of affinities to restric the regHitsDf to, before compiling  clusters.
        plot : TYPE, optional
            Plot if true.

        Returns
        -------
        ovlpClus_TotalMeanOcc_AffScan.

        """
        
        nonOvlp_totalMeanOcc = pd.DataFrame(columns=['affThresh', 'Total_mean_occ'])
        ovlp_totalMeanOcc = pd.DataFrame(columns=['affThresh', 'Total_mean_occ'])
        
        for i in range(len(affThreshList)):
            
            resRegHitsDf = self.misc_process_df(regHitsDf, threshold=affThreshList[i], resetIndex=True)
            
            self.generate_reg_elements_clusters(resRegHitsDf, maxGap=-1, siteLength=siteLen)
            
            clusDF = self.cluster_reglist_calc_meanOcc(self.unpClusterList, conc=16, consensusKd=16, concO=0.6)
        
            clusSumMeanOcc = clusDF['meanOcc'].sum()
            
            ssDF = self.ss_reglist_calc_meanOcc(self.non_clus_reg_hits, conc=16, consensusKd=16, concO=0.6)
            
            ssSumMeanOcc = ssDF['meanOcc'].sum()
            
            nonOvlp_totalMeanOcc = nonOvlp_totalMeanOcc.append({'affThresh': affThreshList[i], 'Total_mean_occ': ssSumMeanOcc}, ignore_index=True)
            ovlp_totalMeanOcc = ovlp_totalMeanOcc.append({'affThresh': affThreshList[i], 'Total_mean_occ': clusSumMeanOcc}, ignore_index=True)

        
        plt.plot(nonOvlp_totalMeanOcc['affThresh'], nonOvlp_totalMeanOcc['Total_mean_occ'], label='Everything else', color='black', marker='o', linestyle='dashed',linewidth=2, markersize=5, alpha=0.7)
        plt.plot(ovlp_totalMeanOcc['affThresh'], ovlp_totalMeanOcc['Total_mean_occ'], label='Overlapping sites', color='red', marker='o', linestyle='dashed',linewidth=2, markersize=5, alpha=0.7)
        plt.legend(loc='best')
        plt.ylabel('Total Mean Occ.  ')
        plt.xlabel('Affinity Cut-off (Kd/Kd$_(cons)$)')
        fname='/Users/transcend/SynologyDrive/Python_Stuff/python_scripts/cluster_scan/plots/mean_occupancy_distributions/ovlpCluster_total_meanOcc_forDiffAffinityCutOffs.jpeg'
        plt.savefig(fname, dpi=1200, quality=95)
        plt.show()  
            
@dataclass    
class RegEle: #Regulatory Elements
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
class AnnotRegEle(RegEle): # AnnotRegEle is a child class of RegEle, with the additional attribute geneID
    geneID: str = 'UNKNOWN' #This default value should actually never get triggered


















