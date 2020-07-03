

# This project uses some material from Yu-Cheng Lin's Github project: https://github.com/linyc74/pwm_scan (project accessed Feb. 28th 2020)

# https://en.wikipedia.org/wiki/Position_weight_matrix#Creation for good information on PWM calculations

# PWM Scan
Position-weight-matrix (PWM) scan through a genome.

## Install

```
pip install pwm_scan
```

## Getting started

At least two things are needed to perform a position-weight-matrix (PWM) scan:

- A TEXT (.txt) file containing the collection of target sites, in order to build the PWM.
- A genome sequence file in FASTA (.fna) format.

In addition, if you have the genome annotation in GTF format (.gtf), genes adjacent to the PWM-detected sites will be included in the output file.

```

scan = PWMScan()

# Load target sites to generate a position weight matrix
scan.load_Kd_Kdref_pwm('example_target_sites.txt') 
                                          ## A - input: a bunch of target sites, output: PWM-probability of base (max 1), based on 
                                          ## the frequency of occurences of that base in the target sites. Note: probably want to use
                                          ## load_PWM_alternate for now, which I created... or my version for delta delta G 

# Load genome sequence
scan.load_sequence('example_genome_sequence.fna') 
                                                  ## A - input: either ACGTs, or FASTA file. 
                                                  ## output: sequence of "ACGTs" but encoded as 0123

# Load annotation
scan.load_annotation('example_genome_annotation.gtf') 
                                                        ## A - skip annotations for now

# Launch scan and output the result in a .csv file
scan.launch_scan(filename='output.csv', threshold=12) 
                                                      ## A - convert PWM into PSM (using background GC content to calculate
                                                      ## the relative likelihood)
                                                      ## Note --  I removed a period after the 10 in argument threshold
```
------- FOR THE Restrict-to-annot stuff:

#Load the conversion dict
scan.create_ID_conversion_dict("/Users/transcend/Python_Stuff/python_scripts/cluster_scan/input_data/convert_EPDnewID_ENSEMBLgeneID.txt")



-----------
Note that the Zif268 mouse and HUMAN PWM is similar. I believe this applies across vertebrates, so it makes sense to scan it against a human genome.
#Typical Run for Humans with Zif268 (affinity matrix from MITOMI): 
runfile('/Users/transcend/Python_Stuff/python_scripts/cluster_scan/pwm_scan/main.py', wdir='/Users/transcend/Python_Stuff/python_scripts/cluster_scan/pwm_scan')
scan = PWMScan()
scan.load_Kd_Kdref_pwm('/Users/transcend/Python_Stuff/python_scripts/cluster_scan/input_data/Zif268_AAA_pwm.txt',9)

scan.load_sequence('GRCh38.fna', 'FASTA', 0.01) 
scan.launch_scan(threshold=50)
scan.generate_clusters(35)
scan.generate_overlapping_clusters()


#alternatively, skip load_sequence, and use this for multifasta/regulatory seq
scan.launch_scan_multifasta('/Users/transcend/Python_Stuff/python_scripts/cluster_scan/input_data/hg38_EPDnewSeqExtract.fa')
scan.generate_reg_elements_clusters(scan.reg_hits, 35)


#Restrict reg_hits to promoters with ChIP peaks
scan.restrict_hits_annot(reg_hits, "/Users/transcend/Python_Stuff/python_scripts/cluster_scan/input_data/EGR1_ENCODE_chip_peaks_gene_annotations.csv")

----------

To view the list of Clusters objects as a dataframe: pd.DataFrame([vars(s) for s in scan.ovlpClusterList])

----------
>Restrict reg_hits to the ChIP peaks
>restrict the df to a lower threshold
>use restricted_reg_hits to generate clusters/regulatory elements
>unpack the cluster list into a df
>feed this new df into the plot method
----------


FASTA and GTF are very common formats for genome sequence and genome annotations, respectively.

The FASTA format for genome sequence: <https://en.wikipedia.org/wiki/FASTA_format>

The gene transfer format (GTF) for genome annotation: <https://en.wikipedia.org/wiki/Gene_transfer_format>

## Dependency

`numpy`, `pandas`, `re`

---------------------RANDOM DOCUMENTATION ON THE APPROACH---------------------
Note: The original scan took the following approach to annotations: Use the major 
.fna file for the organism (e.g. https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39),
then find CDS start sites (with the corresponding gene name/ID), identify the 
position range corresponding to -500~bp (or however long) upstream of this TSS 
(an idealistic promoter region), and see if anything in self.hits is within that position range. 
One problem I envision with this approach, is that some .fna files have multiple 
labels (i.e. not just CDS; exon, etc.), and I was not sure which of these have 
upstream regions that correspond to promoters. Also, this approach is a bit rough.
==>An approach that is seemingly more straightforward, that I'm currently using, 
relies on the EPD (eukaryotic promoter database, it's actually curated by the 
SIB + EPFL). This database identifies transcription start sites (TSSs) in differet organisms, which have been
validated experimentally (the data is actually maintained for the purpose of 
identifying promoters), and through the website you can obtain whatever length 
of sequence before and after the TSS in a BED file (I believe this is a state-of-the-art 
'high-throughput' way to identify promoters) -- https://epd.epfl.ch/get_promoters.ph
==>I can then simply run the scan on these extracted sequences (to generate self.hits
from the beginning), save an extra column with the gene name, an extra column with 
the distance to the TSS, and then collate self.hits into clusters, preserving both 
of these columns. 
==>A benefit of this approach is also that there is dramatically less sequence to 
scan through -- e.g. 18 million bases for the human genome. However, it still makes sense 
to do both this approach and scan the entire organism's genome, because promoters aren't
everything... there is action-at-a-distance in gene regulation (we're not even extracting
the enhancers). Unfortunately, for the time being I'm unaware of any solid 
enhancer databases (I don't think this is well known), but this can be incorporated at a later date. 


Note: Currently I'm using EPDnew sequence extraction, and they don't provide much info. I'm trusting
that what they return is -1000 and +100 of the TSS when looking at it 5' to 3' (this needs to be
confirmed). However, they do not provide whether the sequence is the + or - strand in this file. For now
in these DataFrames, I can include strand information as C (coding) and N (noncoding), and in the future, if I want to 
label it as being the + or - strand objectively, then I would need to match the EPDnew ID against FPS file
(https://epd.epfl.ch/get_promoters.php) and retrieve the + or - from there. 


In order to build the refined annotated list of promoters, I used this paper: https://www.frontiersin.org/articles/10.3389/fnbeh.2017.00035/full,
who did an analysis of EGR-1 (Zif268) ChIP-seq peaks (from ENCODE), and identified genes nearby: 
# Indeed, as the ENCODE project included EGR1 as part of the tier 1 chromatin immunoprecipitation followed by deep sequencing 
# (ChIP-seq), a wealth of information regarding EGR1 DNA binding characteristics and target genes has been made available 
# (ENCODE Project Consortium, 2012). In particular, we are thus able to analyze and compare the binding pattern of 161 
# transcription factors across 91 cell types and a total of 4,380,444 genomic regions, among which 44,985 correspond to 
# an EGR1 binding event. Out of the 15,872 genes thus annotated, 8552 (53.9%) contain at least one EGR1 binding region 
# (peak) within 3 kb of their transcription start site (TSS), which indicates that across several human cell types, EGR1 
# can bind a very large number of genes and thus potentially regulate a very large gene expression profile (see full 
# annotated list in Supplementary Table S1). As previously reported, EGR1 binds in close vicinity to the TSS (Project 
# Kubosaki et al., 2009; ENCODE Project Consortium, 2012), but even though 41.6% of all EGR1 peaks are located within 
# the promoter region, 26.4% are located within intronic regions.
See csv file: EGR1_ENCODE_chip_peaks_gene_annotations.csv

Note: file convert_EPDnewID_ENSEMBLgeneID.txt contains conversion from EPDnewID to ENSEMBL, which is provided in the excel file

