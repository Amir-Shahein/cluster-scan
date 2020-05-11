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


Note: alternatively, skip load_sequence, and use 
scan.launch_scan_multifasta('/Users/transcend/Python_Stuff/python_scripts/cluster_scan/input_data/hg38_EPDnewSeqExtract.fa')
----------

To view the list of Clusters objects as a dataframe: pd.DataFrame([vars(s) for s in scan.ovlpClusterList])

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