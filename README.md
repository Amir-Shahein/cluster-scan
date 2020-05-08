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
scan.load_Kd_Kdref_pwm('Zif268_AAA_pwm.txt',9)
scan.load_sequence('GRCh38.fna', 'FASTA', 0.01) 
scan.launch_scan(threshold=50)
scan.generate_clusters(35)
----------


FASTA and GTF are very common formats for genome sequence and genome annotations, respectively.

The FASTA format for genome sequence: <https://en.wikipedia.org/wiki/FASTA_format>

The gene transfer format (GTF) for genome annotation: <https://en.wikipedia.org/wiki/Gene_transfer_format>

## Dependency

`numpy`, `pandas`, `re`
