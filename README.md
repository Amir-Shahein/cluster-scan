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
import pwm_scan

scan = pwm_scan.PWMScan()

# Load target sites to generate a position weight matrix
scan.load_pwm('example_target_sites.txt') ## A - input: a bunch of target sites, output: PWM-probability of base (max 1), based on 
                                          ## the frequency of occurences of that base in the target sites

# Load genome sequence
scan.load_sequence('example_genome_sequence.fna') ## A - input: either ACGTs, or FASTA file. 
                                                  ## output: sequence of "ACGTs" but encoded as 0123

# Load annotation
scan.load_annotation('example_genome_annotation.gtf') ## A - skip annotations for now

# Launch scan and output the result in a .csv file
scan.launch_scan(filename='output.csv', threshold=12) ## A - convert PWM into PSM (using background GC content to calculate
                                                      ## the relative likelihood)
                                                      ## Note --  I removed a period after the 10 in argument threshold
```

# import pwm_scan
# scan = pwm_scan.PWMScan()
# scan.load_pwm_alternate('MA01623.txt', 11)

## File format

The format of the input TEXT (.txt) file of target sites is very flexible. It can be comma, space or tab delimited. For example, the following three formats are all acceptable. Also note that all target sites should have the same length.

```
TTGATTCCGGTCAA,TTGACTTTCATCAA,TTGATTGCCATCAA,TTGACCGGAATCAA,TTGACGGCCGTCAA
```

```
TTGATTCCGGTCAA TTGACTTTCATCAA TTGATTGCCATCAA TTGACCGGAATCAA TTGACGGCCGTCAA
```

```
TTGATTCCGGTCAA <tab> TTGACTTTCATCAA <tab> TTGATTGCCATCAA <tab> TTGACCGGAATCAA <tab> TTGACGGCCGTCAA
```

FASTA and GTF are very common formats for genome sequence and genome annotations, respectively.

The FASTA format for genome sequence: <https://en.wikipedia.org/wiki/FASTA_format>

The gene transfer format (GTF) for genome annotation: <https://en.wikipedia.org/wiki/Gene_transfer_format>

## Dependency

`numpy`, `pandas`, `re`
