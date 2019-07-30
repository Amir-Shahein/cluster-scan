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
scan.load_pwm('example_target_sites.txt')

# Load genome sequence
scan.load_sequence('example_genome_sequence.fna')

# Load annotation
scan.load_annotation('example_genome_annotation.gtf')

# Launch scan and output the result in a .csv file
scan.launch_scan(filename='output.csv', threshold=12)
```

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
