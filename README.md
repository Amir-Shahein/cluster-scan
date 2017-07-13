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
scan.launch_scan('output.csv', thres=12)
```

## Dependency

`numpy`, `pandas`, `re`
