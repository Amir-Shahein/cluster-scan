import pwm_scan


if __name__ == '__main__':

    scan = pwm_scan.PWMScan()

    # Load binding sites to generate a position weight matrix
    scan.load_pwm(filename='example_target_sites.txt')

    # Load genome sequence (fasta)
    scan.load_sequence(seq='example_genome_sequence.fna')

    # Load annotation (csv)
    scan.load_annotation(filename='example_genome_annotation.gtf')

    # Launch scan and output a .csv file
    scan.launch_scan(filename='output.csv', threshold=12)
