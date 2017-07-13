import pwm_scan


if __name__ == '__main__':

    scan = pwm_scan.PWMScan()

    # Load binding sites to generate a position weight matrix
    scan.load_pwm('example_target_sites.txt')

    # Load genome sequence (fasta)
    scan.load_sequence('example_genome_sequence.fna')

    # Load annotation (csv)
    scan.load_annotation('example_genome_annotation.gtf')

    # Launch scan and output a .csv file
    scan.launch_scan('output.csv', thres=12)
