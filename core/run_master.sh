#!/bin/sh

# run_master.sh
# Master program runs the entire subtraction.

# TO DO
# * clean up any *.sam and *.bam from previous iteration
# * clean up any *.fastq.gz left from previous run
# * use same script for all run_bowtie
# * parameterize the input filenames and number of files
# * parameterize the output filenames

# * Call subtract_mapped_reads.sh from here instead of within run_bowtie

# Order of operations subject to change

run_bowtie_align_to_sdb.sh
run_trimmomatic.sh
run_bowtie_align_to_virus.sh
run_bowtie_align_to_univec.sh
run_bowtie_align_to_mycoplasma.sh

