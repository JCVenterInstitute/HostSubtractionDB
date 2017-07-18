#!/bin/sh

FILTER=./fastq-filter-by-name.pl
R1=../NextSeq/nextseq.cutadapt.R1.fastq
R2=../NextSeq/nextseq.cutadapt.R2.fastq
OUT1=cutadapt.unmapped.R1.fastq
OUT2=cutadapt.unmapped.R2.fastq
SAMTOOLS=/usr/local/bin/samtools
SAMFILE=nextseq_to_canuArrow2.bam
MAPPED=cutadapt.mapped_pairs_mapq5.ids
echo INPUTS
ls -l ${R1}
ls -l ${R2}
ls -l ${MAPPED}

READER="gunzip -c"  # use this on zipped files
READER="cat"        # use this on fastq files

# This filter requires mapQ>=5 on both reads.
# For a bowtie2 mapping, hardly any reads pass.
# This assumes sorted inputs only writes the second.
#    awk '{if ($5>=5) {if ($1==p) print $1; p=$1;}}' \

# This filter requires the "each properly aligned" flag.
# In bowtie2 output, this matches intutiion.
# It requires both reads map to same transcript.
# It allows very few mismatches or indels.
# If pairs share an ID, this will write the ID twice
# (put the perl script will uniqify it by using a hash).
#    awk 'BEGIN {mask=02;} {if (and($2,mask)>0) print $1;}}' \

echo "FILTER BAM, WRITE IDS"
date
echo SAMTOOLS
${SAMTOOLS} view ${SAMFILE} |\
    awk 'BEGIN {mask=02;} {if (and($2,mask)>0) print $1;}' \
	> ${MAPPED}
date
echo "USE IDS TO FILTER FASTQ"
echo FILTER 1
${READER} ${R1} | ${FILTER} -v ${MAPPED} > cutadapt.unmapped.R1.fastq &
echo FILTER 2
${READER} ${R2} | ${FILTER} -v ${MAPPED} > cutadapt.unmapped.R2.fastq &
echo RUNNING
