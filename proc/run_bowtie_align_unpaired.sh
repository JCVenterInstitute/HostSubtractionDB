#!/bin/sh
# Run bowtie2.
# Align given reads.
# Align to given reference index.
# Generate given alignment file.

READSFILE=$1  # do include .fastq or .fastq.gz
INDEXNAME=$2  # do not include .bt2
SAM=$3        # do not include .sam or .bam

if [ $# != 3 ]
   then
   echo "Usage: 3 parameters"
   exit 1
fi
echo READSFILE $READSFILE
echo INDEXNAME $INDEXNAME
echo SAM $SAM

# This is the command to align reads using the index.
# The command has 'small' and 'large' binaries.
# Use small unless the target genome is > 4GB.
# Our target genome is about 3 GB.
export BOWTIE_ALIGN=/usr/local/bin/bowtie2
export BOWTIE_BUILD=/usr/local/bin/bowtie2-build
#export BOWTIE_ALIGN=/usr/local/bin/bowtie2-align-s
#export BOWTIE_ALIGN=/usr/local/bin/bowtie2-align-l
export BOWTIE_VERSION="/usr/local/bin/bowtie2 --version"
# Enable python, especially on the grid
source /usr/local/common/env/python.sh
SAMTOOLS=/usr/local/bin/samtools

echo BOWTIE VERSION
${BOWTIE_VERSION}

function runit () {
    echo "COMMAND"
    echo ${CMD}
    date
    nice ${CMD}
    echo -n $?;    echo " exit status"
    date
}

# Note bowtie2 can take fastq or fastq.gz. Assume all gz.
# If given a pair, bowtie will insist that both reads map and the pair is concordant.
# By default, bowtie uses "mixed mode", aligning separately if pair align fails.
# For maximum sensitivity, map reads separately, not as pairs.

THREADS="-p 4"
#THREADS="-p 4 -mm"
#ALIGNMENT="--end-to-end "
ALIGNMENT="--sensitive-local"
FASTQ="-q  --phred33"

# We have dectected a bug in bowtie2.
# If we use the "--un" option with cutadapt files (which can have zero-length reads),
# then bowtie writes fastq files that are corrupt in subtle ways.

UNALIGNED="--un R1.${MYBASE}.u.fq"   # warning: bug
UNALIGNED="--un R2.${MYBASE}.u.fq"
UNALIGNED="--no-unal"                # keep unaligned out of the sam file

echo "RUN ALIGNER"
CMD="${BOWTIE_ALIGN} ${UNALIGNED} ${THREADS} ${ALIGNMENT} ${FASTQ} -x ${INDEXNAME} -U ${READSFILE} -S ${SAM}.sam"
runit

echo "CONVERT SAM TO BAM"
CMD="${SAMTOOLS} view -h -b -o ${SAM}.bam ${SAM}.sam "
runit
rm ${SAM}.sam

echo "DONE"
