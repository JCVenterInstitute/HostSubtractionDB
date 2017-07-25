#!/bin/sh
# Run bowtie2.
# Given reads and a reference, generate given alignment file.
# CPU optimiation: Align R1, then subtract mates of mapped reads from R2, then align R2.

SCRIPTDIR=$1  # path to other scripts
READSFILE1=$2  # do include .fastq or .fastq.gz
READSFILE2=$3  # do include .fastq or .fastq.gz
INDEXNAME=$4   # do not include .bt2
BAMBASE=$5     # do not include .sam or .bam 
THREADS=$6     # usually 4

if [ $# != 5 ];   then
    echo "ERROR. WRONG NUMBER OF PARAMETERS"
    echo "Usage: $0 <scriptdir> <reads1> <reads2> <index> <threads>"
   exit 1
fi
echo READSFILE1 $READSFILE1
echo READSFILE2 $READSFILE2
echo INDEXNAME $INDEXNAME
echo THREADS $THREADS

FILTER=fastq-filter-by-name.pl
echo FILTER $FILTER
if [ ! -f "${SCRIPDIR}/${FILTER}" ]; then
    echo "ERROR. REQUIRED SCRIPT NOT FOUND"
    echo ${SCRIPTDIR}/${FILTER}
    exit 2
fi
    
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
echo "LOCATION FOR SCRIPTS ${DIR}"
SAMTOOLS=/usr/local/bin/samtools
echo SAMTOOLS $SAMTOOLS
SAM1="${BAMBASE}.R1.sam"
BAM1="${BAMBASE}.R1.bam"
SAM2="${BAMBASE}.R2.sam"
BAM2="${BAMBASE}.R2.bam"
TMP_R2_FILE="tmp.${READSFILE2}"
TMP_IDS_FILE="tmp.${READSFILE1}.mapped_ids"

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
ALIGNMENT="--sensitive-local"
FASTQ="-q  --phred33"

# We have dectected a bug in bowtie2.
# If we use the "--un" option with cutadapt files (which can have zero-length reads),
# then bowtie writes fastq files that are corrupt in subtle ways.
UNALIGNED="--no-unal"                # keep unaligned out of the sam file

echo "RUN ALIGNER R1"
CMD="${BOWTIE_ALIGN} ${UNALIGNED} -p ${THREADS} ${ALIGNMENT} ${FASTQ} -x ${INDEXNAME} -U ${READSFILE1} -S ${SAM1}"
runit

echo "CONVERT SAM TO BAM"
CMD="${SAMTOOLS} view -h -b -o ${BAM1} ${SAM1} "
runit
CMD="rm -v ${SAM1}"
runit

# Do not filter. For subtraction, we found that any filter is too restrictive.

echo "SAMTOOLS TO EXTRACT MAPPED R1 READ IDs..."
${SAMTOOLS} view ${BAM1} | cut -f 1 >  ${TMP_IDS_FILE}
echo -n $?; echo " exit status"

EXCLUDE="-v"  # option to exclude the named reads
echo "MAKE ${TMP_R2_FILE}"
${SCRIPTDIR}/${FILTER} ${EXCLUDE} ${TMP_IDS_FILE} < ${READSFILE2} > ${TMP_R2_FILE}
echo -n $?; echo " exit status"

echo "MAP REMAINING READS FROM R2"

echo "RUN ALIGNER R1"
CMD="${BOWTIE_ALIGN} ${UNALIGNED} -p ${THREADS} ${ALIGNMENT} ${FASTQ} -x ${INDEXNAME} -U ${TMP_R2_FILE} -S ${SAM2}"
runit

echo "CONVERT SAM TO BAM"
CMD="${SAMTOOLS} view -h -b -o ${BAM2} ${SAM2} "
runit
CMD="rm -v ${SAM2}"
runit

