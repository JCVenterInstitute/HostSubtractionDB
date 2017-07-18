#!/bin/sh
# Submit grid jobs for each fastq file in current directory.

# Assuming...
# + reads in *_R1_*.fastq and *_R2_*.fastq
# + mappings in *_R1_*.bam and *_R2_*.bam

# We will...
# + produce new unmapped.*_R1_*.fastq and unmapped.*_R2_*.fastq

if [ -z "$1" ]; then
    echo ERROR MISSING PARAMETER
    echo "Usage: $0 <output_prefix>"
    exit 1
fi
HERE=`pwd`
echo CURRENT DIRECTORY $HERE
BASE=`basename "$0"`
echo "THIS PROGRAM IS ${BASE}"
DIR=`dirname "$0"`
if [ "${DIR}" == "." ]; then
    echo "ERROR SCRIPT MUST BE INVOKED WITH ABSOLUTE PATH"
    exit 2
fi
echo LOCATION FOR SCRIPTS $DIR
SCRIPT=subtract_mapped_reads.sh
echo SCRIPT TO BE RUN $SCRIPT
if [ -z "$GRIDCMD" ]; then
    echo "GRIDCMD NOT SET, USING DEFAULT"
    # memory requirement depends on fastq size. 2G is common.
    QSUB="qsub -cwd -b n -A DHSSDB -P 8370 -N human1 -l medium -l memory=4g -j y"
    QSUB="qsub -cwd -b n -A DHSSDB -P 8370 -N human1 -l memory=4g -j y"
else
    echo "USING GRIDCMD ENVIRONMENT VARIABLE"
    QSUB=${GRIDCMD}
fi
echo "BASE GRID COMMAND ${QSUB}"

OUT_PREFIX=$1  # examples: unmapped or non-human
BASE_PATTERN="*_R1_*"
READS_SUFFIX="fastq"
MAP_SUFFIX="bam"

echo "GET READS FROM ${BASE_PATTERN}.${READS_SUFFIX}"
echo "GET MAPS FROM ${BASE_PATTERN}.${MAP_SUFFIX}"
echo "WRITE READS TO ${OUT_PREFIX}.${BASE_PATTERN}.${READS_SUFFIX}"

function runit () {
    echo "COMMAND"
    echo ${CMD}
    date
    nice ${CMD}
    echo -n $?;    echo " exit status"
    date
}

for FF in ${BASE_PATTERN}.${READS_SUFFIX} ; do
    R1=${FF}
    R2=`echo ${R1} | sed 's/_R1_/_R2_/'`
    BAM1=`echo ${R1} | sed "s/.${READS_SUFFIX}/.${MAP_SUFFIX}/"`
    BAM2=`echo ${R2} | sed "s/.${READS_SUFFIX}/.${MAP_SUFFIX}/"`
    echo R1 $R1 R2 $R2
    echo BAM1 $BAM BAM2 $BAM2
    CMD="${QSUB} -o ${HERE} ${DIR}/${SCRIPT} ${DIR} ${R1} ${R2} ${BAM1} ${BAM2} ${OUT_PREFIX}"
    echo "GRID SUBMIT"
    echo $CMD
    runit
done

echo "DONE"
