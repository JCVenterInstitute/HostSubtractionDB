#!/bin/sh
# Submit grid jobs for each fastq file in current directory.

FILESUFFIX=$1 # use fastq or fq or fq.gz but omit leading period
INDEXNAME=$2  # do not include .bt2

if [ $# != 2 ]
then
    echo "Usage: $0 <reads_suffix> <index_prefix>"
    exit 1
fi
echo FILESUFFIX $FILESUFFIX
echo INDEXNAME $INDEXNAME

# Special case handling required for timmomatic outputs:
# Process the R1 and R2 that are still paired.
# Process the R1 and R2 that are no longer paired.

HERE=`pwd`
echo CURRENT DIRECTORY $HERE
BASE=`basename "$0"`
echo THIS PROGRAM IS $BASE
DIR=`dirname "$0"`
echo LOCATION FOR SCRIPTS $DIR
SCRIPT=run_bowtie_align_unpaired.sh
echo SCRIPT TO BE RUN $SCRIPT
THREADS=4
echo THREADS $THREADS
if [ -z "$GRIDCMD" ]; then
    echo "GRIDCMD NOT SET, USING DEFAULT"
    QSUB="qsub -cwd -b n -A DHSSDB -P 8370 -N human1 -pe threaded ${THREADS} -l medium -l memory=1g -j y"
else
    echo "USING GRIDCMD ENVIRONMENT VARIABLE"
    QSUB=${GRIDCMD}
fi
echo BASE GRID COMMAND $QSUB

function runit () {
    echo "COMMAND"
    echo ${CMD}
    date
    nice ${CMD}
    echo -n $?;    echo " exit status"
    date
}

ls -l *.${FILESUFFIX}
for FASTQ in *.${FILESUFFIX} ; do
    echo FASTQ $FASTQ
    PREFIX=`echo ${FASTQ} | sed "s/.${FILESUFFIX}//" `
    echo PREFIX $PREFIX
    CMD="${QSUB} -o $HERE ${DIR}/${SCRIPT} ${FASTQ} ${INDEXNAME} ${PREFIX}"
    echo "GRID SUBMIT"
    echo $CMD
    runit
done

echo "DONE"
