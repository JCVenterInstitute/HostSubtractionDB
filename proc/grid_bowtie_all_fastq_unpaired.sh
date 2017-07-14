#!/bin/sh
# Submit grid jobs for each fastq file in current directory.

FILESUFFIX=$1 # use fastq or fq or fq.gz but omit leading period
INDEXNAME=$2  # do not include .bt2

if [ $# != 2 ]
then
    echo "Usage: 2 parameters"
    exit 1
fi
echo FILESUFFIX $FILESUFFIX
echo INDEXNAME $INDEXNAME

# Special case handling required for timmomatic outputs:
# Process the R1 and R2 that are still paired.
# Process the R1 and R2 that are no longer paired.

HERE=`pwd`
BASE=`basename "$0"`
DIR=`dirname "$0"`
SCRIPT=run_bowtie_align_unpaired.sh
echo DIR $DIR
echo BASE $BASE
echo SCRIPT $SCRIPT

function runit () {
    echo "COMMAND"
    echo ${CMD}
    date
    nice ${CMD}
    echo -n $?;    echo " exit status"
    date
}

THREADS=4
echo WHICH QSUB
which qsub
QSUB="qsub -cwd -b n -A DHSSDB -P 8370 -N human1 -pe threaded ${THREADS} -l medium -l memory=1g -j y"

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
