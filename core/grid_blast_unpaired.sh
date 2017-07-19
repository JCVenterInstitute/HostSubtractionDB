#!/bin/sh
# Submit blast job to grid for each fastq in current directory.
# Blast requires FASTA but reads are always in FASTQ.
# For convenience, we accept FASTQ but write a FASTA version.
# A post-process could clean up the FASTA files.

if [ $# != 1 ]; then
    echo "Usage: $0 <suffix>"
    exit 1
fi
READS_SUFFIX=$1    # use fq or fastq but omit the period
echo READS_SUFFIX $READS_SUFFIX

HERE=`pwd`
echo "CURRENT DIRECTORY ${HERE}"
BASE=`basename "$0"`
echo "THIS PROGRAM IS ${BASE}"
DIR=`dirname "$0"`
if [ "${DIR}" == "." ]; then
    echo "ERROR SCRIPT MUST BE INVOKED WITH ABSOLUTE PATH"
    exit 2
fi
echo LOCATION FOR SCRIPTS $DIR
SCRIPT=run_blast.sh
echo SCRIPT TO BE RUN $SCRIPT
THREADS=4
echo THREADS $THREADS
if [ -z "$GRIDCMD" ]; then
    echo "GRIDCMD NOT SET, USING DEFAULT"
    # BLAST requires ~16GB RAM. We reserve 16 GB by 4 threads * 4g each.
    QSUB="qsub -cwd -b n -A DHSSDB -P 8370 -N human1 -l memory=4g -j y"
    QSUB="qsub -cwd -b n -A DHSSDB -P 8370 -N human1 -pe threaded ${THREADS} -l medium -l memory=4g -j y"
else
    echo "USING GRIDCMD ENVIRONMENT VARIABLE"
    QSUB=${GRIDCMD}
fi
echo "BASE GRID COMMAND ${QSUB}"

function runit () {
    echo "COMMAND"
    echo ${CMD}
    date
    nice ${CMD}
    echo -n $?;    echo " exit status"
    date
}

for FF in *.${READS_SUFFIX} ; do
    READS=${FF}
    echo READS $READS
    FASTA=`echo $READS | sed "s/.$READS_SUFFIX/.fasta/"`
    awk '{if (++C ==5)C=1; if (C==1) print ">" substr($1,2); if (C==2) print $0;}' $READS > $FASTA
    CMD="${QSUB} -o ${HERE} ${DIR}/${SCRIPT} ${FASTA}"
    echo "GRID SUBMIT"
    runit
done
echo DONE

    
