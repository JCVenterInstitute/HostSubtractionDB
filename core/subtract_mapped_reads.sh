#!/bin/sh

if [ $# != 4 ]; then echo "Usage: $0 <R1> <R2> <BAM> <out>"; exit 1; fi
IN_R1=$1
IN_R2=$2
IN_BAM=$3
OUT_PREFIX=$4
OUT_R1=${OUT_PREFIX}.${IN_R1}
OUT_R2=${OUT_PREFIX}.${IN_R2}

echo "GET READS FROM ${IN_R1} AND ${IN_R2}"
echo "GET MAPS FROM ${IN_BAM}"
echo "WRITE READS TO ${OUT_R1} AND ${OUT_R2}"

HERE=`pwd`
echo CURRENT DIRECTORY $HERE
BASE=`basename "$0"`
echo THIS PROGRAM IS $BASE
DIR=`dirname "$0"`
if [ "${DIR}" == "." ]; then
    echo "ERROR SCRIPT MUST BE INVOKED WITH ABSOLUTE PATH"
    exit 2
fi
echo LOCATION FOR SCRIPTS $DIR
SAMTOOLS=/usr/local/bin/samtools
FILTER=fastq-filter-by-name.pl
echo SAMTOOLS $SAMTOOLS
echo FILTER $FILTER

function runit () {
    echo "COMMAND"
    echo "${CMD} > $OUTFILE"
    date
    ${CMD} > ${OUTFILE}
    echo -n $?;    echo " exit status"
    date
}

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

# For subtraction, we found that any filter is too restrictive.

TMP_IDS_FILE="${OUT_PREFIX}.${IN_BAM}.ids"

CMD="${SAMTOOLS} view ${IN_BAM} | cut -f 1 "
OUTFILE="${OUT_PREFIX}.${IN_BAM}.ids"
# runit
# The command above complains the bam must be indexed, but only if run by expanding CMD.
# The invocation below is a work-around.
echo "SAMTOOLS TO EXTRACT MAPPED READ IDs..."
${SAMTOOLS} view ${IN_BAM} | cut -f 1 > ${OUT_PREFIX}.${IN_BAM}.ids

EXCLUDE="-v"  # option to exclude the named reads
CMD="cat ${IN_R1} | ${DIR}/${FILTER} ${EXCLUDE} ${TMP_IDS_FILE} "
OUTFILE="${OUT_R1}"
runit
CMD="cat ${IN_R2} | ${DIR}/${FILTER} ${EXCLUDE} ${TMP_IDS_FILE} "
OUTFILE="${OUT_R2}"
runit
echo DONE
