#!/bin/sh

# Assuming...
# + reads in *.fastq.gz (ok if R1 and R2 are in separate files)
# + list of read IDs (assume both reads of a pair have the same name)

# We will...
# + extract from fastq the unmapped reads

# Dependency...
# + unix
# + cut
# + perl
# + gzip
# + ${SDB_HOME} environment variable
# + ${SDB_HOME}/util/fastq-filter-by-name.pl

BASEDIR=$(dirname "$0")
SDB_HOME="${BASEDIR}/.."
echo SDB_HOME $SDB_HOME
UTILITY=${SDB_HOME}/util/fastq-filter-by-name.pl 
echo UTILITY $UTILITY

date
READS_IN=$1
NAMES_IN=$2
READS_OUT=$3
echo READS_IN $READS_IN
echo READS_OUT $READS_OUT
ls -l $READS_IN
echo NAMES_IN $NAMES_IN
ls -l $NAMES_IN

date
echo "Get reads..."
EXCLUDE="-v"  # option to exclude the named reads
gunzip -c ${READS_IN} | ${UTILITY} ${EXCLUDE} ${NAMES_IN} > ${READS_OUT}
echo -n $?; echo " exit status"
ls -l ${READS_OUT}
date
