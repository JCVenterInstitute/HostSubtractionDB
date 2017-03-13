#!/bin/sh

# Assuming...
# + reads in *.fastq.gz (ok if R1 and R2 are in separate files)
# + mapping in *.bam (lists reads mapped to SDB; sam format ok)

# We will...
# + extract from bam files a list of mapped read IDs
# + extract from fastq the unmapped reads

# Dependency...
# + unix
# + cut
# + perl
# + samtools
# + gzip
# + ${SDB_HOME} environment variable
# + ${SDB_HOME}/util/fastq-filter-by-name.pl

BASEDIR=$(dirname "$0")
SDB_HOME="${BASEDIR}/.."
echo SDB_HOME $SDB_HOME
UTILITY=${SDB_HOME}/util/fastq-filter-by-name.pl 
echo UTILITY $UTILITY

SAMTOOLS=/usr/local/bin/samtools
echo SAMTOOLS $SAMTOOLS

date
READS_IN=$1
MAPS_IN=$2
READS_OUT=$1.SDB_filtered.fastq
echo READS_IN $READS_IN
echo READS_OUT $READS_OUT
ls -l $READS_IN
echo MAPS_IN $MAPS_IN
ls -l $MAPS_IN
TEMPFILE=${MAPS_IN}.read_names
echo TEMPFILE $TEMPFILE

date
echo "Get read names..."
${SAMTOOLS} view ${MAPS_IN} | cut -f 1 > ${TEMPFILE}
echo -n $?; echo " exit status"
ls -l ${TEMPFILE}
date
echo "Get reads..."
EXCLUDE="-v"  # option to exclude the named reads
gunzip -c ${READS_IN} | ${UTILITY} ${EXCLUDE} ${TEMPFILE} > ${READS_OUT}
echo -n $?; echo " exit status"
ls -l ${READS_OUT}
date
