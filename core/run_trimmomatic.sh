#!/bin/sh

# Trim the paired reads

# Trimmomatic expects adapter sequence in current directory.
# Choices include: TruSeq3-SE.fa TruSeq3-PE.fa NexteraPE-PE.fa
# TruSeq2 was for GA II. TruSeq3 is for HiSeq, NextSeq, MiSeq.
# The TruSeq3-PE-2.fa file contains TruSeq3-PE plus more. 

JPATH=/usr/local/packages/trimmomatic-0.35
APATH="adapters"
ADAPT="TruSeq3-PE-2.fa"
echo "LOOKING FOR PROGRAM JAR FILE AND ADAPTER FILE"
echo JPATH ${JPATH}
echo APATH ${APATH}
echo ADAPT ${ADAPT}

cp -v ${JPATH}/${APATH}/${ADAPT} .

PATTERN="*_R1_*.fastq"
echo "WORKING ON FILES LIKE"
echo PATTERN ${PATTERN}

ENCODING="-phred33"
ENCODING=""     # trimmomatic will determine this automatically
THREADS="-threads 4"

for FF in ${PATTERN} ;
do
    
    INFILE1=$FF
    INFILE2=`echo $INFILE1 | sed 's/_R1_/_R2_/'`
    OUTPAIR1=trim.pair.${INFILE1}
    OUTPAIR2=trim.pair.${INFILE2}
    OUTSING1=trim.sing.${INFILE1}
    OUTSING2=trim.sing.${INFILE2}

    date
    CMD="java -jar ${JPATH}/trimmomatic-0.35.jar PE ${ENCODING} ${INFILE1} ${INFILE2} ${OUTPAIR1} ${OUTSING1} ${OUTPAIR2} ${OUTSING2} ILLUMINACLIP:${ADAPT}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ${THREADS}"
    echo $CMD
    $CMD
    echo -n $?; echo " exit status"
    date
done
echo "DONE"
