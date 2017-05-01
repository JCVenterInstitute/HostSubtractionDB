#!/bin/sh

JPATH=/usr/local/packages/trimmomatic-0.35
APATH="adapters"
ADAPT="TruSeq3-SE.fa"
echo "LOOKING FOR PROGRAM JAR FILE AND ADAPTER FILE"
echo JPATH ${JPATH}
echo APATH ${APATH}
echo ADAPT ${ADAPT}
# Trimmomatic expects adapter sequence in current directory.
# Choices include: TruSeq3-SE.fa TruSeq3-PE.fa NexteraPE-PE.fa

cp -v ${JPATH}/${APATH}/${ADAPT} .

PATTERN="*.fastq.gz"
echo "WORKING ON FILES LIKE"
echo PATTERN ${PATTERN}

for FF in *.fastq.gz;
do
    
    INFILE=$FF
    OUTFILE=trim.${FF}

    date
    CMD="java -jar ${JPATH}/trimmomatic-0.35.jar SE -phred33 $INFILE $OUTFILE ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
    echo $CMD
    $CMD
    echo -n $?; echo " exit status"
    date
done
echo "DONE"
