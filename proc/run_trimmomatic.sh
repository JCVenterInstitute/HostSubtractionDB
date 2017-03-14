#!/bin/sh

JPATH=/usr/local/packages/trimmomatic-0.35

#/usr/local/packages/trimmomatic-0.35/adapters/TruSeq3-SE.fa
#/usr/local/packages/trimmomatic-0.35/adapters/TruSeq3-PE.fa
#/usr/local/packages/trimmomatic-0.35/adapters/NexteraPE-PE.fa
cp /usr/local/packages/trimmomatic-0.35/adapters/TruSeq3-SE.fa .
																    
for FF in *.fastq;
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
