#!/bin/sh

IN[0]=GCA_000001405.22_GRCh38.p7_genomic.fasta
IN[1]=Mycoplasma.remove_blanks.fasta
IN[2]=UniVec.fasta
IN[3]=phiX.fasta
IN[4]=clean.hepg2_unmapped_contigs.fasta
IN[5]=clean.huh7_unmapped_contigs.fasta
IN[6]=clean.jurkat_unmapped_contigs.fasta
IN[7]=all_Sendai.fasta
IN[8]=all_Zika.fasta
IN[9]=all_Newcastle.fasta
IN[10]=all_Ebola.fasta
IN[11]=all_Ecoli.fasta

NAME[0]=HUM
NAME[1]=MYC
NAME[2]=VEC
NAME[3]=PHI
NAME[4]=HEP
NAME[5]=HUH
NAME[6]=JUR
NAME[7]=SEN
NAME[8]=ZIK
NAME[9]=NEW
NAME[10]=EBO
NAME[11]=ECO

rm SDB.fasta
for ii in {0..11}; do
    date
    echo ITERATION ${ii}
    echo INPUT ${IN[${ii}]}
    echo CATEGORY ${NAME[${ii}]}
    cat ${IN[${ii}]} | awk -v X=${NAME[${ii}]} '{if(substr($0,1,1)==">")print ">"X"."substr($0,2);else print $0;}' >> SDB.fasta
done
date
echo DONE
