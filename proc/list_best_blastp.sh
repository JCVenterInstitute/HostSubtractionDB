#!/bin/sh

if [ 0 -eq 1 ]; then
    echo PROTEIN NAMES
    echo AgamP4
    # Convert >AGAP000008-PA DnaJ  to  AGAP000008
    grep '^>' Anopheles-gambiae-PEST_PEPTIDES_AgamP4.4.fa \
	| awk '{print substr($0,2,10) ",\"" substr($0,16) "\"";}' \
	|sort > AgamP4.names.csv
    echo AaegL3
    grep '^>' Aedes-aegypti-Liverpool_PEPTIDES_AaegL3.3.fa \
	| awk '{print substr($0,2,10) ",\"" substr($0,16) "\"";}' \
	|sort > AaegL3.names.csv
    echo AaloF1
    grep '^>' Aedes-albopictus-Foshan_PEPTIDES_AaloF1.1.fa \
	| awk '{print substr($0,2,10) ",\"" substr($0,16) "\"";}' \
	|sort > AaloF1.names.csv
    echo C636
    # convert >gi|1121079054|ref|XP_019524832.1| PREDICTED  to  XP_019524832.1
    grep '^>' C636.RefSeq_Proteins_v2.fa \
	| awk -F '|' '{D=substr($5,2);print $4",\""D"\"";}' \
	|sort > C636.names.csv
fi
    
if [ 1 -eq 1 ]; then
    echo PAIR PER MOSQUITO
    echo best AgamP4
    ./best_blast_hit.pl       < report.blastp.Ref-C636.Qry-AgamP4.tsv |\
         tr '|' '\t' | awk '{print $5,$1;}' > best.C636.perAgamP4.protein.txt 
    echo list AgamP4
    ./list_best_blastp.pl \
	        C636.names.csv AgamP4.names.csv \
	          best.C636.perAgamP4.protein.txt \
	   > pair.best.C636.perAgamP4.csv

    echo best AaloF1
    ./best_blast_hit.pl        < report.blastp.Ref-C636.Qry-AaloF1.tsv |\
          tr '|' '\t' | awk '{print $5,$1;}' > best.C636.perAaloF1.protein.txt 
    echo list AaloF1
    ./list_best_blastp.pl \
	        C636.names.csv AaloF1.names.csv \
	          best.C636.perAaloF1.protein.txt \
	   > pair.best.C636.perAaloF1.csv
    
    echo best AaegL3
    ./best_blast_hit.pl        < report.blastp.Ref-C636.Qry-AaegL3.tsv |\
	  tr '|' '\t' | awk '{print $5,$1;}' > best.C636.perAaegL3.protein.txt  
    echo list AaegL3
    ./list_best_blastp.pl \
	        C636.names.csv AaegL3.names.csv \
	          best.C636.perAaegL3.protein.txt \
	   > pair.best.C636.perAaegL3.csv
fi

if [ 1 -eq 1 ]; then
    echo PAIR PER C636
    echo best AgamP4
    ./best_blast_hit.pl         < report.blastp.Ref-AgamP4.Qry-C636.tsv |\
	  tr '|' '\t' | awk '{print $4,$5;}' > best.AgamP4.perC636.protein.txt  
    echo list AgamP4
    ./list_best_blastp.pl \
	       C636.names.csv AgamP4.names.csv \
	                 best.AgamP4.perC636.protein.txt \
                  > pair.best.AgamP4.perC636.csv
    
    echo best AaegL3
    ./best_blast_hit.pl         < report.blastp.Ref-AaloF1.Qry-C636.tsv |\
	  tr '|' '\t' | awk '{print $4,$5;}' > best.AaloF1.perC636.protein.txt  
    echo list AaegL3
    ./list_best_blastp.pl \
	       C636.names.csv AaloF1.names.csv \
	                 best.AaloF1.perC636.protein.txt \
                  > pair.best.AaloF1.perC636.csv
    
    echo best AaloF1
    ./best_blast_hit.pl         < report.blastp.Ref-AaegL3.Qry-C636.tsv |\
	  tr '|' '\t' | awk '{print $4,$5;}' > best.AaegL3.perC636.protein.txt  
    echo list AaloF1
    ./list_best_blastp.pl \
	       C636.names.csv AaegL3.names.csv \
	                 best.AaegL3.perC636.protein.txt \
                  > pair.best.AaegL3.perC636.csv
fi

echo DONE
date
