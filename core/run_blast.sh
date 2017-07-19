#!/bin/sh
# Run blast.
# This requires ~16GB RAM

if [ $# != 1 ]; then
    echo "Usage: $0 <reads>"
    exit 1
fi
READS=$1
echo READS $READS

export BLASTN="/usr/local/packages/ncbi-blast+/bin/blastn"
echo BLASTN $BLASTN
echo BLAST VERSION
${BLASTN} -version
export BLASTDB=/usr/local/projdata/99999/IFX/CommonDB/taxdb/taxdb_current/:/usr/local/projdata/99999/IFX/CommonDB/nt/nt_current
echo BLASTDB $BLASTDB
THREADS=4
TARGETS_PER_QUERY=1

# ask for the taxon ID
# at most 1 hit per read
# Note: With >1 hit per read, we see that hits are non-specific. One read hit all sorts of parasites.

date
echo COMMAND
OUTFMT="7 qseqid sseqid pident length mismatch gapopen evalue bitscore staxids"
echo BLAST OUTPUT FORMAT $OUTFMT
echo ${BLASTN} -db nt -query $READS -out ${READS}.blast -outfmt "${OUTFMT}" -num_threads ${THREADS} -max_target_seqs ${TARGETS_PER_QUERY} 
${BLASTN} -db nt -query $READS -out ${READS}.blast -outfmt "${OUTFMT}" -num_threads ${THREADS} -max_target_seqs ${TARGETS_PER_QUERY} 

echo -n $?; echo " exit status"
date
echo DONE

