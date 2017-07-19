#!/bin/sh

if [ $# != 1 ]; then
    echo "Usage: $0 <blast_report>"
    exit 1
fi
BLAST_REPORT=$1
echo BLAST_REPORT $BLAST_REPORT
OUTFMT="7 qseqid sseqid pident length mismatch gapopen evalue bitscore staxids"
echo EXPECTED BLAST OUTPUT FORMAT $OUTFMT

HERE=`pwd`
echo "CURRENT DIRECTORY ${HERE}"
BASE=`basename "$0"`
echo "THIS PROGRAM IS ${BASE}"
DIR=`dirname "$0"`
if [ "${DIR}" == "." ]; then
    echo "ERROR SCRIPT MUST BE INVOKED WITH ABSOLUTE PATH"
    exit 2
fi
echo "SCRIPT DIRECTORY IS ${DIR}"
# This is a derivative of the NCBI file, containing ID TAXON, sorted, possibly with duplicate IDs
# cat names.dmp | awk -F '|' '{if ($1 != p) print $1, $2 ; p=$1;}' | sort -k1,1 > ncbi_taxonomy.taxdb.id_name.txt
TAXDB=ncbi_taxonomy.taxdb.id_name.txt

TEMPFILE="${BLAST_REPORT}.tmp"
OUTFILE="${BLAST_REPORT}.taxa"
echo MAKE TEMP FILE $TEMPFILE
cat $BLAST_REPORT |\
    grep -v "^# 0 hits found" |\
    grep -A 1 "^# .* hits found" |\
    grep -v "hits found" | grep -v "^--" |\
    cut -f 9 | sort | uniq -c | sort -k2,2 > $TEMPFILE
echo MAKE OUT FILE $OUTFILE
join -1 2 ${TEMPFILE} ${DIR}/${TAXDB} | sort -k2,2n > $OUTFILE

echo "OUTPUT FIELDS: TAXON_ID READ_COUNT TAXON_NAME"
echo done

