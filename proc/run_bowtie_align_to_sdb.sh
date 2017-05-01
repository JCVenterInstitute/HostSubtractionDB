#!/bin/sh
# Run bowtie2.
# Align raw reads.
# Align to all references.

# Special case handling required for timmomatic outputs:
# Process the R1 and R2 that are still paired.
# Process the R1 and R2 that are no longer paired.

HERE=`pwd`
THIS=`basename "$0"`
SDB_HOME=/local/ifs2_projdata/8370/projects/DHSSDB/GitHubRepo/HostSubtractionDB
SDB_UTIL=${SDB_HOME}/proc/subtract_mapped_reads.sh

# This is the command to align reads using the index.
# The command has 'small' and 'large' binaries.
# Use small unless the target genome is > 4GB.
# Our target genome is about 3 GB.
export BOWTIE_ALIGN=/usr/local/bin/bowtie2
export BOWTIE_BUILD=/usr/local/bin/bowtie2-build
#export BOWTIE_ALIGN=/usr/local/bin/bowtie2-align-s
#export BOWTIE_ALIGN=/usr/local/bin/bowtie2-align-l

# This is the command for grid submissions
QSUB=/usr/local/sge_current/bin/lx-amd64/qsub

# Enable python on the grid
source /usr/local/common/env/python.sh

# This is the command for sorting bam files
SAMTOOLS=/usr/local/bin/samtools

function runit () {
    echo "RUN COMMAND"
    echo ${CMD}
    date
    nice ${CMD}
    echo -n $?;    echo " exit status"
    date
}

if [ -z ${SGE_TASK_ID} ]; then
    # This means we are not on the SGE grid.
    # When run for the first time, create the index and submit alignments to the grid.
    JOB=0
else
    # This means we are running on the grid.
    # When run on the grid, start an alignment job.
    JOB=${SGE_TASK_ID}
fi

INDEX=SDB
MYINDEX=${INDEX}

if [ ${JOB} -eq 0 ]; then
    # Create the index, submit alignments to the grid, and stop.
    echo "Bowtie index"
    CMD="${BOWTIE_BUILD} phiX.fasta,UniVec.fasta,clean.hepg2_unmapped_contigs.fasta,clean.huh7_unmapped_contigs.fasta,clean.jurkat_unmapped_contigs.fasta,Mycoplasma.remove_blanks.fasta,GCA_000001405.22_GRCh38.p7_genomic.fasta ${INDEX}"
    runit
    echo "Bowtie align on grid"
    CMD="${QSUB} -cwd -b n -A DHSSDB -P 8370 -N human1 -pe threaded 4 -l medium -l memory=1g -t 1-12 -j y -o $HERE ${HERE}/${THIS}"
    echo $CMD
    $CMD
    echo "Jobs sent to grid"
    exit
fi

# Do this when using the -t option on qsub.
echo SGE_TASK_ID $SGE_TASK_ID
JOB=$SGE_TASK_ID
echo JOB $JOB

# Note bowtie2 can take fastq or fastq.gz. Assume all gz.
# If given a pair, bowtie will insist that both reads map and the pair is concordant.
# By default, bowtie uses "mixed mode", aligning separately if pair align fails.
# For maximum sensitivity, map reads separately, not as pairs.

R1[1]=5JURKATRIBO_S11_R1_001.fastq
R2[1]=5JURKATRIBO_S11_R2_001.fastq
R1[2]=5JURKAT_S5_R1_001.fastq
R2[2]=5JURKAT_S5_R2_001.fastq
R1[3]=6JURKATCELLSSEVRIBOR_S12_R1_001.fastq
R2[3]=6JURKATCELLSSEVRIBOR_S12_R2_001.fastq
R1[4]=6JURKATSEV_S6_R1_001.fastq
R2[4]=6JURKATSEV_S6_R2_001.fastq

R1[5]=1HEPG2RIBO_S7_R1_001.fastq
R2[5]=1HEPG2RIBO_S7_R2_001.fastq
R1[6]=1HEPG2_S1_R1_001.fastq
R2[6]=1HEPG2_S1_R2_001.fastq
R1[7]=2HEPG2SEVRIBO_S8_R1_001.fastq
R2[7]=2HEPG2SEVRIBO_S8_R2_001.fastq
R1[8]=2HEPG2SEV_S2_R1_001.fastq
R2[8]=2HEPG2SEV_S2_R2_001.fastq

R1[9]=3HUH7RIBO_S9_R1_001.fastq
R2[9]=3HUH7RIBO_S9_R2_001.fastq
R1[10]=3HUH7_S3_R1_001.fastq
R2[10]=3HUH7_S3_R2_001.fastq
R1[11]=4HUH7SEVRIBO_S10_R1_001.fastq
R2[11]=4HUH7SEVRIBO_S10_R2_001.fastq
R1[12]=4HUH7SEV_S4_R1_001.fastq
R2[12]=4HUH7SEV_S4_R2_001.fastq

BASE[1]=JUR.RIBO.MOCK
BASE[2]=JUR.NONE.MOCK
BASE[3]=JUR.RIBO.SEND
BASE[4]=JUR.NONE.SEND
BASE[5]=HEP.RIBO.MOCK
BASE[6]=HEP.NONE.MOCK
BASE[7]=HEP.RIBO.SEND
BASE[8]=HEP.NONE.SEND
BASE[9]=HUH.RIBO.MOCK 
BASE[10]=HUH.NONE.MOCK
BASE[11]=HUH.RIBO.SEND
BASE[12]=HUH.NONE.SEND

MYBASE=${BASE[${JOB}]}
SAM=${MYBASE}.sam
BAM=${MYBASE}.bam
MYR1=${R1[${JOB}]}.gz
MYR2=${R2[${JOB}]}.gz
MYNAMES=${MYBASE}.read_names
echo MYINDEX $MYINDEX
echo MYR1 $MYR1
echo MYR2 $MYR2
echo MYBASE $MYBASE
echo SAM $SAM
echo BAM $BAM
echo MYNAMES $MYNAMES

THREADS="-p 4"
#THREADS="-p 4 -mm"
#ALIGNMENT="--end-to-end "
ALIGNMENT="--sensitive-local"
FASTQ="-q  --phred33"

# We have dectected a bug in bowtie2.
# If we use the "--un" option with cutadapt files (which can have zero-length reads),
# then bowtie writes fastq files that are corrupt in subtle ways.

UNALIGNED="--un R1.${MYBASE}.u.fq"   # warning: bug
UNALIGNED="--un R2.${MYBASE}.u.fq"
UNALIGNED="--no-unal"                # keep unaligned out of the sam file
CMD="${BOWTIE_ALIGN} ${UNALIGNED} ${THREADS} ${ALIGNMENT} ${FASTQ} -x ${MYINDEX} -U ${MYR1} -S R1.${SAM}"
runit
CMD="${BOWTIE_ALIGN} ${UNALIGNED} ${THREADS} ${ALIGNMENT} ${FASTQ} -x ${MYINDEX} -U ${MYR2} -S R2.${SAM}"
runit

echo "CONVERT SAM TO BAM"
THREADS="-@ 4"
TEMP=TEMP.${JOB}
CMD="${SAMTOOLS} sort -T ${TEMP} -l 9 -m 1G ${THREADS} -o R1.${BAM} R1.${SAM} "
CMD="${SAMTOOLS} view -h -b -o R1.${BAM} R1.${SAM} "
runit
CMD="${SAMTOOLS} sort -T ${TEMP} -l 9 -m 1G ${THREADS} -o R2.${BAM} R2.${SAM} "
CMD="${SAMTOOLS} view -h -b -o R2.${BAM} R2.${SAM} "
runit

echo "GET MAPPED READ NAMES"
# If necessary, use a command like this
# sed 's/\/1$//'
# to strip off any read-specific suffix like /1 or /2
# ${SAMTOOLS} view R1.${BAM} | cut -f 1  > ${MYNAMES}

echo "SUBTRACT MAPPED READS"
${SAMTOOLS} view R1.${BAM} | cut -f 1  > ${MYNAMES}
CMD="${SDB_UTIL} ${MYR1} ${MYNAMES} nonhost.${MYBASE}.R1.fastq"
runit
${SAMTOOLS} view R2.${BAM} | cut -f 1  > ${MYNAMES}
CMD="${SDB_UTIL} ${MYR2} ${MYNAMES} nonhost.${MYBASE}.R2.fastq"
runit

echo "COMPRESS FASTQ"
CMD="gzip -v nonhost.${MYBASE}.R1.fastq"
runit
CMD="gzip -v nonhost.${MYBASE}.R2.fastq"
runit

echo "DONE"
