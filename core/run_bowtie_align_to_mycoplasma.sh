#!/bin/sh
# Run bowtie2 to align reads to MYCOPLASMA.
# Assume the R1 and R2 reads of a pair have the same read name
# Use this assumption to remove a read if its mate mapped to the MYCOPLASMA.

HERE=`pwd`
THIS=`basename "$0"`
OUTPUT="nonmyco"
rm -v ${OUTPUT}.*.fastq.gz

export BOWTIE_ALIGN=/usr/local/bin/bowtie2
export BOWTIE_BUILD=/usr/local/bin/bowtie2-build

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

if [ ${JOB} -eq 0 ]; then
    # Create the index, submit alignments to the grid, and stop.
    date
    echo "Bowtie index"
    CMD="${BOWTIE_BUILD} Mycoplasma.genomes.fasta Mycoplasma"
 ##   runit
    
    echo "Bowtie align on grid"
    CMD="${QSUB} -cwd -b n -A DHSSDB -P 8370 -N humanCL -pe threaded 4 -l medium -l memory=1g -t 1-12 -j y -o $HERE ${HERE}/${THIS}"
    echo $CMD
    $CMD
    echo "Jobs sent to grid"
    exit
fi

# Do this when using the -t option on qsub.
echo SGE_TASK_ID $SGE_TASK_ID
JOB=$SGE_TASK_ID
echo JOB $JOB

# Note bowtie2 can take fastq or fastq.gz. Assume gz.
# If given a pair, bowtie will insist that both reads map and the pair is concordant.
# By default, bowtie uses "mixed mode", aligning separately if pair align fails.
# For maximum sensitivity, map reads separately, not as pairs.

R1[1]=nonvec.HEPG2NONE.R1.fastq
R1[2]=nonvec.HEPG2RIBO.R1.fastq
R1[3]=nonvec.HEPG2SEVNONE.R1.fastq
R1[4]=nonvec.HEPG2SEVRIBO.R1.fastq
R1[5]=nonvec.HUH7NONE.R1.fastq
R1[6]=nonvec.HUH7RIBO.R1.fastq
R1[7]=nonvec.HUH7SEVNONE.R1.fastq
R1[8]=nonvec.HUH7SEVRIBO.R1.fastq
R1[9]=nonvec.JURKATCELLSSEVRIBO.R1.fastq
R1[10]=nonvec.JURKATNONE.R1.fastq
R1[11]=nonvec.JURKATRIBO.R1.fastq
R1[12]=nonvec.JURKATSEVNONE.R1.fastq

R2[1]=nonvec.HEPG2NONE.R2.fastq
R2[2]=nonvec.HEPG2RIBO.R2.fastq
R2[3]=nonvec.HEPG2SEVNONE.R2.fastq
R2[4]=nonvec.HEPG2SEVRIBO.R2.fastq
R2[5]=nonvec.HUH7NONE.R2.fastq
R2[6]=nonvec.HUH7RIBO.R2.fastq
R2[7]=nonvec.HUH7SEVNONE.R2.fastq
R2[8]=nonvec.HUH7SEVRIBO.R2.fastq
R2[9]=nonvec.JURKATCELLSSEVRIBO.R2.fastq
R2[10]=nonvec.JURKATNONE.R2.fastq
R2[11]=nonvec.JURKATRIBO.R2.fastq
R2[12]=nonvec.JURKATSEVNONE.R2.fastq

BASE[1]=HEPG2NONE
BASE[2]=HEPG2RIBO
BASE[3]=HEPG2SEVNONE
BASE[4]=HEPG2SEVRIBO
BASE[5]=HUH7NONE
BASE[6]=HUH7RIBO
BASE[7]=HUH7SEVNONE
BASE[8]=HUH7SEVRIBO
BASE[9]=JURKATCELLSSEVRIBO
BASE[10]=JURKATNONE
BASE[11]=JURKATRIBO
BASE[12]=JURKATSEVNONE

INDEX[1]=Mycoplasma
INDEX[2]=Mycoplasma
INDEX[3]=Mycoplasma
INDEX[4]=Mycoplasma
INDEX[5]=Mycoplasma
INDEX[6]=Mycoplasma
INDEX[7]=Mycoplasma
INDEX[8]=Mycoplasma
INDEX[9]=Mycoplasma
INDEX[10]=Mycoplasma
INDEX[11]=Mycoplasma
INDEX[12]=Mycoplasma

MYBASE=${BASE[${JOB}]}
SAM=${MYBASE}.sam
BAM=${MYBASE}.bam
MYINDEX=${INDEX[${JOB}]}
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
# Later, use a command like this
# sed 's/\/1$//'
# to strip off any read-specific suffix like /1 or /2
${SAMTOOLS} view R1.${BAM} | cut -f 1  > ${MYNAMES}
${SAMTOOLS} view R2.${BAM} | cut -f 1 >> ${MYNAMES}

echo "SUBTRACT MAPPED READS"
SDB_HOME=/local/ifs2_projdata/8370/projects/DHSSDB/GitHubRepo/HostSubtractionDB
SDB_UTIL=${SDB_HOME}/proc/subtract_mapped_reads.sh
CMD="${SDB_UTIL} ${MYR1} ${MYNAMES} ${OUTPUT}.${MYBASE}.R1.fastq"
runit
CMD="${SDB_UTIL} ${MYR2} ${MYNAMES} ${OUTPUT}.${MYBASE}.R2.fastq"
runit

echo "COMPRESS FASTQ"
CMD="gzip -v ${OUTPUT}.${MYBASE}.*.fastq"
runit

echo "DONE"
