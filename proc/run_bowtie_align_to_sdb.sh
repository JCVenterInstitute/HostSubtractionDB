#!/bin/sh
# Run bowtie2 to align reads to SDB.
# SDB = Subtraction DataBase = reference human + cell-line-specific-contigs + phiX.

HERE=`pwd`

# This is the command to align reads using the index.
# The command has 'small' and 'large' binaries.
# Use small unless the target genome is > 4GB.
# Our target genome is about 3 GB.
export BOWTIE_ALIGN=/usr/local/bin/bowtie2
#export BOWTIE_ALIGN=/usr/local/bin/bowtie2-align-s
#export BOWTIE_ALIGN=/usr/local/bin/bowtie2-align-l

# This is the command for grid submissions
QSUB=/usr/local/sge_current/bin/lx-amd64/qsub

# Enable python on the grid
source /usr/local/common/env/python.sh

# This is the command for sorting bam files
SAMTOOLS=/usr/local/bin/samtools

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
    CMD="${BOWTIE_BUILD} ${MYTARGET} ${MYINDEX}"
    #    echo CMD $CMD
    #    COMMENDTED OUT BECAUSE WE ALREADY HAVE THE INDEXES
    #    ${CMD}
    #    echo -n $?; echo " exit status"
    date
    echo "Bowtie align on grid"
    CMD="${QSUB} -cwd -b n -A DHSSDB -P 8370 -N humanCL -pe threaded 4 -l medium -l memory=1g -t 1-12 -j y -o $HERE ${HERE}/run_bowtie_align_to_sdb.sh"
    echo $CMD
    $CMD
    echo "Jobs sent to grid"
    exit
fi

# Do this when using the -t option on qsub.
echo SGE_TASK_ID $SGE_TASK_ID
JOB=$SGE_TASK_ID
echo JOB $JOB

# Note bowtie2 can take fastq or fastq.gz
# If given a pair, bowtie will insist that both reads map and the pair is concordant.
# By default, bowtie uses "mixed mode", aligning separately if pair align fails.
# For maximum sensitivity, map reads separately, not as pairs.

R1[1]=cutadapt.5JURKATRIBO_S11_R1_001.fastq.gz
R2[1]=cutadapt.5JURKATRIBO_S11_R2_001.fastq.gz
R1[2]=cutadapt.5JURKAT_S5_R1_001.fastq.gz
R2[2]=cutadapt.5JURKAT_S5_R2_001.fastq.gz
R1[3]=cutadapt.6JURKATCELLSSEVRIBOR_S12_R1_001.fastq.gz
R2[3]=cutadapt.6JURKATCELLSSEVRIBOR_S12_R2_001.fastq.gz
R1[4]=cutadapt.6JURKATSEV_S6_R1_001.fastq.gz
R2[4]=cutadapt.6JURKATSEV_S6_R2_001.fastq.gz

R1[5]=cutadapt.1HEPG2RIBO_S7_R1_001.fastq.gz
R2[5]=cutadapt.1HEPG2RIBO_S7_R2_001.fastq.gz
R1[6]=cutadapt.1HEPG2_S1_R1_001.fastq.gz
R2[6]=cutadapt.1HEPG2_S1_R2_001.fastq.gz
R1[7]=cutadapt.2HEPG2SEVRIBO_S8_R1_001.fastq.gz
R2[7]=cutadapt.2HEPG2SEVRIBO_S8_R2_001.fastq.gz
R1[8]=cutadapt.2HEPG2SEV_S2_R1_001.fastq.gz
R2[8]=cutadapt.2HEPG2SEV_S2_R2_001.fastq.gz

R1[9]=cutadapt.3HUH7RIBO_S9_R1_001.fastq.gz
R2[9]=cutadapt.3HUH7RIBO_S9_R2_001.fastq.gz
R1[10]=cutadapt.3HUH7_S3_R1_001.fastq.gz
R2[10]=cutadapt.3HUH7_S3_R2_001.fastq.gz
R1[11]=cutadapt.4HUH7SEVRIBO_S10_R1_001.fastq.gz
R2[11]=cutadapt.4HUH7SEVRIBO_S10_R2_001.fastq.gz
R1[12]=cutadapt.4HUH7SEV_S4_R1_001.fastq.gz
R2[12]=cutadapt.4HUH7SEV_S4_R2_001.fastq.gz

INDEX[1]=JURKAT
INDEX[2]=JURKAT
INDEX[3]=JURKAT
INDEX[4]=JURKAT
INDEX[5]=HEPG2
INDEX[6]=HEPG2
INDEX[7]=HEPG2
INDEX[8]=HEPG2
INDEX[9]=HUH7
INDEX[10]=HUH7
INDEX[11]=HUH7
INDEX[12]=HUH7

BASE[1]=JURKATRIBO
BASE[2]=JURKATNONE
BASE[3]=JURKATCELLSSEVRIBO
BASE[4]=JURKATSEVNONE
BASE[5]=HEPG2RIBO
BASE[6]=HEPG2NONE
BASE[7]=HEPG2SEVRIBO
BASE[8]=HEPG2SEVNONE
BASE[9]=HUH7RIBO 
BASE[10]=HUH7NONE
BASE[11]=HUH7SEVRIBO
BASE[12]=HUH7SEVNONE

MYBASE=${BASE[${JOB}]}
SAM=${MYBASE}.sam
BAM=${MYBASE}.bam

MYINDEX=${INDEX[${JOB}]}
MYR1=${R1[${JOB}]}
MYR2=${R2[${JOB}]}
echo MYINDEX $MYINDEX
echo MYR1 $MYR1
echo MYR2 $MYR2
echo MYBASE $MYBASE
echo SAM $SAM
echo BAM $BAM

function runit () {
    echo "RUN COMMAND"
    echo ${CMD}
    date
    nice ${CMD}
    echo -n $?;    echo " exit status"
    date
}

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
##runit
CMD="${BOWTIE_ALIGN} ${UNALIGNED} ${THREADS} ${ALIGNMENT} ${FASTQ} -x ${MYINDEX} -U ${MYR2} -S R2.${SAM}"
##runit

echo "CONVERT SAM TO BAM"
THREADS="-@ 4"
TEMP=TEMP.${JOB}
CMD="${SAMTOOLS} sort -T ${TEMP} -l 9 -m 1G ${THREADS} -o R1.${BAM} R1.${SAM} "
CMD="${SAMTOOLS} view -h -b -o R1.${BAM} R1.${SAM} "
runit
CMD="${SAMTOOLS} sort -T ${TEMP} -l 9 -m 1G ${THREADS} -o R2.${BAM} R2.${SAM} "
CMD="${SAMTOOLS} view -h -b -o R2.${BAM} R2.${SAM} "
runit

SDB_HOME=/local/ifs2_projdata/8370/projects/DHSSDB/GitHubRepo/HostSubtractionDB
SDB_UTIL=${SDB_HOME}/proc/subtract_mapped_reads.sh
CMD=${SDB_UTIL} ${MYR1} R1.${BAM}
runit
CMD=${SDB_UTIL} ${MYR2} R2.${BAM}
runit

echo "OK TO DELETE SAM ONCE BAM HAS TESTED OK"
echo "DONE"
