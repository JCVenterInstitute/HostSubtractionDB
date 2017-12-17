#!/bin/sh

# Here is a project-specific pipeline in script form.
# This script provides an example but it is specific to one dataset and one compute environment.
# The dataset is a subset of the JCVI data using Sendai virus.
# The compute environment is the SGE grid at JCVI.

SCRIPTDIR=/local/ifs2_projdata/8370/projects/DHSSDB/GitHubRepo/HostSubtractionDB/core
SCRIPTDIR=$1
echo "ASSUME CORE SCRIPT DIRECTORY IS ${SCRIPTDIR}"
if [ -z "$1" ]
then
    echo "Usage: $0 <core>"
    exit 1
fi

SDB=subtraction.HUH7.contigs.fasta  # for the demo
INDEX=SDB-HUH7                      # for the demo

TRIMMER="run_trimmomatic.sh"
BOWTIE_BUILD="/usr/local/bin/bowtie2-build"
BOWTIE_ALIGN=/usr/local/bin/bowtie2
SAMTOOLS=/usr/local/bin/samtools
FILTER=fastq-filter-by-name.pl

date
echo "Clean up any output from previous runs."
rm -v trim.pair.*.fastq
rm -v trim.sing.*.fastq
rm -v nonSDB.*.fastq
rm -v sdb.*.sam sdb.*.bam

echo
echo "Index the subtraction database."
echo "This uses bowtie2."
echo "This requires RAM 2 x database size."
echo "This reads FASTA and generates *.bt2 files."
ls ${BOWTIE_BUILD} 
${BOWTIE_BUILD} ${SDB} ${INDEX} > build.out 2>&1
echo -n $?; echo " exit status"
echo RESULTS
ls ${INDEX}
date

echo
echo "Trim the RNAseq reads."
echo "This uses trimmomatic which requires Java."
echo "This requires an adapter file corresponding to the Illumina library prep."
echo "Several adapter files are provided with trimmomatic."
echo "To avoid re-trimming, this will clean up trim.*.fastq files from previous runs."
echo "This runs on all files named *_R1_*.fastq and looks for the corresponding R2."
echo "This generates files named trim.pair.*.fastq which get mapped next."
echo "This also generates files named trim.sing.*.fastq which get ignored."
ls ${SCRIPTDIR}/${TRIMMER}
${SCRIPTDIR}/${TRIMMER} > trim.out 2>&1
echo RESULTS
ls trim.pair.*.fastq
date

echo
echo "Looop over pairs of trimmed read files."
ls ${BOWTIE_ALIGN} 
THREADS="-p 4"
ALIGNMENT="--sensitive-local"
FASTQ="-q  --phred33"
UNALIGNED="--no-unal"                # keep unaligned out of the sam file
EXCLUDE="-v"  # option to exclude the named reads
for FF in trim.pair.*_R1_*.fastq; do    
    MYR1=${FF}
    MYR2=` echo $MYR1 | sed 's/_R1_/_R2_/'  `
    SAM="${MYR1}.sam"
    BAM="${MYR1}.bam"
    TEMP_IDS="${MYR1}.tmp"
    echo
    echo "Map trimmed reads to SDB."
    echo "Writinng $BAM"
    CMD="${BOWTIE_ALIGN} ${UNALIGNED} ${THREADS} ${ALIGNMENT} ${FASTQ} -x ${INDEX} -1 ${MYR1} -2 ${MYR2} -S sdb.${SAM}"
    echo ${CMD}
    nice ${CMD} > map.out 2>&1
    echo -n $?;    echo " exit status"
    CMD="${SAMTOOLS} view -h -b -o sdb.${BAM} sdb.${SAM} "
    echo ${CMD}
    nice ${CMD} 
    echo -n $?;    echo " exit status"
    echo "Write read ID for every read mapped."
    echo "These are the IDs to subtract."
    echo "Assume reads 1 and 2 have the same read ID."
    echo "If your reads have the /1 and /2 suffix, please change that and start again."
    ${SAMTOOLS} view sdb.${BAM} | cut -f 1 >  ${TEMP_IDS}
    perl ${SCRIPTDIR}/${FILTER} ${EXCLUDE} ${TEMP_IDS} < ${MYR1} > nonSDB.${MYR1} 2> subtract.R1.out
    perl ${SCRIPTDIR}/${FILTER} ${EXCLUDE} ${TEMP_IDS} < ${MYR2} > nonSDB.${MYR2} 2> subtract.R2.out
done
echo
echo RESULTS
ls nonSDB.*
date

exit



