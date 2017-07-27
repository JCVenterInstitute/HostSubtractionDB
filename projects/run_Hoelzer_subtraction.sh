#!/bin/sh

# Here is a project-specific pipeline in script form.
# This script provides an example but it is specific to one dataset and one compute environment.
# The dataset is the Hoelzer data from HUH7 cells infected with Ebola and Marv.
#    Hoelzer et al (2016) Scientific Reports.  DOI: 10.1038/srep34589 
#    "Differential transcriptional responses to Ebola and Marburg virus infection in bat and human cells"
# The compute environment is the SGE grid at JCVI.

HUMAN_INDEX=HumanRef
CELL_INDEX=subtraction.huh7.hepg2.jurkat
CONTAM_INDEX=ecoli.myco.univec.phix
PWD=`pwd`

if [ $# -ne 1 ]
then
    echo "ERROR. MISSING PARAMETER: 1=human, 2=cell-line, 3=contam"
    exit 1
fi
PARAMETER=$1
echo PARAMETER $PARAMETER

if [ "${PARAMETER}" == 1 ]
   then
       if [ ! -f "${HUMAN_INDEX}.1.bt2" ]
       then
	   echo "ERROR. Missing human index"
	   exit 1
       fi
fi

if [ "${PARAMETER}" == 2 ]
then
    if [ ! -f "${CELL_INDEX}.1.bt2" ]
    then
	echo "ERROR. Missing cell-line index"
	exit 1
    fi
fi

if [ "${PARAMETER}" == 3 ]
then
    if [ ! -f "${CONTAM_INDEX}.1.bt2" ]
    then
	echo "ERROR. Missing contam index"
	exit 1
    fi
fi

SCRIPTDIR=/local/ifs2_projdata/8370/projects/DHSSDB/GitHubRepo/HostSubtractionDB/core
SCRIPT=${SCRIPTDIR}/run_bowtie_align_R1thenR2.sh
THREADS=4
echo SCRIPT $SCRIPT
echo THREADS $THREADS
PATTERN="*_R1_*.fastq"
QSUB="qsub -cwd -b n -A DHSSDB -P 8370 -N subtr -pe threaded ${THREADS} -l medium -l memory=4g -j y -o ${PWD}"

for FF in ${PATTERN}; do
    # get basename without fastq suffix
    R1=` echo ${FF} | sed 's/.fastq//' `
    # get basename for read#2
    R2=` echo ${R1} | sed 's/_R1_/_R2_/' `
    echo R1 $R1 R2 $R2
    if [ "${PARAMETER}" == 1 ]
    then	
	${QSUB} ${SCRIPT} ${SCRIPTDIR} ${R1} ${R2} ${HUMAN_INDEX}  "human"     ${THREADS} 
    elif [ "${PARAMETER}" == 2 ]
    then
	${QSUB} ${SCRIPT} ${SCRIPTDIR} ${R1} ${R2} ${CELL_INDEX}   "cell-line" ${THREADS} 
    elif [ "${PARAMETER}" == 3 ]
    then
	${QSUB} ${SCRIPT} ${SCRIPTDIR} ${R1} ${R2} ${CONTAM_INDEX} "contam"    ${THREADS}
    fi    
done

echo "NEW FASTQ WILL BE CREATED."
echo "MOVE OLD FASTQ BEFORE STARTING NEXT ROUND"
