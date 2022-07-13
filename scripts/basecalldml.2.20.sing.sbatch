#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -J bcl2fastq.2.20.sing
#SBATCH --qos vip
#SBATCH -n 20

set -uE
trap ' echo Error $? occured on $LINENO && exit 1' ERR

#load the nextflow module
module load picard/2.26.1
#NXF_OPTS='-Xms1g -Xmx4g'
#WKDIR=/var/scratch/$USER/work_nano
#NXF_SINGULARITY_CACHEDIR=${WKDIR}/singularity
#export NXF_SINGULARITY_CACHEDIR=${WKDIR}/singularity

#Setting work Directory in HPC
INDIR=$PWD
src_dir=${INDIR##*/}
DATADIR=${INDIR}/../core_data
WORKDIR=/var/scratch/$USER/"`date +"%Y.%m.%d"`_${src_dir}"
src_dir1=${WORKDIR##*/}
until [[ -d ${WORKDIR} ]]
do
	echo "Creating output directory '${WORKDIR}'"
	mkdir -p ${WORKDIR}/demultiplexed/Stats
done
cd $WORKDIR

echo -e "\tUsing $WORKDIR/demultiplexed on $SLURMD_NODENAME"
#echo -e "\tSetting up the working environment installations"
#rsync -u --size-only "${ENVDIR}" "${ENVDIR_WKDIR}"
#source ${INDIR}/../../scripts/process_files.sh

#generating the samplesheet.csv file: Please provide the forward and reverse read suffixes below:

#run basecalling job on N (SBATCH -ni) threads/cores

#java -jar /export/apps/picard/2.26.1/picard.jar 	


singularity run ${INDIR}/../../scripts/bcl2fastq/singularity_bcl2fastq_latest.sif \
	-i ${INDIR}/220322_NB552490_0015_AHWY2JBGXK/Data/Intensities/BaseCalls/ \
	--runfolder-dir ${INDIR}/220322_NB552490_0015_AHWY2JBGXK/ \
	--loading-threads 8 \
	--processing-threads 4 \
	--writing-threads 4 \
	--minimum-trimmed-read-length 35 \
	--mask-short-adapter-reads 22 \
	--adapter-stringency 0.9 \
	--create-fastq-for-index-reads


if [ $? -eq 0 ]
then
        echo -e "\tThe results have been stored in '${WORKDIR}' in ${SLURMD_NODENAME}.\n\tCopying the results to ${INDIR} ..."
        cp -rf ${WORKDIR}/ ${INDIR}/
        if [ $? -eq 0 ]
        then
                echo -e "\tCopying successful..."
        elif [ $? -ne 0 ]
        then
                echo -e "\tCopying was NOT successful..."
        fi
fi
