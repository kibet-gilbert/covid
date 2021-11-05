#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH --qos vip
#SBATCH -J basecallalbacore_sing
#SBATCH -n 30

set -uE
trap ' echo Error $? occured on $LINENO && exit 1' ERR

#load the nextflow module
#NXF_OPTS='-Xms1g -Xmx4g'
#WKDIR=/var/scratch/$USER/work_nano
#NXF_SINGULARITY_CACHEDIR=${WKDIR}/singularity
#export NXF_SINGULARITY_CACHEDIR=${WKDIR}/singularity

#Setting work Directory in HPC
INDIR=$PWD
src_dir=${INDIR##*/}
DATADIR=${INDIR}/../core_data
WORKDIR=${INDIR}/data/EXP-NBD196_basecall
src_dir1=${WORKDIR##*/}
until [[ -d ${WORKDIR} ]]
do
	echo "Creating output directory '${WORKDIR}'"
	mkdir -p ${WORKDIR}
done

#echo -e "\tUsing $WORKDIR/fast5_pass on $SLURMD_NODENAME"
#echo -e "\tSetting up the working environment installations"
#rsync -u --size-only "${ENVDIR}" "${ENVDIR_WKDIR}"
#source ${INDIR}/../../scripts/process_files.sh

#generating the samplesheet.csv file: Please provide the forward and reverse read suffixes below:

#run basecalling job on N (SBATCH -ni) threads/cores
singularity run --app read_fast5_basecaller.py ${INDIR}/albacore.sif \
	-f FLO-MIN106 \
	-k SQK-LSK109 \
	-i ${INDIR}/data/EXP-NBD196 \
	-t 30 \
	-s ${INDIR}/data/EXP-NBD196_basecall \
	-o fastq \
	-b \
	--disable_filtering \
	--barcoding \
	--recursive
