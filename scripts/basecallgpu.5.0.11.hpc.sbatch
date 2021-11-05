#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -J guppygpu.5.0.11
#SBATCH --qos vip
#SBATCH --gres=gpu:v100:1
#SBATCH --nodelist=compute06
#SBATCH -n 10

set -uE
trap ' echo Error $? occured on $LINENO && exit 1' ERR

module load guppy/5.0.11
#load the nextflow module
#NXF_OPTS='-Xms1g -Xmx4g'
#WKDIR=/var/scratch/$USER/work_nano
#NXF_SINGULARITY_CACHEDIR=${WKDIR}/singularity
#export NXF_SINGULARITY_CACHEDIR=${WKDIR}/singularity

#Setting work Directory in HPC
INDIR=${PWD}
src_dir=${INDIR##*/}
DATADIR=${INDIR}/../core_data
WORKDIR=/var/scratch/$USER/"`date +"%Y.%m.%d"`_${src_dir}"
src_dir1=${WORKDIR##*/}
#FLOWCELL="FLO-FLG001"
FLOWCELL="FLO-MIN106"
until [[ -d ${WORKDIR} ]]
do
	echo "Creating output directory '${WORKDIR}'"
	mkdir -p ${WORKDIR}/fast5_pass
done
cd ${WORKDIR}

echo -e "\tUsing ${WORKDIR}/fast5_pass on ${SLURMD_NODENAME}:${SLURM_JOBID}"
#echo -e "\tSetting up the working environment installations"
#rsync -u --size-only "${ENVDIR}" "${ENVDIR_WKDIR}"
#source ${INDIR}/../../scripts/process_files.sh

#generating the samplesheet.csv file: Please provide the forward and reverse read suffixes below:

#run basecalling job on N (SBATCH -ni) threads/cores

guppy_basecaller --device "cuda:0" \
	-i ${INDIR}/fast5/ \
        --save_path ${INDIR}/basecall/ \
	--flowcell ${FLOWCELL} \
	--kit SQK-LSK109 \
	--calib_detect \
	--num_callers 12 \
	--gpu_runners_per_device 32 \
	--compress_fastq \
	--post_out #\
#	--resume

guppy_barcoder --device "cuda:0" \
        -i ${INDIR}/basecall/ \
	--recursive \
        --save_path ${INDIR}/demultiplex/ \
	--allow_inferior_barcodes \
	--detect_mid_strand_barcodes \
	--trim_barcodes \
	--compress_fastq \
	--barcode_kits "EXP-NBD196 EXP-NBD104 EXP-NBD114 SQK-NBD110-24 SQK-NBD110-96"

if [ $? -eq 0 ]
then
	echo -e "\tThe results have been stored in '${WORKDIR}' in ${SLURMD_NODENAME}.\n\tCopying the results to ${INDIR} ..."
	cp -rf ${WORKDIR}/fast5_pass ${INDIR}/
	if [ $? -eq 0 ]
	then
		echo -e "\tCopying successful..."
		mkdir ${INDIR}/demultiplex_fast5
    		cd ${INDIR}/demultiplex_fast5
		ls ${INDIR}/demultiplex | grep "barcode" | xargs mkdir
		mkdir unclassified
		gunzip ${INDIR}/*/*/*fastq.gz
	elif [ $? -ne 0 ]
	then
		echo -e "\tCopying was NOT successful..."
	fi
fi
