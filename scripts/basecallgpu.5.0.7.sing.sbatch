#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -J guppygpu.5.0.7.sing
#SBATCH --qos vip
#SBATCH --gres=gpu:v100:1
#SBATCH --nodelist=compute06
#SBATCH -n 1

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
WORKDIR=/var/scratch/$USER/"`date +"%Y.%m.%d"`_${src_dir}"
src_dir1=${WORKDIR##*/}
until [[ -d ${WORKDIR} ]]
do
	echo "Creating output directory '${WORKDIR}'"
	mkdir -p ${WORKDIR}/fast5_pass
done
cd $WORKDIR

echo -e "\tUsing $WORKDIR/fast5_pass on $SLURMD_NODENAME"
#echo -e "\tSetting up the working environment installations"
#rsync -u --size-only "${ENVDIR}" "${ENVDIR_WKDIR}"
#source ${INDIR}/../../scripts/process_files.sh

#generating the samplesheet.csv file: Please provide the forward and reverse read suffixes below:

#run basecalling job on N (SBATCH -ni) threads/cores
singularity run --app guppy_basecaller --nv ${INDIR}/../../scripts/guppy/ont_guppy_GPU-5.0.7.sif \
	--device "cuda:0" \
	-i ${INDIR}/fast5_skip/ \
        --save_path ${INDIR}/basecall/ \
	--flowcell FLO-MIN106 \
	--kit SQK-NBD110-96 \
	--calib_detect \
	--gpu_runners_per_device 32 \
	--num_callers 12 \
	--compress_fastq \
	--post_out #\
#	--resume 

singularity run --app guppy_barcoder --nv ${INDIR}/../../scripts/guppy/ont_guppy_GPU-5.0.7.sif \
	--device "cuda:0" \
	--recursive \
        -i ${INDIR}/basecall/ \
        --save_path ${INDIR}/demultiplex/ \
	--require_barcodes_both_ends \
	--compress_fastq \
	--barcode_kits SQK-NBD110-96

if [ $? -eq 0 ]
then
	echo -e "\tThe results have been stored in '${WORKDIR}' in ${SLURMD_NODENAME}.\n\tCopying the results to ${INDIR} ..."
	cp -rf ${WORKDIR}/fast5_pass ${INDIR}/
	if [ $? -eq 0 ]
	then
		echo -e "\tCopying successful..."
	elif [ $? -ne 0 ]
	then
		echo -e "\tCopying was NOT successful..."
	fi
fi
