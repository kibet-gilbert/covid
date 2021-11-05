#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -J Guppy_gpu
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

singularity exec --nv ${INDIR}/../../scripts/guppy/ont-guppy.5.0.7.sif guppy_basecaller --device "cuda:0" \
	-i ${INDIR}/fast5_skip/ \
        --save_path ${INDIR}/basecall/ \
	--config ${INDIR}/../../scripts/guppy/data/dna_r9.4.1_450bps_hac.cfg \
	--barcode_kits SQK-NBD110-96 \
	--num_callers 12 \
	--as_gpu_runners_per_device 32 \
	--compress_fastq \
	--post_out \
	--fast5_out \
	--resume 

singularity exec --nv ${INDIR}/../../scripts/guppy/ont-guppy.5.0.7.sif guppy_basecaller --device "cuda:0" \
        -i ${INDIR}/basecall/ \
        --save_path ${INDIR}/demultiplex/ \
	--require_barcodes_both_ends \
	--barcode_kits SQK-NBD110-96 \
	--arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg barcode_arrs_nb96.cfg"


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
