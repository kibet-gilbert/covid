#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -J nanoseq_nf
#SBATCH -n 18

set -uE
trap ' echo Error $? occured on $LINENO && exit 1' ERR

#load the nextflow module
module load nextflow/21.04.1
NXF_OPTS='-Xms1g -Xmx4g'
WKDIR=/var/scratch/$USER/work_nano
NXF_SINGULARITY_CACHEDIR=${WKDIR}/singularity
export NXF_SINGULARITY_CACHEDIR=${WKDIR}/singularity

#Setting work Directory in HPC
INDIR=$PWD
src_dir=${INDIR##*/}
DATADIR=${INDIR}/../core_data
WORKDIR=/var/scratch/$USER/"`date +"%Y.%m.%d"`_${src_dir}"
src_dir1=${WORKDIR##*/}
mkdir -p $WORKDIR
cd $WORKDIR

echo -e "\tUsing $WORKDIR on $SLURMD_NODENAME"
#echo -e "\tSetting up the working environment installations"
#rsync -u --size-only "${ENVDIR}" "${ENVDIR_WKDIR}"
source ${INDIR}/../../scripts/process_files.sh

#generating the samplesheet.csv file: Please provide the forward and reverse read suffixes below:

#run basecalling job on N (SBATCH -ni) threads/cores

nextflow run nf-core/nanoseq -r 1.1.0 \
	-w $WKDIR \
	--input ${INDIR}/samplesheet_bc_dx.csv \
	--protocol cDNA \
    	--input_path ${INDIR}/fast5_skip/ \
	--flowcell FLO-MIN106 \
	--barcode_kit SQK-NBD11096 \
	--kit SQK-LSK109 \
	--skip_quantification \
	--skip_differential_analysis \
	--max_cpus 18 \
	--max_memory '40.GB' \
	-profile singularity -resume |& tee ${INDIR}/nf_viralrecon.out.log
if [ $? -eq 0 ]
then
	echo -e "\tThe results have been stored in '${WORKDIR}' in ${SLURMD_NODENAME}.\n\tCopying the results to ${INDIR} ..."
	cp -rf ${WORKDIR}/results ${INDIR}/
	if [ $? -eq 0 ]
	then
		echo -e "\tCopying successful..."
	elif [ $? -ne 0 ]
	then
		echo -e "\tCopying was NOT successful..."
	fi
