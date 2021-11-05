#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -J VarCallTest_nf
#SBATCH -n 7

set -uE
trap ' echo Error $? occured on $LINENO && exit 1' ERR

#load the nextflow module
module load nextflow/19.10.0
NXF_OPTS='-Xms1g -Xmx4g'


#Setting work Directory in HPCi
INDIR=$PWD
src_dir=${INDIR##*/}

ENVDIR=${INDIR}/../../work/conda/nf-core-viralrecon-1.1.0-5283800d5f631c4fbcac00bc34786af1
WORKDIR=/var/scratch/$USER/"${src_dir}-`date +"%d.%m.%Y"`"
mkdir -p $WORKDIR
WKDIR=/var/scratch/$USER/work/
mkdir -p $WKDIR
cd $WORKDIR
echo -e "Setting up conda environment"
conda env create --prefix /var/scratch/${USER}/work/conda/nf-core-viralrecon-1.1.0-5283800d5f631c4fbcac00bc34786af1 --file /home/${USER}/bioinformatics/github/covid/scripts/covid_env.yml

echo -e "\tUsing $WORKDIR on $SLURMD_NODENAME"

#run variant calling job on N (SBATCH -n) threads/cores
# 1. Test code : will help set-up/install the tools and run a test run.
# Uses: -profile conda
nextflow run nf-core/viralrecon -r 1.1.0 -w $WKDIR --input ${INDIR}/samplesheet_test_amplicon.csv --max_memory '40.GB' --max_cpus 6 -profile test,conda -resume

# Uses: -profile singularity
#nextflow run nf-core/viralrecon -r 1.1.0 -w $WKDIR -profile test,singularity -resume

