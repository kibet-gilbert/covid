#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -J VarCall_nf
#SBATCH --qos vip
#SBATCH -n 10

set -uE
trap ' echo Error $? occured on $LINENO && exit 1' ERR

#load the nextflow module
module load iqtree/2.2.0
WKDIR=/var/scratch/$USER/work

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
source ${INDIR}/../../../scripts/process_files.sh

#Nextclade analysis:
#singularity run ${INDIR}/../../scripts/nextclade/nextclade_latest.sif nextclade run -i ${INDIR}/alignment/${src_dir1}.consensus.all.fasta --output-dir ${INDIR}/nextclade --output-basename ${src_dir1}.nextclade --output-tree ${INDIR}/nextclade/${src_dir1}.nextclade.json --output-tsv ${INDIR}/nextclade/${src_dir1}.nextclade.tsv -r ${DATADIR}/sar-cov-2-MN908947.3.fa -a ${DATADIR}/tree.json --input-qc-config ${DATADIR}/qc.json --input-gene-map ${DATADIR}/genemap.gff |& tee -a ${INDIR}/nextclade/${src_dir1}.nextclade.log

#Phylohenetic inference
if [ ! -d ${INDIR}/phylogeny ]; then mkdir -p ${INDIR}/phylogeny; fi
iqtree2 -s ${INDIR}/DRC-74_global-399_open_metadata.filtered.fasta -m GTR+F+R3 --date ${INDIR}/DRC-74_global-399_open_metadata.filtered.dates.tsv -T AUTO -te ${INDIR}/phylogeny/DRC-74_global-399_open_metadata.filtered.nextclade.aligned.ufboot.treefile -pre ${INDIR}/phylogeny/DRC-74_global-399_open_metadata.filtered.nextclade.aligned.ufboot.timetree --date-ci 100 
