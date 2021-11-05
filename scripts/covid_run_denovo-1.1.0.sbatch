#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -J VarCall-denovo_nf
#SBATCH --qos vip
#SBATCH -n 40

set -uE
trap ' echo Error $? occured on $LINENO && exit 1' ERR

#load the nextflow module
module load nextflow/20.10.0
module load nextclade/0.12.0
NXF_OPTS='-Xms1g -Xmx4g'
NXF_SINGULARITY_CACHEDIR=/var/scratch/gkibet/work/singularity

#Setting work Directory in HPC
INDIR=$PWD
src_dir=${INDIR##*/}
#Testing the current working dir
#regex="^/home/${USER}/*/covid/data/${src_dir}$"
#if [[ $INDIR =~ $regex ]]
#then
#	continue
#else
#	echo -e "Current working directory is not in the right location/path. Please ensure that your working/input data directory:($PWD) is in /home/${USER}/*/covid/data/<working directory>"
#	exit 1
#fi

DATADIR=${INDIR}/../core_data
ENVDIR=${INDIR}/../../work/conda/nf-core-viralrecon-1.1.0-5283800d5f631c4fbcac00bc34786af1
WORKDIR=/var/scratch/$USER/"${src_dir}-`date +"%d.%m.%Y"`"
mkdir -p $WORKDIR
WKDIR=/var/scratch/$USER/work
ENVDIR_WKDIR=${WKDIR}/conda/nf-core-viralrecon-1.1.0-5283800d5f631c4fbcac00bc34786af1
mkdir -p $ENVDIR_WKDIR
cd $WORKDIR

echo -e "\tUsing $WORKDIR on $SLURMD_NODENAME"
#echo -e "\tSetting up the working environment installations"
#rsync -u --size-only "${ENVDIR}/*" "${ENVDIR_WKDIR}"
source ${INDIR}/../../scripts/process_files.sh

#generating the samplesheet.csv file: Please provide the forward and reverse read suffixes below:
FRead_suffix="L001_R1_001.fastq.gz"
RRead_suffix="L001_R2_001.fastq.gz"
samplesheet_gen -f ${FRead_suffix} -r ${RRead_suffix}

#run variant calling job on N (SBATCH -n) threads/cores

# The actual data run : for the analyisis
# NOTES: Assuming that you have n CPUs allocated to the job here are some key points on CPU allocation/usage to consider
#       Limit the CPU availed to every process to at least 1/2 of n (--max_cpus)
#       If maximum memory available "mem" Gbs, it is best to limit it to half that i.e 24GB (--max_memory)
#

nextflow run nf-core/viralrecon -r 1.1.0 \
	-w $WKDIR \
	--input ${INDIR}/samplesheet.csv \
	--cut_mean_quality 20 \
	--qualified_quality_phred 20 \
	--min_trim_length 50 \
	--unqualified_percent_limit 10 \
	--mpileup_depth 0 \
	--min_base_qual 20 \
	--ivar_trim_min_len 20 \
	--ivar_trim_min_qual 20 \
	--ivar_trim_window_width 4 \
	--min_allele_freq 0.25 \
	--max_allele_freq 0.75 \
	--amplicon_left_suffix '_LEFT' \
	--amplicon_right_suffix '_RIGHT' \
	--fasta ${DATADIR}/GCA_009858895.3_ASM985889v3_genomic.200409.fna.gz \
	--gff ${DATADIR}/GCA_009858895.3_ASM985889v3_genomic.200409.gff.gz \
	--protocol amplicon \
	--amplicon_bed ${DATADIR}/nCoV-2019.artic.V3.bed \
	--amplicon_fasta ${DATADIR}/nCoV-2019.artic.V3.primer.fasta \
	--assemblers "spades" \
	--skip_blast \
	--callers ivar \
	--max_cpus 40 \
	--max_memory '40.GB' \
	-profile singularity -resume
if [ $? -eq 0 ]
then
	echo -e "\tThe results have been stored in '${WORKDIR}' in ${SLURMD_NODENAME}.\n\tCopying the results to ${INDIR} ..."
	mkdir -p ${INDIR}/results_denovo
	cp -rf ${WORKDIR}/results/* ${INDIR}/results_denovo
		if [ $? -eq 0 ]
		then
			echo -e "\tCopying successful..."
			echo -e "\tGenerating variant annotation (*.snpEff.vcf.tsv) files"
			results=${INDIR}/results_denovo
			variant_calls
		elif [ $? -ne 0 ]
		then
			echo -e "\tCopying was NOT successful..."
		fi
	echo -e "Performing nextclade analysis"
	#cat results/variants/ivar/consensus/*.fa >> results/variants/ivar/consensus/covid_11-02-2021_con.fa
	src_dir1=${WORKDIR##*/}
	touch ${INDIR}/${results}/variants/ivar/consensus/${src_dir1}_con.fa
	cat ${INDIR}/${results}/variants/ivar/consensus/*.consensus.fa ${DATADIR}/sar-cov-2-MN908947.3.fa >> ${INDIR}/${results}/variants/ivar/consensus/${src_dir1}_con.fa
	nextclade.js -i ${INDIR}/${results}/variants/ivar/consensus/${src_dir1}_con.fa --output-tree ${src_dir1}_con.json #--input-tree --input-root-seq ${DATADIR}/sar-cov-2-MN908947.3.fa --input-gene-map

fi


