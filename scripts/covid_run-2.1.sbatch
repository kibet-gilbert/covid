#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -J VarCall_nf
#SBATCH --qos vip
#SBATCH -n 40

set -uE
trap ' echo Error $? occured on $LINENO && exit 1' ERR

#load the nextflow module
module load nextflow/21.04.1
module load nextclade/0.12.0
module load muscle/3.8.31
module load iqtree/1.6.12
NXF_OPTS='-Xms1g -Xmx4g'
WKDIR=/var/scratch/$USER/work
NXF_SINGULARITY_CACHEDIR=${WKDIR}/singularity
export NXF_SINGULARITY_CACHEDIR=${WKDIR}/singularity

#Setting work Directory in HPC
INDIR=$PWD
FASTQDIR=${INDIR}/fastq
#FASTQDIR=${INDIR}
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
FRead_suffix="con_R1_001.fastq.gz"
RRead_suffix="con_R2_001.fastq.gz"
#samplesheet_gen -f ${FRead_suffix} -r ${RRead_suffix}

#run variant calling job on N (SBATCH -n) threads/cores

# The actual data run : for the analyisis
# NOTES: Assuming that you have n CPUs allocated to the job here are some key points on CPU allocation/usage to consider
#	Limit the CPU availed to every process to between n-1(best) and 2/3(least) of n (--max_cpus)
#	If maximum memory available "mem" Gbs, limit it to half that i.e 1/2 mem (--max_memory)
#

nextflow run nf-core/viralrecon -r 2.1 \
	-w $WKDIR \
	--input ${FASTQDIR}/samplesheet.csv \
	--platform illumina \
	--skip_markduplicates true \
	--filter_duplicates false \
	--min_mapped_reads 1000 \
        --primer_left_suffix '_LEFT' \
        --primer_right_suffix '_RIGHT' \
    	--fasta ${DATADIR}/GCA_009858895.3_ASM985889v3_genomic.200409.fna.gz\
	--gff ${DATADIR}/GCA_009858895.3_ASM985889v3_genomic.200409.gff.gz \
	--protocol amplicon \
	--primer_fasta ${DATADIR}/nCoV-2019.artic.V3.primer.fasta \
	--primer_bed ${DATADIR}/nCoV-2019.artic.V3.bed \
	--primer_set_version 3 \
	--callers ivar \
	--kraken2_db_name 'human' \
	--skip_assembly \
	--skip_abacas \
	--skip_plasmidid \
	--skip_cutadapt \
	--skip_assembly_quast \
	--max_cpus 40 \
	--max_memory '40.GB' \
	-profile singularity -resume |& tee ${INDIR}/nf_viralrecon.out.log
if [ $? -eq 0 ]
then
	echo -e "\tThe results have been stored in '${WORKDIR}' in ${SLURMD_NODENAME}.\n\tCopying the results to ${INDIR} ..."
	cp -rf ${WORKDIR}/results ${INDIR}/
	if [ $? -eq 0 ]
	then
		echo -e "\tCopying successful..."
		echo -e "\tGenerating variant annotation (*.snpEff.vcf.tsv) files"
		results=results
		snpeff_dir="variants/ivar/snpeff"
		#Variants
		variant_calls
		echo -e "\tCopying amplicon/genome coverage plots into ${INDIR}/plots/"
		#Genome coverage
		mkdir -p $INDIR/plots/amplicon $INDIR/plots/genome
		cp ${INDIR}/${results}/variants/bowtie2/mosdepth/amplicon/*.pdf ${INDIR}/plots/amplicon/
		cp ${INDIR}/${results}/variants/bowtie2/mosdepth/genome/*.pdf ${INDIR}/plots/genome/
		convert ${INDIR}/plots/amplicon/*.pdf $INDIR/plots/genome/*.pdf ${INDIR}/plots/${src_dir1}.amplicon.genome.coverage.pdf
	elif [ $? -ne 0 ]
	then
		echo -e "\tCopying was NOT successful..."
	fi

	#Post Variant-Calling analysis analysis
	if [ -f ${INDIR}/${results}/variants/ivar/consensus/${src_dir1}_con.fa ]; then rm ${INDIR}/${results}/variants/ivar/consensus/${src_dir1}_con.fa; fi &&	cat ${INDIR}/${results}/variants/ivar/consensus/*.consensus.fa >> ${INDIR}/${results}/variants/ivar/consensus/${src_dir1}_con.fa
	if [ ! -d ${INDIR}/alignment ]; then mkdir -p ${INDIR}/alignment; fi
	if [ -f ${INDIR}/alignment/${src_dir1}_con.fasta ]; then rm ${INDIR}/alignment/${src_dir1}_con.fasta; fi && cat ${INDIR}/${results}/variants/ivar/consensus/${src_dir1}_con.fa ${DATADIR}/sar-cov-2-MN908947.3.fa >> ${INDIR}/alignment/${src_dir1}_con.fasta


	#Nextclade analysis
	echo -e "Performing nextclade analysis"
	if [ ! -d ${INDIR}/nextclade ]; then mkdir -p ${INDIR}/nextclade; fi
	echo -e "\tBeginning nextclade analysis....\n\tDownloading 'nextclade datasets': genemap.gff, qc.json, sequences.fasta, tree.json, primers.csv, reference.fasta, tag.json"
	singularity run ./nextclade_latest.sif nextclade dataset get --name='sars-cov-2' --reference='MN908947' --output-dir='../../data/core_data'
	echo -e "\tProceeding with Nextclade analysis."
	singularity run ${INDIR}/../../scripts/nextclade/nextclade_latest.sif nextclade run -i ${INDIR}/${results}/variants/ivar/consensus/${src_dir1}_con.fa --output-dir ${INDIR}/nextclade --output-basename ${src_dir1}.nextclade --output-tree ${INDIR}/nextclade/${src_dir1}.nextclade.json --output-tsv ${INDIR}/nextclade/${src_dir1}.nextclade.tsv -r ${DATADIR}/sar-cov-2-MN908947.3.fa -a ${DATADIR}/tree.json --input-qc-config ${DATADIR}/qc.json --input-gene-map ${DATADIR}/genemap.gff |& tee -a ${INDIR}/nextclade/${src_dir1}.nextclade.log
	#Nextclade V 0.12.0 HPC
	#nextclade.js -i ${INDIR}/${results}/variants/ivar/consensus/${src_dir1}_con.fa --output-tree ${INDIR}/nextclade/${src_dir1}_con.json --output-tsv ${INDIR}/nextclade/${src_dir1}_con.tsv #--input-tree --input-root-seq ${DATADIR}/sar-cov-2-MN908947.3.fa --input-gene-map
	
	#Pangolin analysis
	if [ ! -d ${INDIR}/pangolin ]; then mkdir -p ${INDIR}/pangolin; fi
	#echo -e "Performing Pangolin analysis"

	#Multiple Sequence Alignment
	if [ -f ${INDIR}/alignment/${src_dir1}_con.aln ]; then rm ${INDIR}/alignment/${src_dir1}_con.aln; fi && cat ${INDIR}/nextclade/${src_dir1}.nextclade.aligned.fasta ${DATADIR}/sar-cov-2-MN908947.3.fa >> ${INDIR}/alignment/${src_dir1}_con.aln
	#muscle -in ${INDIR}/alignment/${src_dir1}_con.fasta -fastaout ${INDIR}/alignment/${src_dir1}_con.aln -maxiters 2

	#Phylohenetic inference
	if [ ! -d ${INDIR}/phylogeny ]; then mkdir -p ${INDIR}/phylogeny; fi
	iqtree -s ${INDIR}/alignment/${src_dir1}_con.aln -m MFP -bb 1000 -bnni -alrt 1000 -nt AUTO -pre ${INDIR}/phylogeny/${src_dir1}_con
fi
