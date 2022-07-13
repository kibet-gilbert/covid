#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -J MetagenMag_nf
#SBATCH --qos vip
#SBATCH -n 30

set -uE
trap ' echo Error $? occured on $LINENO && exit 1' ERR

#load the nextflow module
module load nextflow/21.10
NXF_OPTS='-Xms1g -Xmx4g'
WKDIR=/var/scratch/$USER/work_mag
NXF_SINGULARITY_CACHEDIR=${WKDIR}/singularity
export NXF_SINGULARITY_CACHEDIR=${WKDIR}/singularity

#Setting work Directory in HPC
INDIR=$PWD
FASTQDIR=${INDIR}/fastq
src_dir=${INDIR##*/}
DATADIR=${INDIR}/../core_data
WORKDIR=/var/scratch/$USER/${src_dir}_analysis
src_dir1=${WORKDIR##*/}
mkdir -p $WORKDIR
cd $WORKDIR

echo -e "\tUsing $WORKDIR on $SLURMD_NODENAME"
#echo -e "\tSetting up the working environment installations"
#rsync -u --size-only "${ENVDIR}" "${ENVDIR_WKDIR}"
source ${INDIR}/../../scripts/process_files.sh

#generating the samplesheet.csv file: Please provide the forward and reverse read suffixes below: Allow for grouped samples or long-reads inclution.

#MiSeq
#FRead_suffix="L001_R1_001.fastq.gz"
#RRead_suffix="L001_R2_001.fastq.gz"

#NextSeq
#FRead_suffix="con_R1_001.fastq.gz"
#RRead_suffix="con_R2_001.fastq.gz"

#samplesheet_gen -f ${FRead_suffix} -r ${RRead_suffix}

#run variant calling job on N (SBATCH -n) threads/cores

# The actual data run : for the analyisis
# NOTES: Assuming that you have n CPUs allocated to the job here are some key points on CPU allocation/usage to consider
#	Limit the CPU availed to every process to between n-1(best) and 2/3(least) of n (--max_cpus)
#	If maximum memory available "mem" Gbs, limit it to half that i.e 1/2 mem (--max_memory)
#*_con_R{1,2}_001.fastq.gz

#nextflow pull nf-core/mag

nextflow run nf-core/mag -r 2.1.1 \
	--input $INDIR/fastq/samplesheet.csv \
	--email g.kibet@cgiar.org \
	-w $WKDIR \
	--single_end false \
	--outdir $WORKDIR/results \
	--igenomes_base 's3://ngi-igenomes/igenomes' \
	--megahit_fix_cpu_1 \
	--spades_fix_cpus 4 \
	--spadeshybrid_fix_cpus 4 \
	--metabat_rng_seed 1 \
	--save_trimmed_fail false \
	--fastp_qualified_quality 20 \
	--fastp_cut_mean_quality 20 \
	--host_genome GRCh37 \
	--host_removal_verysensitive false \
	--host_removal_save_ids \
	--skip_adapter_trimming false \
	--centrifuge_db $DATADIR/metagen/hpvc.tar.gz \
	--cat_db_generate false \
	--kraken2_db $DATADIR/metagen/k2_pluspf_20210517.tar.gz \
	--skip_krona false \
	--cat_db $DATADIR/metagen/CAT_prepare_20210107.tar.gz \
	--cat_db_generate false \
	--save_cat_db false \
	--gtdb $DATADIR/metagen/gtdbtk_r202_data.tar.gz \
	--gtdbtk_min_completeness 50 \
	--gtdbtk_max_contamination 10 \
	--gtdbtk_min_perc_aa 10 \
	--gtdbtk_min_af 0.65 \
	--gtdbtk_pplacer_cpus 4 \
	--gtdbtk_pplacer_scratch false \
	--coassemble_group false \
	--skip_spades false \
	--skip_spadeshybrid true \
	--skip_megahit false \
	--skip_quast false \
	--skip_prodigal false \
	--binning_map_mode 'group' \
	--skip_binning true \
	--min_contig_size 1500 \
	--min_length_unbinned_contigs 1000000 \
	--max_unbinned_contigs 100 \
	--skip_prokka false \
	--skip_busco false \
	--busco_download_path $DATADIR/metagen/ \
	--busco_auto_lineage_prok false \
	--save_busco_reference false \
	--max_cpus 20 \
	--max_memory '100.GB' \
	-profile singularity -resume |& tee ${INDIR}/nf_mag.out.log
if [ $? -eq 0 ]
then
	echo -e "\tThe results have been stored in '${WORKDIR}' in ${SLURMD_NODENAME}.\n\tCopying the results to ${INDIR} ..."
	cp -rf ${WORKDIR}/results ${INDIR}/
	if [ $? -eq 0 ]
	then
		echo -e "\tCopying successful..."
		results=results
		#Genome coverage
	elif [ $? -ne 0 ]
	then
		echo -e "\tCopying was NOT successful..."
	fi

fi
