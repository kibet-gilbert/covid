#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -J artic
#SBATCH --qos vip
#SBATCH -n 20

set -uE
trap ' echo Error $? occured on $LINENO && exit 1' ERR

#load the nextflow module
WKDIR=/var/scratch/$USER/work_nano

#Setting work Directory in HPC
INDIR=$PWD
src_dir=${INDIR##*/}
DATADIR=${INDIR}/../core_data
artic=${INDIR}/../../scripts/fieldbioinformatics/artic.sif
filter=${INDIR}/filter
variants=${INDIR}/variants
primer_schemes=${INDIR}/../../scripts/fieldbioinformatics/test-data/primer-schemes

#Read Filtering
ROOTFOLDER=${INDIR}/demultiplex
FOLDERS="barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07 barcode08 barcode09 barcode10 barcode11 barcode12 barcode13 barcode14 barcode15 barcode16 barcode17 barcode18 barcode19 barcode20 barcode21 barcode22 barcode23 barcode24 barcode25 barcode26 barcode27 barcode28 barcode29 barcode30 barcode31 barcode32 barcode33 barcode34 barcode35 barcode36 barcode37 barcode38 barcode39 barcode40 barcode41 barcode42 barcode43 barcode44 barcode45 barcode46 barcode47 barcode48 barcode49 barcode50 barcode51 barcode52 barcode53 barcode54 barcode55 barcode56 barcode57 barcode58 barcode59 barcode60 barcode61 barcode62 barcode63 barcode64 barcode65 barcode66 barcode67 barcode68 barcode69 barcode70 barcode71 barcode72 barcode73 barcode74 barcode75 barcode76 barcode77 barcode78 barcode79 barcode80 barcode81 barcode82 barcode83 barcode84 barcode85 barcode86 barcode87 barcode88 barcode89 barcode90 barcode91 barcode92 barcode93 barcode94 barcode95 barcode96 unclassified"
until [[ -d ${filter} ]]
do
	echo "Creating output directory '${filter}'"
	mkdir ${filter}
done

for folder in ${FOLDERS}
do
	echo "filtering ${folder}"
	singularity exec ${artic} artic guppyplex \
		--min-length 400 \
		--max-length 700 \
		--directory ${ROOTFOLDER}/${folder} \
		--prefix ${filter}/filter
done

#MinION Pipeline
#artic minion --medaka --normalise 200 --threads 4 --scheme-directory ~/artic-ncov2019/primer_schemes --read-file run_name_barcode01.fastq nCoV-2019/V1 samplename
until [[ -d ${variants} ]]
do
	echo "Creating output directory '${variants}'"
	mkdir ${variants}
done

for folder in ${FOLDERS}
do
	folder_dir=${variants}/${folder}
	until [[ -d ${folder_dir} ]]
	do
		echo "Creating output directory '${folder_dir}"
		mkdir -p ${folder_dir}
	done
	echo "variant calling ${folder}"
	cd ${folder_dir}
	singularity exec ${artic} artic minion \
		--medaka \
		--normalise 200 \
		--threads 20 \
		--scheme-directory ${primer_schemes} \
		--read-file ${filter}/filter_${folder}.fastq \
		nCoV-2019/V3 ${folder}
done 

#awk 'BEGIN{}; /^Creating output/{split($0,b,"/"); filen=b[10]".log"}; /^variant calling/{split($0,a," "); filen=a[3]".log"}; {print $0 > filen}' slurm-721007.out
