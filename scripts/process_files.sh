#!/usr/bin/env bash

AWK_EXEC=$( which gawk )
PYTHON_EXEC=$( which python )

realpath() { #
	${PYTHON_EXEC} -c "import os,sys; print(os.path.realpath(sys.argv[1]))" $1
}

if [ -z ${INDIR} ]
then
	INDIR=$PWD
	echo -e "\tInput directory 'INDIR' has been set to ${INDIR}"
fi
deml=${INDIR}/../../scripts/deML/src/deML
#set -u

#INDIR=$(pwd)
rename() { #Takes input file name and extracts the prefix of the name ("prefix"), The suffix (".suffix"), Absolute path to the file and the name of the Directory within which the file is found.
	input_filename=`basename -- "$@"`
	output_filename=${input_filename%.*}
	filename_ext=${input_filename##*.}
	full_path=$(realpath "$@")
	src_dir_path=`dirname ${full_path}`
	src_dir=${src_dir_path##*/}
	#echo -e "input is $i \ninput_filename is $input_filename \noutput_filename is $output_filename \nfilename_ext is $filename_ext \nsrc_dir_path is $src_dir_path \nsrc_dir is $src_dir"
}

#===============================================================================================================================================================
samplesheet_gen() { #Generate a samplesheet.csv
	#Error trapping
	trap ' echo Error $? occured on $LINENO && return 1' ERR
	if [ $# -ne 4 ]
	then
		echo "Input error..."
		echo "Usage: ${FUNCNAME[0]} -f <R1.fastq.gz> -r <R2.fastq.gz>"
		return 1
	fi

	#Parsing arguments provided to the samplesheet_gen function
        local OPTIND=1
        while getopts 'f:r:' key
        do
                case "${key}" in
                        f)
				if [[ ! $OPTARG =~ .*(\.fastq\.gz)$ ]]
				then
					echo -e "\tinput error: $OPTARG is not a vialble fastq.gz file extension!\nExample: ${FUNCNAME[0]} -f <R1.fastq.gz> -r <R2.fastq.gz> ::Replace <Rn.fastq.gz>"
					return 1
				elif [[ ( "$OPTARG" =~ .*(\.fastq\.gz)$ ) && ( "$(ls -A $FASTQDIR/*$OPTARG)" ) ]]
				then
					f_ext=$OPTARG
					fastqgz_suf=$f_ext
					#Generating an array containing all forward read .fastq.gz files in the working dir. based on the provided suffix
					array=($FASTQDIR/*$fastqgz_suf)
					f_array=("${array[@]}")
				else
					echo -e "\tinput error: $OPTARG is not a vialble fastq.gz file extension!"
					echo "Example: ${FUNCNAME[0]} -f <R1.fastq.gz> -d <R2.fastq.gz> ::Replace <Rn.fastq.gz>"
					return 1
				fi
				;;
			r)
				if [[ ( ! $OPTARG =~ .*(\.fastq\.gz)$ ) || ( ${f_ext} =~ ${OPTARG} ) ]]
				then
					echo -e "\tinput error: $OPTARG is not a vialble fastq.gz file extension!"
					echo "Example: ${FUNCNAME[0]} -f <R1.fastq.gz> -d <R2.fastq.gz> ::Replace <Rn.fastq.gz>"
					return 1
				elif [[ ( $OPTARG =~ .*(\.fastq\.gz)$ ) && ( "$(ls "$FASTQDIR"/*"$OPTARG")" ) ]]
                                then
					r_ext=$OPTARG
                                        fastqgz_suf=$r_ext
					#Generating an array containing all reverse read *fastq.gz files from the working dir. based on the provided suffix
					array=($FASTQDIR/*$fastqgz_suf)
					r_array=("${array[@]}")
				else
					echo -e "\tinput error: $OPTARG is not a vialble fastq.gz file extension!"
					echo "Example: ${FUNCNAME[0]} -f <R1.fastq.gz> -d <R2.fastq.gz> ::Replace <Rn.fastq.gz>"
					return 1
                                fi
                                ;;
                        ?)
				echo "Forward/Reverse read file suffixes have not been provided...i.e *R?(1|2)_?[0-9]*.fastq.qz "
                                echo "Example: ${FUNCNAME[0]} -f <R1.fastq.gz> -d <R2.fastq.gz> ::Replace <Rn.fastq.gz>"
                                return 1
                                ;;
                esac
        done

	index=0
	echo "sample,fastq_1,fastq_2" > "$FASTQDIR/samplesheet.csv"
	#echo "sample,group,short_reads_1,short_reads_2,long_reads" > "$FASTQDIR/samplesheet.csv"
	for file in ${f_array[@]}
	do
		ffilename=`basename -- ${file}`
		#Testing if the provided forward and reverse read suffices are similar.
		if [[ ${f_ext} =~ ${r_ext} ]]
		then
			echo "ERROR!: forward and reverse read extensions provided are similar"
			return 1
		fi

		#All conditions passed, we proceed to generate the samplesheet
		f_prefix=${ffilename%_${f_ext}}
		for file2 in ${r_array[@]}
		do
			rfilename=`basename -- ${file2}`
			r_prefix=${rfilename%_${r_ext}}
			if [[ "$f_prefix" =~ "$r_prefix" ]]
			then

				echo "${f_prefix}",$file,$file2 >> "$FASTQDIR/samplesheet.csv"
				#echo "${f_prefix}",,$file,$file2, >> "$FASTQDIR/samplesheet.csv"
				index=`expr ${index} + 1`
			else
				continue
			fi
		done
	done
	if [ $? -eq 0 ]
	then
		echo -e "\t${index} Samples have been added to the samplesheet.csv for processing"
	fi
}

#===============================================================================================================================================================
variant_calls() { #This file compiles a list of variant annotations in a .tsv file.
	#Error trapping
	trap ' echo Error $? occured on $LINENO && return 1' ERR
	if [ -z ${results} ]
	then
		results=results
		echo -e "\tresults* directory 'results' has been set to './${results}'"
	fi

	if [ -z ${snpeff_dir} ]
	then
		snpeff_dir="variants/ivar/snpeff"
		echo -e "\tresults* directory 'results' has been set to './${results}'"
	fi
	variants_dir=$INDIR/${results}/${snpeff_dir}
	snpeff_suf=*.snp?ff.vcf.gz
	src_dir=${INDIR##*/}
	#Generating an array containing all annotated variant call *.snpEff,vcf.gz files in the results/*/ivar/sneff dir.
	vcf_array=($variants_dir/$snpeff_suf)
	mutations=$INDIR/mutations
	mkdir -p ${mutations}
	variantstsv=$INDIR/${src_dir}.mutations.tsv
	var=$INDIR/${src_dir}_aa.mutations.tsv
	aa_codes=$INDIR/../../scripts/aa_codes

	echo -e "Sample\tMutation-Annotated[gene:(filter:frequency)]" > $variantstsv
	for i in "${vcf_array[@]}"
	do
		if [ ! -f $i ]
		then
			echo "input error: file '$i' is non-existent!"
			return 1
		elif [[ ( `basename -- "$i"` =~ .*\.AF0\.75\.* ) ]]
		then
			#skipping variant calls given based on 0.75 frequency cut-off. nf-core/viralrecon pipeline generates it automattically.
			continue 1
		elif [[ ( -f $i ) && ( `basename -- "$i"` =~ .*\.snp.ff\.vcf\.gz$ ) ]]
		then
			rename $i
			#Unzipping the *.snpEff.vcf.gz files
			gunzip -c $i > ${src_dir_path}/${output_filename}
			input=${src_dir_path}/${output_filename}
			outfile=`echo "${mutations}/${output_filename}.tsv"`
			#Awk code processing the unziped variant call files to simpler tsv test files.
			${AWK_EXEC} -v outfile=$outfile -v variants=$variantstsv -v sample=$output_filename 'BEGIN{FS="\t";
			OFS="\t"} {
			if (FNR==NR) {y[$1]=$2}
			else {
				if ($0 ~ /^##/) {next}
				else if ($0 ~ /^#CHROM/) {print $1,$2,$4,$5,$7,"Total_Depth","Annotation","Annotation_Impact","Gene_Name","Rank"," HGVS.c","HGVS.p","cDNA.pos/cDNA.length","AA.pos/AA.length","GT:REF_DP:REF_RV:REF_QUAL:ALT_DP:ALT_RV:ALT_QUAL:ALT_FREQ" > outfile } 
				else if ($0 ~ !/^#/) {split($8,a,"|");b=a[1];
					split(b,c,"[=;]"); dp=c[2]; ann=a[2];
					anni=a[3]; gn=a[4]; rk=a[8]; 
					dnamute=a[10]; aamute=a[11];
					dnapos=a[12]; aapos=a[14]; {
					for (key in y) {
						gsub(key,y[key],aamute)
						gsub(/p\./,"",aamute)}};
					split($10,z,":");f=z[8];
					mut[$2]=gn":"aamute"("$7";"f")";
					gsub(".snp.ff.vcf","\t",sample);
					print $1,$2,$4,$5,$7,dp,ann,anni,gn,rk,dnamute,aamute,dnapos,aapos,$10 >> outfile }}};
					END {printf sample >> variants; for (site in mut) {printf "%s ", mut[site] >> variants}; printf "\n" >> variants }' ${aa_codes} ${input}
			
			${AWK_EXEC} -v sample=$output_filename '{
			gsub(".snp.ff.vcf","\t",sample)};
			!/^[[:space:]]*$/{gsub(" "sample," ",$0); print}' ${variantstsv} > ${variantstsv}.file && mv ${variantstsv}.file ${variantstsv}

			if [ $? -eq 0 ]
			then
				echo -e "\tDONE. Variant calls from ${output_filename}.gz are stored in ${outfile}"
			fi
		else
			echo "input file error in `basename $i`: input file should be a .snp[e|E]ff.vcf.gz file format"
			continue
		fi
	done

	echo -e "Sample\torf1ab:[Var(filter:freq)]\tS:[Var(filter:freq)]\tORF3a:[Var(filter:freq)]\tM:[Var(filter:freq)]\tN:[Var(filter:freq)]\tE:[Var(filter:freq)]\tORF6:[Var(filter:freq)]\tORF7a:[Var(filter:freq)]\tORF7b:[Var(filter:freq)]\tORF8:[Var(filter:freq)]\tORF9b:[Var(filter:freq)]\tORF10:[Var(filter:freq)]\tother:[gene:Var(filter:freq)]" > $var
	${AWK_EXEC} -v sample=$output_filename -v var=$var 'BEGIN{FS="\t";
	gene[1]="orf1ab"; gene[2]="S"; gene[3]="ORF3a";
	gene[4]="M"; gene[5]="N"; gene[6]="E"; gene[7]="ORF6";
	gene[8]="ORF7a"; gene[9]="ORF7b"; gene[10]="ORF8";
	gene[11]="ORF9b"; gene[12]="ORF10"; gene[13]="others"}
	FNR>1{delete orf1abA; delete SA; delete ORF3aA; delete MA;
	delete NA; delete EA; delete ORF6A; delete ORF7aA; delete ORF7bA;
	delete ORF8A; delete ORF9bA; delete ORF10A; delete othersA;
	delete aa; delete m;
	{split($2,aa," ")};{ for (mut in aa) { split(aa[mut],m,":");
		if (m[1] ~ /^orf1ab$/) {orf1abA[++n]=m[2]}
		else if (m[1] ~ "S") {SA[++n]=m[2]}
		else if (m[1] ~ "ORF3a") {ORF3aA[++n]=m[2]}
		else if (m[1] ~ "M") {MA[++n]=m[2]}
		else if (m[1] ~ "N") {NA[++n]=m[2]}
		else if (m[1] ~ "E") {EA[++n]=m[2]}
		else if (m[1] ~ "ORF6") {ORF6A[++n]=m[2]}
		else if (m[1] ~ "ORF7a") {ORF7aA[++n]=m[2]}
		else if (m[1] ~ "ORF7b") {ORF7bA[++n]=m[2]}
		else if (m[1] ~ "ORF8") {ORF8A[++n]=m[2]}
		else if (m[1] ~ "ORF9b") {ORF9bA[++n]=m[2]}
		else if (m[1] ~ "ORF10") {ORF10A[++n]=m[2]}
		else {othersA[++n]=m[1]":"m[2]}}};
	{printf $1"\t" >> var; for (g=1; g<=13; g++) {gn=gene[g];printf gn":" >> var;
		if (gn ~ "orf1ab") {
			for (i in orf1abA) {printf orf1abA[i]"," >> var} ; printf "\t" >> var}
		else if (gn ~ "S"){
			for (i in SA) {printf SA[i]"," >> var} ; printf "\t" >> var}
		else if (gn ~ "ORF3a"){
			for (i in ORF3aA) {printf ORF3aA[i]"," >> var} ; printf "\t" >> var}
		else if (gn ~ "M"){
			for (i in MA) {printf MA[i]"," >> var} ; printf "\t" >> var}
		else if (gn ~ "N"){
			for (i in NA) {printf NA[i]"," >> var} ; printf "\t" >> var}
		else if (gn ~ "E"){
			for (i in EA) {printf EA[i]"," >> var} ; printf "\t" >> var}
		else if (gn ~ "ORF6"){
			for (i in ORF6A) {printf ORF6A[i]"," >> var} ; printf "\t" >> var}
		else if (gn ~ "ORF7a"){
			for (i in ORF7aA) {printf ORF7aA[i]"," >> var} ; printf "\t" >> var}
		else if (gn ~ "ORF7b"){
			for (i in ORF7bA) {printf ORF7bA[i]"," >> var} ; printf "\t" >> var}
		else if (gn ~ "ORF8"){
			for (i in ORF8A) {printf ORF8A[i]"," >> var} ; printf "\t" >> var}
		else if (gn ~ "ORF9b"){
			for (i in ORF9bA) {printf ORF9bA[i]"," >> var} ; printf "\t" >> var}
		else if (gn ~ "ORF10"){
			for (i in ORF10A) {printf ORF10A[i]"," >> var} ; printf "\t" >> var}
		else if (gn ~ "others"){
			for (i in othersA) {printf othersA[i]"," >> var}} 
		}; printf "\n" >> var}}' ${variantstsv}
	if [ $? -eq 0 ]
	then
		echo -e "Cleaned Variant calls of amino acids are stored in ${var}"
	fi
}

#===============================================================================================================================================================
demultiplex_ML() { # This function uses deML(https://github.com/grenaud/deML) a machine learning algorithm to demultiplex a set of fastq.gz illumina output files.

	# Generate the index file:
	demldir=${INDIR}/demultiplex
	indexfile=${demldir}/index.txt

	#Error trapping
	trap ' echo Error $? occured on $LINENO ' ERR
	if [[ ($# -lt 2 && $1 =~ ^-(a|d)$) || $# -gt 8 ]]
	then
		:
		#echo "'$1'"
	else
		echo "Input error..."
		echo "Usage: ${FUNCNAME[0]} [-a -d] -f <FRead.R1.fastq.gz> -r <RRead.R2.fastq.gz> -x <FRead_index.R1.fastq.gz> -y <RRead_index.R2.fastq.gz> -s <SampleSheet.csv>"
		return 1
	fi
	
	R1regexp='^.*_L00._R1_001\.fastq\.gz$'
	R2regexp='^.*_L00._R2_001\.fastq\.gz$'
	I1regexp='^.*_L00._L1_001\.fastq\.gz$'
	I2regexp='^.*_L00._L2_001\.fastq\.gz$'
	SSregexp='^SampleSheet\.csv$'
	local OPTIND=1
	unset FRead
	unset RRead
	unset FRead_index
	unset RRead_index
	unset samplesheet

	while getopts 'adf:r:x:y:s:' key
	do
    		case "${key}" in
			a)
				echo -e "\tConcatenating all fastq.gz files."
				rm ${INDIR}/Undetermined_S0_L00a_R1_001.fastq.gz
				rm ${INDIR}/Undetermined_S0_L00a_R2_001.fastq.gz
				rm ${INDIR}/Undetermined_S0_L00a_I1_001.fastq.gz
				rm ${INDIR}/Undetermined_S0_L00a_I2_001.fastq.gz 
				cat ${INDIR}/*_L00?_R1_001.fastq.gz >> ${INDIR}/Undetermined_S0_L00a_R1_001.fastq.gz
				cat ${INDIR}/*_L00?_R2_001.fastq.gz >> ${INDIR}/Undetermined_S0_L00a_R2_001.fastq.gz
				cat ${INDIR}/*_L00?_I1_001.fastq.gz >> ${INDIR}/Undetermined_S0_L00a_I1_001.fastq.gz
				cat ${INDIR}/*_L00?_I2_001.fastq.gz >> ${INDIR}/Undetermined_S0_L00a_I2_001.fastq.gz
				FRead=${INDIR}/Undetermined_S0_L00a_R1_001.fastq.gz
				RRead=${INDIR}/Undetermined_S0_L00a_R2_001.fastq.gz
				FRead_index=${INDIR}/Undetermined_S0_L00a_I1_001.fastq.gz
				RRead_index=${INDIR}/Undetermined_S0_L00a_I2_001.fastq.gz
				samplesheet=${INDIR}/SampleSheet.csv

				if [[ -f $FRead && -f $RRead && -f $FRead_index && -f $RRead_index && -f $samplesheet ]]
				then
					echo -e "\tDefaults will be used::\n\t\tForward read: ${FRead}.\n\t\tReverse reads: ${RRead}.\n\t\tForward read indices: ${FRead_index}.\n\t\tReverse read indices: ${RRead_index}.\n\t\tSample sheet: ${samplesheet}."
				else
					echo "Input error..."
					echo "Usage: ${FUNCNAME[0]} [-a -d] -f <FRead.R1.fastq.gz> -r <RRead.R2.fastq.gz> -x <FRead_index.R1.fastq.gz> -y <RRead_index.R2.fastq.gz> -s <SampleSheet.csv>"
					return 1
				fi
				;;
			d)
				FRead=${INDIR}/Undetermined_S0_L001_R1_001.fastq.gz
				RRead=${INDIR}/Undetermined_S0_L001_R2_001.fastq.gz
				FRead_index=${INDIR}/Undetermined_S0_L001_I1_001.fastq.gz
				RRead_index=${INDIR}/Undetermined_S0_L001_I2_001.fastq.gz
				samplesheet=${INDIR}/SampleSheet.csv

				if [[ -f $FRead && -f $RRead && -f $FRead_index && -f $RRead_index && -f $samplesheet ]]
				then
					echo -e "\tDefaults will be used::\n\t\tForward read: ${FRead}.\n\t\tReverse reads: ${RRead}.\n\t\tForward read indices: ${FRead_index}.\n\t\tReverse read indices: ${RRead_index}.\n\t\tSample sheet: ${samplesheet}."
				else
					echo "Input error..."
					echo "Usage: ${FUNCNAME[0]} [-a -d] -f <FRead.R1.fastq.gz> -r <RRead.R2.fastq.gz> -x <FRead_index.R1.fastq.gz> -y <RRead_index.R2.fastq.gz> -s <SampleSheet.csv>"
					return 1
				fi
				;;
			f)
				if [[ ( -f $OPTARG ) && ( "$OPTARG" =~ $R1regexp ) ]]
				then
					FRead=${OPTARG}
				else
					until [[ "$FRead" =~ $R1regexp ]]
					do
						read -p "Please enter the absolute path to the fastq.gz file for the forward reads:: " FRead
					done
				fi
				;;
			r)
				if [[ ( -f $OPTARG ) && ( "$OPTARG" =~ $R2regexp ) ]]
                                then
					RRead=${OPTARG}
				else
					until [[ "$RRead" =~ $R2regexp ]]
					do
						read -p "Please enter the absolute path to the fastq.gz file for the reverse reads:: " RRead
					done
				fi
				;;
			x)
				if [[ ( -f $OPTARG ) && ( "$OPTARG" =~ $I1regexp ) ]]
				then
					FRead_index=${OPTARG}
				else
					until [[ "$FRead_index" =~ $I1regexp ]]
					do
						read -p "Please enter the absolute path to the fastq.gz file for the indices file for forward reads:: " FRead_index
					done
				fi
				;;
			y)
				if [[ ( -f $OPTARG ) && ( "$OPTARG" =~ $I2regexp ) ]]
				then
					RRead_index=${OPTARG}
				else
					until [[ "$RRead_index" =~ $I2regexp ]]
					do
						read -p "Please enter the absolute path to the fastq.gz file for the indices file for reverse reads:: " RRead_index
					done
				fi
				;;
			s)
				if [[ ( -f $OPTARG ) && ( `basename $OPTARG` =~ $SSregexp ) ]]
				then
					samplesheet=${OPTARG}
				else
					until [[ "$samplesheet" =~ $SSregexp ]]
					do
						read -p "Please enter the absolute path to the sample sheet file containing the indices and sample names of samples loaded into the miseq/nextseq/hiseq:: " samplesheet
					done
				fi
				;;
			?)
				echo "Input error..."
				echo "Usage: ${FUNCNAME[0]} [-a -d] -f <FRead.R1.fastq.gz> -r <RRead.R2.fastq.gz> -x <FRead_index.R1.fastq.gz> -y <RRead_index.R2.fastq.gz> -s <SampleSheet.csv>"
				return 1
				;;
		esac
	done

	until [[ -d ${demldir} ]]
	do
		echo "Creating output directory '${demldir}'"
		mkdir ${demldir}
	done

	#Creating the index file fro the sample sheet
	echo -e "#Index1\tIndex2\tName" > ${indexfile}
	${AWK_EXEC} -v indexfile=${indexfile} 'BEGIN{FS=","; OFS="\t"}; FNR==1,/^Sample_ID/{next}; {print $7,$9,$1 >> indexfile};' ${samplesheet}

	#Demultiplexing
	${deml} -i ${indexfile} -f ${FRead} -r ${RRead} -if1 ${FRead_index} -if2 ${RRead_index} --rgqual 80 --fracconf 20 --wrongness 60 --mm 2 -s ${demldir}/deml.summary -o ${demldir}/ |& tee ${demldir}/deML.log
	if [ $? -eq 0 ]
	then
		echo -e "\tDONE. The output files have been stored in ${demldir}"
	else
		echo -e "\tError occurred while demultiplexing"
	fi
}

#===============================================================================================================================================================
demultiplex_fast5(){ #THis function demultiplexes *.fast5 files based on read_ids of demultiplexed *.fastq files belonging to nanopore data. It uses ont_fast5_api (https://github.com/nanoporetech/ont_fast5_api), "a simple interface to HDF5 files of the Oxford Nanopore .fast5 file format." The ont_fast5_api singularity image is found in */scripts/ont_fast5_api/ont_fast5_api.3.1.6.sif
	
	fqgzregexp='^.*.fastq\.gz$'
	fqregexp='^.*.fastq$'

	INPUTDIR="demultiplex"
	BARCODE_DIRS=("barcode01" "barcode02" "barcode03" "barcode04" "barcode05" "barcode06" "barcode07" "barcode08" "barcode09" "barcode10" "barcode11" "barcode12" "barcode13" "barcode14" "barcode15" "barcode16" "barcode17" "barcode18" "barcode19" "barcode20" "barcode21" "barcode22" "barcode23" "barcode24" "barcode25" "barcode26" "barcode27" "barcode28" "barcode29" "barcode30" "barcode31" "barcode32" "barcode33" "barcode34" "barcode35" "barcode36" "barcode37" "barcode38" "barcode39" "barcode40" "barcode41" "barcode42" "barcode43" "barcode44" "barcode45" "barcode46" "barcode47" "barcode48" "barcode49" "barcode50" "barcode51" "barcode52" "barcode53" "barcode54" "barcode55" "barcode56" "barcode57" "barcode58" "barcode59" "barcode60" "barcode61" "barcode62" "barcode63" "barcode64" "barcode65" "barcode66" "barcode67" "barcode68" "barcode69" "barcode70" "barcode71" "barcode72" "barcode73" "barcode74" "barcode75" "barcode76" "barcode77" "barcode78" "barcode79" "barcode80" "barcode81" "barcode82" "barcode83" "barcode84" "barcode85" "barcode86" "barcode87" "barcode88" "barcode89" "barcode90" "barcode91" "barcode92" "barcode93" "barcode94" "barcode95" "barcode96" "unclassified")
	OUTPUTDIR="demultiplex_fast5/runids"
	for deml_dir in ${BARCODE_DIRS[@]}
	do
		if [ -d ${INDIR}/${INPUTDIR}/${deml_dir} ]
		then
			echo -e "\n\tGenerating a Read ID list for '${INPUTDIR}/${deml_dir}' in '${OUTPUTDIR}/${deml_dir}'"
			mkdir -p ${INDIR}/${OUTPUTDIR}/${deml_dir}
			for fastqfile in ${INPUTDIR}/${deml_dir}/*
			do
				rename ${fastqfile}
				outputfile=${INDIR}/${OUTPUTDIR}/${deml_dir}/${output_filename}.readids
				echo -e "Processing ${input_filename}"
				if [ -f ${outputfile} ]
				then
					rm ${outputfile}
				fi
				if [[ "$fastqfile" =~ $fqgzregexp ]]
				then
					gunzip ${fastqfile}
					${AWK_EXEC} -v outputfile=${outputfile} 'BEGIN{FS=" "}; /^@/{
					if ($2 ~ /^runid/) {gsub(/^@/,"",$1);
						print $1 > outputfile }
					else {next}}' ${fastqfile}
					if [ $? -eq 0 ]
					then
						echo -e "Creating ${output_filename}.readids"
					fi
				elif [[ "$fastqfile" =~ $fqregexp ]]
				then
					${AWK_EXEC} -v outputfile=${outputfile} 'BEGIN{FS=" "}; /^@/{
					if ($2 ~ /^runid/) {gsub(/^@/,"",$1);
						print $1 > outputfile }
					else {next}}' ${fastqfile}
					if [ $? -eq 0 ]
					then
						echo -e "Creating ${output_filename}.readids"
					fi
					else
					continue	
				fi
			done
		fi
	done
}


#===============================================================================================================================================================
concatenate_fasta_seqs() { # This function converts a multiple line FASTA format sequence into a two line record of a header and a sequnce lines.
	#awk '/^>/ {if (FNR==1) print $0; else print "\n" $0; }; !/^>/ {gsub("\n","",$0); printf $0}' renafroCOI_500to700_data-650to660_st22n1006-en1474n3479_head25000_tail125001.afa | less -S


	if [ $# -eq 0 ]
	then
		echo "Input error..."
		echo "Usage: ${FUNCNAME[0]} seqfile1.fasta [seqfile2.fasta seqfile3.fasta ...]"
		return 1
	fi

	for i in "$@"
	do
		if [ ! -f $i ]
		then
			echo "input error: file '$i' is non-existent!"
		elif [[ ( -f $i ) && ( `basename -- "$i"` =~ .*\.(aln|afa|fasta|fa|fst)$ ) ]]
		then
			echo -e "\t...Concatinating sequence lines for each record in `basename -- ${i}`..."
			input_src=`dirname "$( realpath "${i}" )"`
			$AWK_EXEC '/^>/ {if (FNR==1) print $0; else print "\n" $0; }; !/^>/ {gsub("\n","",$0); printf $0}' $i > ${input_src}/outfile.afa #&& mv ${input_src}/outfile.afa ${i}
			if [ $? -eq 0 ]
                        then
				#echo -e "\t moving ${input_src}/outfile.afa to ${input_src}/${i}"
				mv ${input_src}/outfile.afa ${input_src}/${i}
			fi
		else
			echo "input file error in `basename $i`: input file should be a .aln file format"
			continue
		fi
	done
}


#===============================================================================================================================================================
substitute_hdrs() { #This function takes an input file of edited_fasta_format_headers and searches through a fasta_format_sequence file and substitute it's headers if their uniq IDs match

 	if [ $# -eq 0 ]
	then
		echo "Input error..."
		echo -e "\tUsage: ${FUNCNAME[0]} -i 'file1.fasta* [file2.fasta* file3.fasta* ...]' [-h <headers>]"
      		return 1
        fi
	
	unset headers
	unset files

        local OPTIND=1
        while getopts 'i:h:' key
        do
                case "${key}" in
                        i)
                                files=$OPTARG
                                for inputfile in $files
                                do
                                        if [ ! -f $inputfile ]
                                        then
                                                echo -e "\tInput file error in `basename $inputfile`: provided file does not exist"
                                                echo -e "\tUsage: ${FUNCNAME[0]} -i 'file1.fasta* [file2.fasta* file3.fasta* ...]' [-l <INTEGER>] [-m <INTEGERg|G|m|M>]"
                                                return 1
                                        elif [[ ( -f $inputfile ) && ( `basename -- "$inputfile"` =~ .*\.(aln|afa|fasta|fa|fst)$ ) ]]
                                        then
						concatenate_fasta_seqs ${inputfile}
                                                continue
                                        fi
                                done
                                ;;
                        h)
				headers=$OPTARG
				if [ ! -f $headers ]
				then
					echo -e "\tinput file error: the headers file (-h) `basename $inputfile` does not exist"
					until [[ ( -f "$headers" ) && ( `basename -- "$headers"` =~ .*(\.|_)headers$ ) ]]
					do
						echo -e "\nFor the fasta headers input, provide the full path to the file, the filename included."
						read -p "Please enter the file with the FASTA headers:: " headers
					done
				elif [[ ( -f "$headers" ) && ( `basename -- "$headers"` =~ .*(\.|_)headers$ ) ]]
				then
					continue
				else
					echo -e "The file containing headers should be provided via the flag -h <*[\.|_]headers>"
					return 1
				fi
				#$.*/[\]'^
				sed -i "s/\r$//g; s/ /_/g; s/\&/_n_/g; s/\//_/g; s/'//g; s/\[//g; s/\]//g" $headers
				;;
			?)
				echo "Input error..."
				echo 'Usage: ${FUNCNAME[0]} -i "file1.fasta* [file2.fasta* file3.fasta* ...]" [-h <headers>]'
				return 1
				;;
		esac
	done

	
	echo -e "\n\tStarting operation....\n\tPlease wait, this may take a while...."
	for i in "$files"
	do
		unset records
		number_of_replacements=0
		records=$( grep ">" $i | wc -l )
		unset x
		unset y
		unset z
		echo -e "\nProceeding with `basename $i`..."
		for line in `cat ${headers}`
		do
			#x=$( head -10 idrck_headers | tail -1 | awk 'BEGIN { FS="|"; }{print $1;}') && echo $x
			x=`echo "$line" | ${AWK_EXEC} 'BEGIN { RS="\n"; FS="[>|]"; }{ x = $2; print x; }'`
			y=`echo "$line" | ${AWK_EXEC} 'BEGIN { RS="\n"; FS="|"; }{ y = $0; print y; }'`
			#echo -e "\n $x \n $y"

			#Characters to replace from the headers as they will affect the performance of sed: carriage Returns (), white spaces ( ), back slashes (/), and ampersand '&' characters; they greately hamper the next step of header substitution.
			sed -i "s/\r$//g; s/ /_/g; s/\&/_n_/g; s/\//_/g; s/'//g; s/\[//g; s/\]//g" $i

			z=`grep "$x" $i`
			#echo "$z"
			for one_z in `echo -e "${z}"`
			do
				if [ $one_z == $y ]
				then
					echo -e "Change for ${x} already in place..."
					continue
				else
					echo -e "Substituting header for ${x}..."
					sed -i "s/${one_z}/${y}/g" $i
					#sed -i "s/^.*\b${x}\b.*$/${y}/g" $i
				fi
				number_of_replacements=$( expr $number_of_replacements + 1 )
			done
		done
		if [ $? -eq 0 ]
		then
			echo -e "\nDONE. $number_of_replacements replacements done in `basename $i` out of $records records it has"
		fi
	done
	echo -e "\n\tCongratulations...Operation done."
}

#===============================================================================================================================================================
move_matched() { #This is updated version of move_unwanted function 
	# It is used within the filter_seqs function below. It moves and deletes the matching fasta records.
	for i in "$@"
	do
		matching_records=$(grep "$pattern_name" $i | wc -l)
		#echo -e "$i"
		#echo -e "${pattern_name}"
		echo -e "${matching_records} records match the pattern $pattern_name"
		if [ $matching_records -gt 0 ]
		then
			#echo -e "\tSearching headers with '$pattern_name' word-pattern in `basename $i`..."
			$AWK_EXEC -v pat=$pattern_name '/^>/{
			hdr=$0; next} {
			seq=$0 } {
			pattern=pat; gsub(/\|/,"\\|",pattern)} hdr~pattern{
			print hdr; print seq }' $i >> ${input_src}/${output_filename}_${suffix}.${filename_ext}
			
			#Deleting matched records
			regexpno='^n|N|No|NO|no$'
			regexpyes='^y|Y|Yes|YES|yes$'
			if [[ $Choice =~ $regexpyes ]]
			then
				echo -e "Records that match $pattern_name successfully deleted from `basename $i`"
				sed -i "/$pattern_name/,+1 d" $i
			elif [[ ${Choice} =~ $regexno ]]
			then
				echo -e "Records retained in `basename $i`"
			fi
		fi
        done
}


filter_seqs() { # Updated version of delete_unwanted function.
	#this function copys a record that fits a provided pattern, e.g a taxon_name_description; the arguments provided, are the files to be searched for the patterns
	# LOGIC behind the code: To get the list of orders in description_taxon_names and their frequencies, from  which to select the undesired patterns (names), do:
	#grep ">" seqs.fasta | awk 'BEGIN {FS="|"; OFS="|" ; }; {print $2}' |sort | uniq -c > seqs_orders && less seqs_orders
	
	trap ' echo Error $? occured on $LINENO && return 1 ' ERR
	
	if [ $# -eq 0 ]
	then
		echo "Input error..."
		echo "Usage: ${FUNCNAME[0]} file1.*[file2.* file3.* ...]"
		return 1
	fi
        
#	while getopts 'f:r:' key
#	do
#		case "${key}" in
#			f)
	
	for i in "$@"
	do
		 if [ ! -f $i ]
 		 then
 			 echo "input error: file $i is non-existent!"
			 return 1
 		 elif [[ ( -f $i ) && ( `basename -- $i` =~ .*\.(afa|fasta|fa|aln)$ ) ]]
 		 then
	 		 echo -e "\nProceeding with `basename $i`..."
	 		 rename ${i}
	 		 input_src=`dirname "$(realpath '${i}')"`
	 		 concatenate_fasta_seqs ${i}
	 
			 unset options
	 		 echo -e "\nTo extract sequences with specific words in the headers please select one of the options [1|2|3] to proceed or [4] to cancel::\n"
	 		 options[0]="Copy records with word patterns specified in a file into one file"
	 		 options[1]="Copy records with string-patterns specified in a file into individual word-pattern-specific files"
	 		 options[2]="Copy records with specific single string into a file"
	 		 options[3]="Exit"
	 
			 PS3='Select one of the filter methods [1|2|3], or [4] to exit: '
	 		 select option in "${options[@]}"
	 		 do
	 			 unset pattern_name
	 			 unset input_pattern_file
	 			 regexp='^[a-zA-Z0-9/_-\ \|]+$'
	 			 regexp1='^n|y|N|Y|No|Yes|NO|YES|no|yes$'
				 
				 case $option in
	 				 ${options[0]})
	 					 #echo "no error"
	 					 until [[ -f ${input_pattern_file}  ]]
		 				 do
	 						 read -p "Please enter the path to the file with pattern names to be extracted:: " input_pattern_file
	 						 #cat ${input_pattern_file}
	 					 done
						 #deleting matches
	 					 echo -e "\nTo delete records that match the search patterns in '$input_pattern_file' from `basename $i` enter [YES]. Enter [NO] to retain them.\n"
	 					 unset Choice
	 					 read -p "Please enter [Yes] or [NO] to proceed:: " Choice
	 					 until [[ "$Choice" =~ $regexp1 ]]
	 					 do
	 						 echo "INVALID choice $REPLY"
	 					 done
		
						 for line in `cat $(echo ${input_pattern_file})`
	 					 do
							 #echo -e "$line"
							 pattern_name=${line}
	 						 suffix="matched"
	 						 #echo "${pattern_name}"
	 						 move_matched "${i}"
	 					 done
	 					 echo -e "\n\tDONE. All matched records have been stored in '${output_filename}_${suffix}.fasta'"
	 					 break
	 					 ;;
	 				 ${options[1]})
	 					 until [[ -f ${input_pattern_file} ]]
	 					 do
	 						 read -p "Please enter path to the file with pattern names to be extracted:: " input_pattern_file
	 					 done
	 
						 echo -e "\nTo delete records that match the search patterns in '$input_pattern_file' from `basename $i` enter [YES]. Enter [NO] to retain them.\n"
	 					 unset Choice
	 					 read -p "Please enter [Yes] or [NO] to proceed:: " Choice
	 					 until [[ "$Choice" =~ $regexp1 ]]
	 					 do
	 						 echo "INVALID choice $REPLY"
	 					 done
	 
						 for line in `cat ${input_pattern_file}`
	 					 do
	 						 pattern_name=${line}
	 						 suffix=${pattern_name}
		 					 move_matched "${i}"
	 						 echo -e "\n\tDONE. All matched records have been stored in '${output_filename}_${suffix}.fasta'"
	 					 done
	 					 break
	 					 ;;
	 				 ${options[2]})
	 					 until [[ "$pattern_name" =~ $regexp ]]
	 					 do
	 						 read -p "Please enter string pattern to be searched:: " pattern_name
	 					 done
	 
						 echo -e "\nTo delete records that match the search pattern '$pattern_name' from `basename $i` enter [YES]. Enter [NO] to retain them.\n"
	 					 unset Choice
	 					 read -p "Please enter [Yes] or [NO] to proceed:: " Choice
	 					 until [[ "$Choice" =~ $regexp1 ]]
	 					 do
	 						 echo "INVALID choice $REPLY"
	 					 done
			 
	 					 suffix=`echo ${pattern_name} | sed 's/|/\./g'`
	 					 move_matched "${i}"
	 					 echo -e "\n\tDONE. All matched records have been stored in '${output_filename}_${suffix}.fasta'"
	 					 break
	 					 ;;
	 				 ${options[3]})
	 					 echo -e "Exiting filtering of pattern-matched sequences..."
	 					 break
	 					 ;;
	 				 *)
	 					 echo "INVALID choice: Select option [1|2|3] to delete, or [4] to exit: "
	 					 ;;
	 			 esac
	 		 done
		 else
 			 echo "Input file error in `basename $i`: input file should be a .afa file format"
 			 continue
 		 fi
	done
}

#===============================================================================================================================================================

ouso_output() { #This function reorganizes the output as needed by d.ouso's pipeline.
 	if [ $# -eq 0 ]
	then
		echo "Input error..."
		echo -e "\tUsage: ${FUNCNAME[0]} [-i -o ]\n\t-i\tIllumina reads analysis\n\t-o\tONT reads analysis"
      		return 1
	fi
	
	unset gcoverage_path
	unset quast_path
	unset snpeff_path
	unset con_path
	unset mqc_path

	local OPTIND=1
	while getopts 'io' key
	do
		case "${key}" in
			i)
				illumina_regex='fastqc.*$'
				illumina_test=$(grep "fastqc" ${INDIR}/results/pipeline_info/software_versions.*)
				if [[ "$illumina_test" =~ $illumina_regex ]]
				then
					echo -e "\tProceeding with generating aggregation of nf-core/viralrecon results from illumina-reads analysis"
					gcoverage_path=${INDIR}/results/variants/bowtie2/mosdepth
					#echo $gcoverage_path
					quast_path=${INDIR}/results/variants/ivar/quast
					#echo $quast_path
					snpeff_path=${INDIR}/results/variants/ivar/snpeff
					#echo $snpeff_path
					con_path=${INDIR}/results/variants/ivar/consensus
					#echo $con_path
					mqc_path=${INDIR}/results/multiqc
					#echo $mqc_path
				else
					echo -e "\tError!\n\tThe analysis results does not resemble the nf-core/viralrecon analysis output from an illumina sequencing run"
					return 1
				fi
				;;
			o)
				ont_regex='artic.*$'
				ont_test=$(grep "artic" ${INDIR}/results/pipeline_info/software_versions.*)
				if [[ "$ont_test" =~ $ont_regex ]]
				then
					echo -e "\tProceeding with generating aggregation of nf-core/viralrecon results from ONT-reads analysis"
					gcoverage_path=${INDIR}/results/medaka/mosdepth
					quast_path=${INDIR}/results/medaka/quast
					snpeff_path=${INDIR}/results/medaka/snpeff
					mqc_path=${INDIR}/results/multiqc/medaka
					con_path=${INDIR}/results/medaka
				else
					echo -e "\tError!\n\tThe analysis results does not resemble the nf-core/viralrecon analysis output from an ONT sequencing reads"
					return 1
				fi
				;;
			?)
				echo "Input error..."
				echo -e "\tUsage: ${FUNCNAME[0]} [-i -o ]\n\t-i\tIllumina reads analysis\n\t-o\tONT reads analysis"
				return 1
				;;
		esac
	done
	src_dir=${INDIR##*/}
	ousodir=${INDIR}/${src_dir}	
	#echo $ousodir
	if [ -d ${ousodir} ]
	then
		rm -rf ${ousodir}
		echo -e "\t${ousodir} was found and has been deleted"
	fi
	echo -e "\tCreating ${src_dir} and subdirs: nxt png snpEff plt dpt/amplicon dpt/genome qcs var"
	mkdir -p ${ousodir}/nxt ${ousodir}/png ${ousodir}/snpEff ${ousodir}/plt ${ousodir}/dpt/amplicon ${ousodir}/dpt/genome ${ousodir}/qcs ${ousodir}/var
	cp ${gcoverage_path}/amplicon/all_samples.mosdepth.heatmap.pdf ${ousodir}/plt
	echo -e "\tDone copying to */plt"
	cp ${gcoverage_path}/amplicon/*tsv ${ousodir}/dpt/amplicon
	echo -e "\tDone copying to */dpt/amplicon"
	cp ${gcoverage_path}/genome/*tsv ${ousodir}/dpt/genome
	echo -e "\tDone copying to */dpt/genome"
	cp ${quast_path}/transposed_report.tsv ${ousodir}/qcs
	echo -e "\tDone copying to */qcs"
	cp ${INDIR}/pangolin/*.csv ${ousodir}/png
	echo -e "\tDone copying to */png"
	cp ${snpeff_path}/*vcf ${ousodir}/snpEff
	echo -e "\tDone copying to */snpEff"
	cp ${mqc_path}/summary_variants_metrics_mqc.csv ${ousodir}/qcs
	echo -e "\tDone copying to */qcs"
	cp ${con_path}/*.consensus.all.fasta ${ousodir}/png
	echo -e "\tDone copying to */png"
	if [ $? -eq 0 ]
	then
		echo -e "\t${ousodir} has been successfully generated"
	else
		echo -e "\tError copying results. Kindly check"
	fi
	unset INDIR
}
