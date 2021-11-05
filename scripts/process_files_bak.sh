#!/usr/bin/env bash

AWK_EXEC=$( which gawk )
#set -u

#INDIR=$(pwd)
rename() { #Takes input file name and extracts the prefix of the name ("prefix"), The suffix (".suffix"), Absolute path to the file and the name of the Directory within which the file is found.
	input_filename=`basename -- ${i}`
	output_filename=${input_filename%.*}
	filename_ext=${input_filename##*.}
	src_dir_path=`dirname $(realpath ${i})`
	src_dir=${src_dir_path##*/}
	#echo -e "input is $i \ninput_filename is $input_filename \noutput_filename is $output_filename \nfilename_ext is $filename_ext \nsrc_dir_path is $src_dir_path \nsrc_dir is $src_dir"
}

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
				elif [[ ( "$OPTARG" =~ .*(\.fastq\.gz)$ ) && ( "$(ls -A $INDIR/*$OPTARG)" ) ]]
				then
					f_ext=$OPTARG
					fastqgz_suf=$f_ext
					#Generating an array containing all forward read .fastq.gz files in the working dir. based on the provided suffix
					array=($INDIR/*$fastqgz_suf)
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
				elif [[ ( $OPTARG =~ .*(\.fastq\.gz)$ ) && ( "$(ls "$INDIR"/*"$OPTARG")" ) ]]
                                then
					r_ext=$OPTARG
                                        fastqgz_suf=$r_ext
					#Generating an array containing all reverse read *fastq.gz files from the working dir. based on the provided suffix
					array=($INDIR/*$fastqgz_suf)
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
	echo "sample,fastq_1,fastq_2" > "$INDIR/samplesheet.csv"
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

				echo "${f_prefix}",$file,$file2 >> "$INDIR/samplesheet.csv"
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

variant_calls() { #This file compiles a list of variant annotations in a .tsv file.
	#Error trapping
	trap ' echo Error $? occured on $LINENO && return 1' ERR

	variants_dir=$INDIR/${results}/variants/ivar/snpeff
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
						gsub("p.","",aamute)}};
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
