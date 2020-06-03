#!/usr/bin/env bash

# Get merged plink file for chosen individuals in given Plink datasets
# Samantha Jensen 2020-08-02
# Updated from Orion script written in 2019-05-08 to be agnostic to server. 

exit_on_error() {
    exit_code=$1
    last_command=$2
    if [ $exit_code -ne 0 ]; then
        >&2 echo "\"${last_command}\" command failed with exit code ${exit_code}."
        exit $exit_code
    fi
}

#############################################################

# define script usage

usage_statement="---------------------------------------------------------------------------------------------------------------------
Usage: $(basename $0) --output_directory PATH --output_name STRING --sample_info FILE --datasets_file FILE --snp_names FILE [--reference_dataset STRING] [--build STRING] [--maf [VALUE]] [--geno [VALUE]] [--hwe [VALUE]] [--mind [VALUE]] [--plink PATH] 
       $(basename $0) --help

 -o, --output_directory		Directory all output will be written to. Will be created if doesn't exist.
 
 -n, --output_name		Project or combined dataset name to use in writing final output. File extensions
				will be added.
 
 -s, --sample_info		Path to PED file with at least the following six tab-separated columns: 
 				FID (family ID), IID (individual ID), Paternal, Maternal, Sex, and Phenotype. 
 				Order matters so any additional fields should be added after these columns.
 
 -d, --datasets_file		Path to file with three tab-separated columns containing the dataset name, the
 				absolute path to that Plink dataset, and the build of the reference.  No header.
 
 -f, --snp_names		Absolute path to a file with \"chromosome position rsid reference alternate\" for 
				standardizing SNP naming conventions. See documentation on Hoffman2 (in 
				/u/home/s/samlj/project-geschwind/ACE/genotypes/orion_code/get_dbsnp_names/) 
				for information on how to generate these files. The positions must be aligned to
				to the desired reference build.

 -r, -- reference_dataset	If also merging a reference dataset like HapMap or 1000G, specify the name of 
				this dataset as designated in the datasets_file. (Optional)

 -b, --build			The desired reference genome build of your output merged Plink dataset. If
				not specified, all input Plink files must have the same build. (Optional)

 -m, --maf			Minor allele frequency threshold. Will be passed to Plink to filter any variant
				with allele frequency less than this value. (Optional)

 -g, --geno			Variant missingness threshold. Will be passed to Plink to filter any variants with
				missing genotypes in more than this proportion of samples. (Optional)

 -p, --hwe			Hardy Weinberg Equilibrium p-value threshold. Will be passed to Plink to exclude
				variants with HWE exact test p-values below this point. (Optional)

 -i, --mind			Individual missingness threshold. Will be passed to Plink to filter any individuals
				with missing genotypes for more than this proportion of variants. (Optional)

 --plink			Set path to plink executable. If not set, script will assume plink is in the path.
				(Optional)

-h, --help			Display this help message.

NOTE: Order of parameters does not matter.
---------------------------------------------------------------------------------------------------------------------"

#############################################################

# setting constants

plink=$(type -P plink) # this variable will be used to call Plink and should be modified with the --plink flag if the executable is not in your $PATH
if [[ $plink == "" ]]; then
	echo ""
	echo "WARNING: The Plink executable is not in your path. You must set the location of Plink with the --plink flag."
	echo ""
fi

code_directory=$(dirname $0) # supporting scripts must be in the same directory as this script
find_complementary_flips="${code_directory}/find_complementary_flips.R"
convert_reference="${code_directory}/convert_reference.sh"

#############################################################

# using getopt to parse arguments

options=`getopt --options o:n:s:d:f:r:b:m:g:p:i:h --long output_directory:,output_name:,sample_info:,datasets_file:,snp_names:,reference_dataset:,build:,maf:,geno:,hwe:,mind:,plink:,help -n 'merge_and_filter.sh' -- "$@"`

eval set -- "$options"

for argument in $options ; do
	case "$1" in
		-o|--output_directory) output_directory=$2; shift 2 ;;
		-n|--output_name) output_name=$2; shift 2 ;;	
		-s|--sample_info) sample_info=$2; shift 2 ;;
		-d|--datasets_file) datasets_file=$2; shift 2 ;;
		-f|--snp_names) snp_names=$2; shift 2 ;;
		-r|--reference_dataset) reference_dataset=$2; shift 2 ;;
		-b|--build) build=$2; shift 2 ;;
		-m|--maf) maf=$2; shift 2 ;;
		-g|--geno) geno=$2; shift 2 ;;
		-p|--hwe) hwe=$2; shift 2 ;;
		-i|--mind) mind=$2; shift 2 ;;
		--plink) plink=$2; shift 2 ;;
		-h|--help) echo "$usage_statement"; exit 0 ;;
		--) shift ; break ;;
		*) echo "ERROR: Argument parsing failed! Check getopt usage to see if it has changed." ; exit 1 ;;
	esac
done

#############################################################

# print log header

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "RUNNING MERGE_AND_FILTER"
echo "$(date)"
echo "FROM $(pwd)"
echo "WITH $plink"
echo "" # add space in output
echo "$0 $@" # script call with parameters
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""

#############################################################

echo "#############################################################################################"
echo "STEP 1: Check input and set variables"
echo "#############################################################################################"
echo ""

# check that all required variables were set
if [[ -z $output_directory || -z $output_name || -z $sample_info || -z $datasets_file || -z $snp_names ]]; then
        echo "ERROR: Required variable not set!"
	echo ""
	echo "$usage_statement"
	exit 1
fi

# TODO: more sophisticated checks of parameters

plink_output_directory=${output_directory}/plink_files
mkdir -p $plink_output_directory # create output directories if they don't already exist

dataset_names=($(awk '{print $1}' $datasets_file))
datasets=($(awk '{print $2}' $datasets_file))

echo "Project ${output_name}:"
echo "Combining ${#dataset_names[@]} Plink datasets specified in ${datasets_file}"
echo "Filtering to $(wc -l < $sample_info) samples specified in ${sample_info}"
echo "Writing output in ${output_directory}"

if [ ! -z $reference_dataset ]; then
	echo "Including ${reference_dataset} samples."
fi

echo ""

if [ ! -z $build ]; then
	echo "Converting all datasets to reference build ${build}. Using column three of ${datasets_file} to determine the current build of each input dataset."
	dataset_builds=($(awk '{print $3}' $datasets_file)) # third column should have each file's reference build
else
	echo "Assuming all input datasets have the same reference build."	
fi

echo ""

#############################################################

echo "#############################################################################################"
echo "STEP 2: Create file with family and individual sample IDs for filtering"
echo "#############################################################################################"
echo ""

id_file=${output_directory}/${output_name}_IDs.txt

if [ ! -f $id_file ]; then
	echo "Creating file with QC-passed IDs at ${id_file}."
	awk '{print $1,$2}' $sample_info > $id_file # get first two columns
	sed -i '1D'  $id_file # remove header line
else
	echo "ID file already exists. Skipping this step."
fi
echo ""

#############################################################

echo "#############################################################################################"
echo "STEP 3: Get file that defines all individuals as founders for use with --flip-scan"
echo "#############################################################################################"
echo ""

fake_founders=${output_directory}/${output_name}_no_parents.txt

if [ ! -f $fake_founders ]; then
	echo "Creating fake founders file, ${fake_founders}, for --flip-scan."
	cp $id_file $fake_founders
	# add two columns of 0s
	sed -i 's/$/ 0 0/g' $fake_founders	
else
	echo "Founder file for --flip-scan already exists. Skipping this step."
fi
echo ""

############################################################

echo "#############################################################################################"
echo "STEP 4: Format each dataset before merging."
echo "#############################################################################################"
echo ""

cd ${plink_output_directory} # change to plink output directory so merged Plink files will be in their own subdirectory

for index in ${!datasets[@]}; #getting indices so we can get info from both arrays
do
	dataset_file=${datasets[$index]}
	dataset_name=${dataset_names[$index]}


	echo "----- formatting ${dataset_name} from Plink files (${dataset_file}.*) -----"
	echo ""

	
	echo "---------------------------------------------------------------------------------------------"
	echo "A) Check that all files have the same reference build."
	echo "---------------------------------------------------------------------------------------------"
	echo ""

	if [ -f ${dataset_name}.bed ]; then
		echo "${dataset_name} Plink files with the correct reference build have already been created in ${plink_output_directory}. Moving on to next formatting step."
	else
		if [ ! -z $build ]; then
			dataset_build=${dataset_builds[$index]}
	
			if [[ $dataset_build == $build ]]; then
				echo "$dataset_name already has reference build $build according to ${datasets_file}."
	
				cp ${dataset_file}.bed ${dataset_name}.bed
				cp ${dataset_file}.bim ${dataset_name}.bim
				cp ${dataset_file}.fam ${dataset_name}.fam
			else
				echo "Converting $dataset_name from $dataset_build to reference build ${build}."
				echo ""

				$convert_reference $dataset_file ${dataset_name} $dataset_build $build
				exit_on_error $? convert_reference
			fi
		else
			echo "No --build set. Assuming all input datasets have the same reference build."	
	
			cp ${dataset_file}.bed ${dataset_name}.bed
			cp ${dataset_file}.bim ${dataset_name}.bim
			cp ${dataset_file}.fam ${dataset_name}.fam
		fi
	fi
	
	echo ""
	
	echo "---------------------------------------------------------------------------------------------"
	echo "B) Select only samples of interest."
	echo "---------------------------------------------------------------------------------------------"
	echo ""

	filtered_name=${dataset_name}_filtered

	if [ -f ${filtered_name}.bed ]; then
		echo "${dataset_name} files have already been created in ${plink_output_directory}. Moving to next filter step."

	else 
		if [ "$reference_dataset" = "$dataset_name" ]; then
			echo "${dataset_name} is a reference dataset. Copying all samples to ${plink_output_directory}/${filtered_name}."

			cp ${dataset_name}.bed ${filtered_name}.bed
			cp ${dataset_name}.bim ${filtered_name}.bim
			cp ${dataset_name}.fam ${filtered_name}.fam
		else
 
			$plink --bfile $dataset_name --keep $id_file --allow-no-sex --make-bed --out $filtered_name
		fi
	fi
	
	echo ""

	echo "---------------------------------------------------------------------------------------------"
	echo "C) Keep only mapped variants."
	echo "---------------------------------------------------------------------------------------------"
	echo ""
	
	mapped_name=${dataset_name}_mapped

	if [ ! -f ${mapped_name}.bed ]; then
		$plink --bfile $filtered_name --not-chr 0 --make-bed --out $mapped_name # 0 is unknown chromosome code
	else
		echo "$dataset_name already filtered to only mapped variants. See ${mapped_name}. Moving to next filter step."
	fi

	echo ""

	echo "---------------------------------------------------------------------------------------------"
	echo "C) Remove duplicate variants."
	echo "---------------------------------------------------------------------------------------------"
	echo ""

	unduplicated_name=${dataset_name}_unduplicated
	
	if [ ! -f ${unduplicated_name}.bed ]; then

		duplicates=${mapped_name}_duplicates

		echo "Outputting list of duplicated variants to ${duplicates}.dupvar."
		echo ""
	
		$plink --bfile $mapped_name --list-duplicate-vars ids-only --out $duplicates
		echo ""

		awk '{print $1}' ${duplicates}.dupvar > ${duplicates}.out # first column has kgp Illumina style IDs
		number_duplicates=$(wc -l < ${duplicates}.out)
						
		if [ $number_duplicates -gt 0 ]; then	
			echo "Excluding $number_duplicates Illumina style IDs (see ${duplicates}.out for IDs)."
			echo ""

			$plink --bfile $mapped_name --exclude ${duplicates}.out --allow-no-sex --make-bed --out $unduplicated_name
		else
			echo "No synonymous duplicates in ${mapped_name}."
		
			cp ${mapped_name}.bed ${unduplicated_name}.bed
			cp ${mapped_name}.bim ${unduplicated_name}.bim
			cp ${mapped_name}.fam ${unduplicated_name}.fam
		fi

		echo ""

	else
		echo "Variants included more than once in Plink files already removed. See ${unduplicated_name}. Moving to next filter step."
	fi

	echo ""

	echo "---------------------------------------------------------------------------------------------"
	echo "D) Remove all other sites with more than one mapped variant."
	echo "---------------------------------------------------------------------------------------------"
	echo ""

	unique_name=${dataset_name}_unique_positions

	if [ ! -f ${unique_name}.bed ]; then
		multivariant_sites=${unduplicated_name}_multivariant_sites

		echo "Checking for other variants sharing the same position."
		echo ""
		
		# bim file format: 1=CHR, 2=SNP_NAME, 3=CM_POSITION, 4=BP_POSITION, 5=ALLELE_1, 6=ALLELE_2
		awk 'n=x[$1,$3,$4]{print n"\t"$2"\t"$5"\t"$6;} {x[$1,$3,$4]=$0;}' ${unduplicated_name}.bim > ${multivariant_sites}.txt # output duplicated sites by position on one row per pair
		awk '{print $2"\n"$7}' ${multivariant_sites}.txt > ${multivariant_sites}.out # output all IDs of duplicated sites on their own rows
		number_multivariant_sites=$(wc -l < ${multivariant_sites}.out)

		if [ $number_multivariant_sites -gt 0 ]; then	
			echo "Removing $number_multivariant_sites other variants sharing the same position. See ${multivariant_sites}.txt for full details and ${multivariant_sites}.out for the list of excluded SNPs."
			echo ""

			$plink --bfile $unduplicated_name --exclude ${multivariant_sites}.out --allow-no-sex --make-bed --out $unique_name
		else
			echo "No other multivariant sites in ${unduplicated_name}."
			
			cp ${unduplicated_name}.bed ${unique_name}.bed
			cp ${unduplicated_name}.bim ${unique_name}.bim
			cp ${unduplicated_name}.fam ${unique_name}.fam
		fi
	else
		echo "${unique_name} Plink files already exist in ${plink_output_directory}. Skipping filtering step for this dataset."
	fi
	echo ""

	echo "---------------------------------------------------------------------------------------------"
	echo "E) Standardize variant names and alleles."
	echo "---------------------------------------------------------------------------------------------"
	echo ""

	batch_names=${dataset_name}_${build}_names_and_alleles.txt	

	if [ ! -f ${batch_names} ]; then
		echo "Reference file ${batch_names} has not yet been created."
		echo "Generating a reference for all $(wc -l ${unique_name}.bim) SNPs genotyped in ${dataset_name} using reference file ${snp_names}."
		echo ""

		snp_positions=${dataset_name}_snp_positions.txt
		awk '{print $1"\t"$4"\t"}' ${unique_name}.bim > ${snp_positions}
		#awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' ${unique_name}.bim > ${dataset_name}_snps.txt
		matching_snps="${dataset_name}_${build}_matching_snps.txt"
		grep -F -f $snp_positions $snp_names > $matching_snps
		#TODO: process matching snps and get correct names and alleles
	else
		echo "Reference file ${batch_names} already exists. We will use this to rename all variants in ${dataset_name}."
	fi

	echo ""

	final_name=${dataset_name}_final

#	if [ ! -f ${final_name}.bed ]; then
#		echo "Generating new variant names. Each variant will be assigned CHR:POSITION as its variant name to simplify merging."
#		echo ""
#		awk '{print $2"\t"$1":"$4}' ${unique_name}.bim > ${unique_name}_new_names.txt
#
#		echo "Renaming all variants."
#		echo ""
#
#		$plink --bfile $unique_name --update-name ${unique_name}_new_names.txt --allow-no-sex --make-bed --out $final_name
#	else
#		echo " ${final_name} Plink files already exist in ${plink_output_directory}. Skipping variant renaming step for this dataset."
#	fi
#	echo ""
done
#
##############################################################
#
#echo "#############################################################################################"
#echo "STEP 5: Combine all datasets"
#echo "#############################################################################################"
#
#echo ""
#
#number_datasets=${#dataset_names[@]}
#
#if [[ $number_datasets > 0 ]]; then
#	combined=${dataset_names[0]}
#	for (( index=1; index<=$number_datasets-1; index++)); 
#	do
#		dataset1=${combined}_final
#		dataset2=${dataset_names[$index]}_final
#		combined="${combined}_${dataset_names[$index]}"
#		
#		echo "----- combining" $dataset1 "and" $dataset2 "-----"
#		echo "" 
#
#		# check if has been run before
#		if [ ! -f ${combined}_final.bed ]; then
#			echo "---------------------------------------------------------------------------------------------"
#			echo "A) Try to merge"
#			echo "---------------------------------------------------------------------------------------------"
#			echo ""
#
#			$plink --bfile $dataset1 --bmerge $dataset2 --merge-mode 2 --allow-no-sex --make-bed --out ${combined}_initial_merge # NOTE: can't use --merge-equal-pos option because it chooses whichever name is lexicographically first (in this case it will be the kgp ID from Illumina)
#
#			echo ""
#
#			if [ ! -f ${combined}_initial_merge-merge.missnp ]; then # no errors
#				echo "No errors: combined Plink files copied to ${combined}_final."
#				echo ""
#
#				cp ${combined}_initial_merge.bed ${combined}_final.bed
#				cp ${combined}_initial_merge.bim ${combined}_final.bim
#				cp ${combined}_initial_merge.fam ${combined}_final.fam	
#			else
#				echo "Strandedness errors likely."
#				echo ""
#				echo "---------------------------------------------------------------------------------------------"
#				echo "B) Check for strandedness problems by flipping the strand of variants with tri-allelic errors and trying to merge again"
#				echo "---------------------------------------------------------------------------------------------"
#				echo ""
#
#				$plink --bfile $dataset1 --flip ${combined}_initial_merge-merge.missnp --allow-no-sex --make-bed --out ${dataset1}_initial_flip
#				echo ""
#
#				$plink --bfile ${dataset1}_initial_flip --bmerge $dataset2 --merge-mode 2 --allow-no-sex --make-bed --out ${combined}_initial_flip_merge
#				echo ""
#				
#				if [ ! -f ${combined}_initial_flip_merge-merge.missnp ]; then # no tri-allelic variants
#					triallelic=0
#
#					cp ${dataset1}_initial_flip.bed ${dataset1}_merge_ready.bed
#					cp ${dataset1}_initial_flip.bim ${dataset1}_merge_ready.bim
#					cp ${dataset1}_initial_flip.fam ${dataset1}_merge_ready.fam
#					
#					cp ${dataset2}.bed ${dataset2}_merge_ready.bed
#					cp ${dataset2}.bim ${dataset2}_merge_ready.bim
#					cp ${dataset2}.fam ${dataset2}_merge_ready.fam
#						
#					cp ${combined}_initial_flip_merge.bed ${combined}_merged.bed
#					cp ${combined}_initial_flip_merge.bim ${combined}_merged.bim
#					cp ${combined}_initial_flip_merge.fam ${combined}_merged.fam
#				
#				else
#					triallelic=$(wc -l ${combined}_initial_flip_merge-merge.missnp | awk '{print $1}')
#
#					$plink --bfile ${dataset1}_initial_flip --exclude ${combined}_initial_flip_merge-merge.missnp --allow-no-sex --make-bed --out ${dataset1}_merge_ready
#					echo ""
#					
#					$plink --bfile $dataset2 --exclude ${combined}_initial_flip_merge-merge.missnp --allow-no-sex --make-bed --out ${dataset2}_merge_ready
#					echo ""
#
#					$plink --bfile ${dataset1}_merge_ready --bmerge ${dataset2}_merge_ready --merge-mode 2 --allow-no-sex --make-bed --out ${combined}_merged
#					echo ""
#				fi
#
#				echo "True tri-allelic variants: $triallelic"
#				echo ""
#				
#				echo "---------------------------------------------------------------------------------------------"
#				echo "C) Check remaining A/T C/G strand inconsistencies"
#				echo "---------------------------------------------------------------------------------------------"
#				echo ""
#
#				merged_length=$(wc -l ${combined}_merged.fam | awk '{print $1}')
#				dataset1_length=$(wc -l ${dataset1}_merge_ready.fam | awk '{print $1}')
#
#				echo "Merged length: $merged_length"
#				echo "Dataset 1 length: $dataset1_length"
#				echo ""
#
#				if (( $merged_length > $dataset1_length )); then
#
#					echo "Using LD to identify strand inconsistencies in A/T and C/G SNPs (with --flip-scan) that aren't caught during the merge."
#					echo ""
#
#					$plink --bfile ${combined}_merged --update-parents $fake_founders --make-pheno ${dataset1}_merge_ready.fam '*' --flip-scan --out ${combined}_potential_remaining_problems
#					echo ""
#
#					sed -i 's/NA *$/NA NA/g' ${combined}_potential_remaining_problems.flipscan # last column only exists if SNPs are in LD	
#
#					Rscript $find_complementary_flips ${combined}_potential_remaining_problems.flipscan 0.5 1 ${combined}_strandedness_problems.snplist
#					remaining_problems=$(wc -l ${combined}_strandedness_problems.snplist | awk '{print $1}')
#				else
#					echo "WARNING: This merge didn't add any new people to the dataset, so we can't look for strand inconsistencies with --flip-scan."
#					echo ""
#					remaining_problems=0
#				fi
#
#				echo "A/T C/G strand inconsistencies: $remaining_problems"
#				echo ""
#
#				if [[ $remaining_problems != 0 ]]; then
#					echo "Flipping remaining problems."
#					echo ""				
#					
#					$plink --bfile ${dataset1}_merge_ready --flip ${combined}_strandedness_problems.snplist --allow-no-sex --make-bed --out ${dataset1}_merge_ready_final_flip
#
#					echo ""
#					echo "Remerging with flipped SNPs."
#					echo ""
#
#					# remerge
#					$plink --bfile ${dataset1}_merge_ready_final_flip --bmerge ${dataset2}_merge_ready --merge-mode 2 --allow-no-sex --make-bed --out ${combined}_final
#
#				else
#					echo "No additional detected strand inconsistencies. Copying ${combined}_merged to ${combined}_final."
#					echo ""
#
#					cp ${combined}_merged.bed ${combined}_final.bed
#					cp ${combined}_merged.bim ${combined}_final.bim
#					cp ${combined}_merged.fam ${combined}_final.fam	
#				fi
#			fi
#			
#		else
#			echo "Plink files for the merge of ${dataset1} and ${dataset2} already exist. Skipping this run."
#		fi
#
#		echo ""
#	done	
#fi
#
#############################################################
#
#echo "#############################################################################################"
#echo "STEP 6: Get new variant names"
#echo "#############################################################################################"
#
#echo ""
#
#if [ ! -f ${snp_names} ];
#	echo "Reference file ${snp_names} has not yet been created. Now generating this file from dbsnp144."
#	echo "This step is time consuming and memory intensive."
#	echo ""
#
#	awk '{print $2}' ${combined}_final.bim > ${combined}_final_snps.txt
#	
#else
#	echo "Reference file ${snp_names} already exists. We will use this to rename all variants."
#fi
#
#echo ""
#
##############################################################
#
#echo "#############################################################################################"
#echo "STEP 6: Rename all variants"
#echo "#############################################################################################"
#
#echo ""
#
#renamed=${combined}_renamed
#
#if [ ! -f ${renamed}.bed ]; then
#	echo "Renaming all variants based on reference ${snp_names}."
#	echo ""
#
#	$plink --bfile ${combined}_final --update-name $snp_names --allow-no-sex --make-bed --out $renamed
#
#else
#	echo "Variant renaming step has already occurred. Skipping this step."
#fi
#
#echo ""
#
#############################################################
#
#echo "#############################################################################################"
#echo "STEP 6: Create phenotype file from sample info"
#echo "#############################################################################################"
#
#echo ""
#
#phenotypes=${output_directory}/${output_name}.pheno
#
#if [ ! -f $phenotypes ]; then
#	echo "Creating phenotype file (${phenotypes}) from sample info file (${sample_info})."
#	echo ""
#
#	awk '{print $1,$2,$6}' $sample_info > $phenotypes # get family ID, IID, and phenotype
#	sed -i '1D'  $phenotypes # remove header line
#else
#	echo "Phenotype file has already been written. Skipping this step."
#fi
#
#echo ""
#
##############################################################
#
#echo "#############################################################################################"
#echo "STEP 7: Add phenotypes to combined Plink files"
#echo "#############################################################################################"
#
#echo ""
#
#pheno_name=${combined}_pheno
#
#if [ ! -f ${pheno_name}.bed ]; then
#	$plink --bfile ${renamed} --make-pheno $phenotypes 2 --allow-no-sex --make-bed --out ${pheno_name} 
#else
#	echo "Phenotypes have already been added. Skipping this step."
#fi
#echo ""
#
##############################################################
#
#echo "#############################################################################################"
#echo "STEP 8: Quality control"
#echo "#############################################################################################"
#
#echo ""
#
#echo "---------------------------------------------------------------------------------------------"
#echo "A) Filter by minor allele frequency"
#echo "---------------------------------------------------------------------------------------------"
#echo ""
#
#
#if [ ! -z $maf ]; then
#
#	echo "Filtering genotypes with a minor allele frequency less than ${maf}."
#	echo ""
#	
#	maf_name=${pheno_name}_maf_$maf
#
#	if [ ! -f ${maf_name}.bed ]; then
#		$plink --bfile ${pheno_name} --maf $maf --allow-no-sex --make-bed --out ${maf_name}
#	else
#		echo "Genotypes with MAF < $maf have already been filtered: see ${maf_name} Plink files."
#	fi
#else
#	maf_name=${pheno_name}
#	echo "--maf not set. Not filtering combined file based on minor allele frequency."
#fi
#
#echo ""
#
#echo "---------------------------------------------------------------------------------------------"
#echo "B) Filter by genotype missingness"
#echo "---------------------------------------------------------------------------------------------"
#echo ""
#
#
#if [ ! -z $geno ]; then
#
#	echo "Filtering genotypes with missing values in more than ${geno} of samples."
#	echo ""
#	
#	geno_name=${maf_name}_geno_$geno
#
#	if [ ! -f ${geno_name}.bed ]; then
#		$plink --bfile ${maf_name} --geno $geno --allow-no-sex --make-bed --out ${geno_name}
#	else
#		echo "Genotypes with more than $geno missing samples have already been filtered: see ${geno_name} Plink files."
#	fi
#else
#	geno_name=${maf_name}
#	echo "--geno not set. Not filtering genotypes in combined file based on missing samples."
#fi
#
#echo ""
#
#echo "---------------------------------------------------------------------------------------------"
#echo "C) Filter by Hardy-Weinburg equilibrium"
#echo "---------------------------------------------------------------------------------------------"
#echo ""
#
#
#if [ ! -z $hwe ]; then
#
#	echo "Filtering genotypes not in HWE (p<$hwe)."
#	echo ""
#	
#	hwe_name=${geno_name}_hwe_$hwe
#
#	if [ ! -f ${hwe_name}.bed ]; then
#		$plink --bfile ${geno_name} --hwe $hwe --allow-no-sex --make-bed --out ${hwe_name}
#	else
#		echo "Genotypes not in HWE have already been filtered: see ${hwe_name} Plink files."
#	fi
#else
#	hwe_name=${geno_name}
#	echo "--hwe not set. Not filtering genotypes in combined file based on HWE."
#
#fi
#
#echo ""
#
#echo "---------------------------------------------------------------------------------------------"
#echo "D) Filter by sample missingness"
#echo "---------------------------------------------------------------------------------------------"
#echo ""
#
#
#if [ ! -z $mind ]; then
#
#	echo "Filtering samples with more than ${mind} missing genotypes."
#	echo ""
#	
#	mind_name=${hwe_name}_mind_$mind
#
#	if [ ! -f ${mind_name}.bed ]; then
#		$plink --bfile ${hwe_name} --mind $mind --allow-no-sex --make-bed --out ${mind_name}
#	else
#		echo "Samples with more than $mind missing genotypes have already been filtered: see ${mind_name} Plink files."
#	fi
#else
#	mind_name=${hwe_name}
#	echo "--mind not set. Not filtering samples in combined file based on missing genotypes."
#
#fi
#
#echo ""
#
##############################################################
#
#echo "#############################################################################################"
#echo "STEP 9: Create VCF"
#echo "#############################################################################################"
#
#echo ""
#
#final_output=${output_directory}/${output_name}
#
#if [ ! -f ${final_output}.vcf ]; then
#	$plink --bfile $mind_name --recode vcf --out ${final_output}
#else
#	echo "VCF with name ${final_output}.vcf has already been created. Skipping this step."
#fi
#
#echo ""
#
##############################################################
#
#echo "#############################################################################################"
#echo "STEP 10: Copy final merged and filtered Plink files to desired output name"
#echo "#############################################################################################"
#
#echo ""
#
#if [ ! -f ${final_output}.bed ]; then
#	cp ${mind_name}.bed ${final_output}.bed
#	cp ${mind_name}.bim ${final_output}.bim
#	cp ${mind_name}.fam ${final_output}.fam
#else
#	echo "Plink files for ${output_name} already exist in ${output_directory}. Skipping this step."
#fi
#
#echo "Final merged and filtered output is in:
# 1. VCF file - ${final_output}.vcf
# 2. Plink files -  ${final_output}.bed, ${final_output}.bim, and ${final_output}.fam"
#echo ""
#
##############################################################
#
## print end of log
#
#echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
#echo "MERGE_AND_FILTER finished"
#echo "$(date)"
#echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
#echo ""
#
