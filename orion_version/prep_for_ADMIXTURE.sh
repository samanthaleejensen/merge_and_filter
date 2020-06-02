#!/usr/bin/env bash

# Get QC-filtered plink and pop files for use with PCA and ADMIXTURE
# Samantha Jensen 2019-05-08 updated 2019-09-13 and 2019-11-12
# (original is /home/sjensen/ACE_841/merge_and_filter.sh)

usage_statement="---------------------------------------------------------------------------------------------------------------------
Usage: $(basename $0) --plink_file PATH --output_directory PATH --output_name STRING 
       $(basename $0) --help

 -p, --plink_file		Merged and filtered Plink file with reference populations and sample genotypes.

 -o, --output_directory		Directory all output will be written to. Will be created if doesn't exist.
 
 -n, --output_name		Project or combined dataset name to use in writing final output. File extensions
				will be added.

 -h, --help			Display this help message.

NOTE: Order of parameters does not matter.
---------------------------------------------------------------------------------------------------------------------"

#############################################################

# using getopt to parse arguments

options=`getopt --options p:o:n:h --long plink_file:,output_directory:,output_name:,help -n 'prep_for_ADMIXTURE.sh' -- "$@"`

eval set -- "$options"

for argument in $options ; do
	case "$1" in
		-p|--plink_file) plink_file=$2; shift 2 ;;
		-o|--output_directory) output_directory=$2; shift 2 ;;
		-n|--output_name) output_name=$2; shift 2 ;;	
		-h|--help) echo "$usage_statement"; exit 0 ;;
		--) shift ; break ;;
		*) echo "ERROR: Argument parsing failed! Check getopt usage to see if it has changed." ; exit 1 ;;
	esac
done

#############################################################

# setting constants

plink=/share/apps/plink-1.9-beta3c/plink

#############################################################

# print log header

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "RUNNING PREP_FOR_ADMIXTURE"
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
if [[ -z $plink_file || -z $output_directory || -z $output_name ]]; then
        echo "ERROR: Required variable not set!"
	echo ""
	echo "$usage_statement"
	exit 1
fi

plink_output_directory=${output_directory}/plink_files
mkdir -p $plink_output_directory # create output directories if they don't already exist

echo "Project ${output_name}:"
echo "Using merged and filtered Plink file ${plink_file}."
echo "Writing output in ${output_directory} and intermediate Plink files in ${plink_output_directory}."

echo ""

#############################################################

echo "#############################################################################################"
echo "STEP 2: PCA"
echo "#############################################################################################"
echo ""

pca_name=${plink_file}.eigenvec

if [ ! -f $pca_name ]; then
	$plink --bfile $plink_file --pca --out $plink_file
else
	"PCA has already been carried out for this Plink dataset: see ${plink_file}.eigenvec and ${plink_file}.eigenval."
fi

echo ""

#############################################################

echo "#############################################################################################"
echo "STEP 3: Filter genotypes for ADMIXTURE analysis"
echo "#############################################################################################"

echo ""

echo "ADMIXTURE (http://software.genetics.ucla.edu/admixture/admixture-manual.pdf)
works best with SNPs genotyped in all samples and not in LD. It also requires
unrelated samples."

cd $plink_output_directory # moving into plink_files so all intermediate output will be written here

echo "---------------------------------------------------------------------------------------------"
echo "A) Get only autosomes"
echo "---------------------------------------------------------------------------------------------"
echo ""

autosome_name=${plink_file}_autosomes

if [ ! -f ${autosomes_name}.bed ]; then
	$plink --bfile ${plink_file} --chr 1-22 --allow-no-sex --make-bed --out ${autosome_name}
else
	echo "Autosomes have already been filtered: ${autosome_name}."
	echo ""
fi

echo "---------------------------------------------------------------------------------------------"
echo "B) Get only SNPs in linkage equilibrium"
echo "---------------------------------------------------------------------------------------------"
echo ""

LD_threshold=0.1 # filter variants with R^2 > 0.1 with each other
LD_name=${autosome_name}_LD_${LD_threshold}

if [ ! -f ${LD_Name}.bed ]; then
	$plink --bfile ${autosome_name} --indep-pairwise 50 100 $LD_threshold --allow-no-sex --make-bed --out ${LD_name}
else
	echo "LD pruning already completed for ${plink_file}."
	echo ""
fi

echo "---------------------------------------------------------------------------------------------"
echo "C) Extract only founder individuals"
echo "---------------------------------------------------------------------------------------------"
echo ""

founders=${LD_name}_founders

if [ ! -f ${founders}.bed ]; then
	$plink --bfile ${LD_name} --filter-founders --allow-no-sex --make-bed --out ${founders}
else
	echo "Founders have already been extracted for admixture analysis in ${plink_file}."
	echo ""
fi

echo "---------------------------------------------------------------------------------------------"
echo "D) Check relatedness of all founder samples"
echo "---------------------------------------------------------------------------------------------"
echo ""

relatedness_cutoff=0.25 # grandparent/grandchild/aunt/uncle/niece/nephew/half-sibling
unrelated_founders=${founders}_unrelated_$relatedness_cutoff

if [ ! -f ${unrelated_founders}.bed ]; then
	$plink --bfile ${founders} --rel-cutoff $relatedness_cutoff --allow-no-sex --make-bed --out $unrelated_founders # remove individuals with greater than 2nd degree relatives in dataset
else
	echo "Unrelated individuals have already been removed for ${plink_file} founders."
	echo ""
fi

echo "---------------------------------------------------------------------------------------------"
echo "E) Extract only cases"
echo "---------------------------------------------------------------------------------------------"
echo ""

cases=${LD_name}_cases

if [ ! -f ${cases}.bed ]; then
	$plink --bfile ${LD_name} --filter-cases --allow-no-sex --make-bed --out ${cases}
else
	echo "Cases have already been extracted for admixture analysis in ${plink_file}."
	echo ""
fi

echo "---------------------------------------------------------------------------------------------"
echo "F) Check relatedness of all cases"
echo "---------------------------------------------------------------------------------------------"
echo ""

unrelated_cases=${cases}_unrelated_$relatedness_cutoff

if [ ! -f ${unrelated_cases}.bed ]; then
	$plink --bfile ${cases} --rel-cutoff $relatedness_cutoff --allow-no-sex --make-bed --out $unrelated_cases # remove individuals with greater than 2nd degree relatives in dataset
else
	echo "Unrelated individuals have already been removed for ${plink_file} cases."
	echo ""
fi

echo "---------------------------------------------------------------------------------------------"
echo "G) Check relatedness of all samples"
echo "---------------------------------------------------------------------------------------------"
echo ""

unrelated_all=${LD_name}_all_unrelated_$relatedness_cutoff

if [ ! -f ${unrelated_all}.bed ]; then
	$plink --bfile ${LD_name} --rel-cutoff $relatedness_cutoff --allow-no-sex --make-bed --out $unrelated_all # remove individuals with greater than 2nd degree relatives in dataset
else
	echo "All unrelated individuals have already been removed for ${plink_file}."
	echo ""
fi

#############################################################

echo "#############################################################################################"
echo "STEP 4: Copy final Plink files to desired output name"
echo "#############################################################################################"

echo ""

# copy file that's got only unrelated founders
cp ${unrelated_founders}.bed ${output_directory}/${output_name}_founders_unrelated.bed
cp ${unrelated_founders}.bim ${output_directory}/${output_name}_founders_unrelated.bim
cp ${unrelated_founders}.fam ${output_directory}/${output_name}_founders_unrelated.fam

# copy file that's got only unrelated cases
cp ${unrelated_cases}.bed ${output_directory}/${output_name}_cases_unrelated.bed
cp ${unrelated_cases}.bim ${output_directory}/${output_name}_cases_unrelated.bim
cp ${unrelated_cases}.fam ${output_directory}/${output_name}_cases_unrelated.fam

# copy file that has everyone
cp ${LD_name}.bed ${output_directory}/${output_name}_all.bed
cp ${LD_name}.bim ${output_directory}/${output_name}_all.bim
cp ${LD_name}.fam ${output_directory}/${output_name}_all.fam

# copy file that has all unrelated samples
cp ${unrelated_all}.bed ${output_directory}/${output_name}_all_unrelated.bed
cp ${unrelated_all}.bim ${output_directory}/${output_name}_all_unrelated.bim
cp ${unrelated_all}.fam ${output_directory}/${output_name}_all_unrelated.fam

echo "Final output for ADMIXTURE analysis is in:
 1. Unrelated founders -  ${output_directory}/${output_name}_founders_unrelated.bed, ${output_directory}/${output_name}_founders_unrelated.bim, and ${output_directory}/${output_name}_founders_unrelated.fam
 2. Unrelated cases - ${output_directory}/${output_name}_cases_unrelated.bed, ${output_directory}/${output_name}_cases_unrelated.bim, ${output_directory}/${output_name}_cases_unrelated.fam 
 3. All ACE samples - ${output_directory}/${output_name}_all.bed, ${output_directory}/${output_name}_all.bim, ${output_directory}/${output_name}_all.fam
 4. Unrelated ACE samples -  ${output_directory}/${output_name}_all_unrelated.bed, ${output_directory}/${output_name}_all_unrelated.bim, ${output_directory}/${output_name}_all_unrelated.fam"
echo ""

#############################################################

echo "#############################################################################################"
echo "STEP 5: Create pop files for supervised learning"
echo "#############################################################################################"

echo ""

echo "Creating pop files for supervised learning with ADMIXTURE."

awk '{print $1}' ${output_directory}/${output_name}_founders_unrelated.fam > ${output_directory}/${output_name}_founders_unrelated.pop
sed -i 's/AU.*$/-/g' ${output_directory}/${output_name}_founders_unrelated.pop

awk '{print $1}' ${output_directory}/${output_name}_cases_unrelated.fam > ${output_directory}/${output_name}_cases_unrelated.pop
sed -i 's/AU.*$/-/g' ${output_directory}/${output_name}_cases_unrelated.pop

awk '{print $1}' ${output_directory}/${output_name}_all.fam > ${output_directory}/${output_name}_all.pop
sed -i 's/AU.*$/-/g' ${output_directory}/${output_name}_all.pop

awk '{print $1}' ${output_directory}/${output_name}_all_unrelated.fam > ${output_directory}/${output_name}_all_unrelated.pop
sed -i 's/AU.*$/-/g' ${output_directory}/${output_name}_all_unrelated.pop

#############################################################

# print end of log

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "PREP_FOR_ADMIXTURE finished"
echo "$(date)"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""
