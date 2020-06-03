#!/usr/bin/env bash

plink_name=$1 # full path to plink dataset
output_name=$2 # full path to desired output location

from=$3 # current reference of file (use hg terminology)
to=${4^} # desired reference (with capitalized first letter)

code_directory=$(dirname $0)
chain_directory=${code_directory}/liftover_chain_files #in order to use this function there must be a folder with the desired chain files in the same directory as this script
liftOver=/geschwindlabshares/HumanGenotypeArray/scripts/liftOver/liftOver

chain_file=${chain_directory}/${from}To${to}.over.chain

plink=/share/apps/plink-1.9-beta3c/plink

if [ -f $chain_file ]; then
	echo "Using ${chain_file} as chain file for conversion."
	
	# converting plink file to bed format
	original_variants=${output_name}_${from}.bed
	converted_variants=${output_name}_${4}.bed
	failed_variants=${output_name}_${from}To${to}_failed.bed

	# TODO: check if chromosome code already has "chr"
	awk '{print "chr"$1"\t"$4-1"\t"$4"\t"$2}' ${plink_name}.bim > $original_variants # UCSC Genome Browser bed format

	# liftOver variants
	echo ""
	$liftOver $original_variants $chain_file $converted_variants $failed_variants
	echo ""

	number_converted=$(wc -l < $converted_variants)
	total_variants=$(wc -l < $original_variants)

	if [ $number_converted -gt 0 ]; then
		echo "${number_converted}/${total_variants} successfully converted from ${from} to ${4}."
		echo ""
		
		successful_variants=${output_name}_${4}.snplist
		grep -v _ $converted_variants | awk '{print $4}' > $successful_variants # get all variants with standard chromosome codes
		
		number_successful=$(wc -l < $successful_variants)
		echo "${number_successful}/${number_converted} mapped to standard chromosomes in ${4}."
		echo ""

		# remove "chr" from beginning of chromosome codes
		sed -e "s/chr//g" -i $converted_variants

		#creating new Plink files with converted positions
		$plink --bfile $plink_name --extract $successful_variants --update-chr $converted_variants 1 4 --update-map $converted_variants 3 4 --make-bed --out $output_name

	else
		echo "ERROR: No variants were successfully converted. Please check your chain file and the reference genome of the Plink files."
	fi

else
	echo "ERROR: No chain file found with the name ${chain_file}. Please ensure the file name has not been changed since downloading from UCSC (their naming convention is hgXToHgY.over.chain) and that it has been unzipped."
fi
