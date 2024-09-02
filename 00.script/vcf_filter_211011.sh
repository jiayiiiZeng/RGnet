#!/usr/bin/bash
dir=$1 #$outputDir/split_bed
input=$2 #vcf file
gene=$3

bcftools query -H -f '%VEP_SYMBOL\t%CHROM\t%POS\t%REF\t%ALT\t%VEP_Feature\t%VEP_HGVSc\t%VEP_HGVSp\t%VEP_IMPACT\t%VEP_Consequence\t%VEP_gnomAD_AFR_AF\t%VEP_gnomAD_AMR_AF\t%VEP_gnomAD_ASJ_AF\t%VEP_gnomAD_EAS_AF\t%VEP_gnomAD_FIN_AF\t%VEP_gnomAD_NFE_AF\t%VEP_gnomAD_SAS_AF[\t%GT]\n' $input > $dir/$gene/$gene'.tsv'
