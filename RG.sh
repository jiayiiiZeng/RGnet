#!/bin/bash

###The following paths and input files need to be modified according to the actual situation###
input=/public/home/zengjiayi/work/RGnet/01.data/GJB2.vep_control_split_csq.vcf.gz ##The input vcf file that has been annotated by VEP.
outputDir=/public/home/zengjiayi/work/RGnet/03.output_control ##The path to the output files.
inputdir=/public/home/zengjiayi/work/RGnet/02.input ##The path to the input files.
codeDir=/public/home/zengjiayi/work/RGnet/00.script ##The path to the scripts.
result=$outputDir/result ##The path to results.
###
inputsample=$inputdir/samplelist.txt ##The name of samples in the cohort.
samplesize=2504 ##Number of samples in the cohort.
inputgene=$inputdir/GJB2.txt ##The name of GJB2.
inputbed=$inputdir/GJB2.bed ##The bed file of GJB2.
inputTranscript=$inputdir/GJB2_trans.txt ##The given transcript of GJB2.
input_p=$inputdir/P_LP_240813.txt ##The information of P sites and LP sites in GJB2.
input_b=$inputdir/B_LB_240813.txt ##The information of B sites and LB sites in GJB2.
log=$outputDir/01.log/GJB2.log ##The log of the process.
###

mkdir -p $outputDir/01.log
touch $log
echo "start at $(date)" > $log 2>&1
mkdir -p $outputDir/random $outputDir/random/source/ $outputDir/random/homo_split_gene/ $outputDir/random/graph/ $outputDir/random/pempirical/
outputDir1=$outputDir/random
outputDir2=$outputDir/random/source/
inputdir2=$outputDir/random/homo_split_gene/
graphdir=$outputDir/random/graph/
mkdir -p $outputDir/random/graph/
mkdir -p $outputDir/split_bed/  $outputDir/split_bed/temp
cat $inputgene|while read line;do rm -rf $outputDir/split_bed/$line;mkdir -p $outputDir/split_bed/$line;done >>$log 2>&1
echo "created files finished at $(date)" > $log 2>&1

cat $inputgene|xargs -P 1 -n 1 python $codeDir/split.bed.py $inputbed $outputDir/split_bed/ >>$log 2>&1
echo " split bed finished at $(date)" >>$log 2>&1

cat $inputgene|xargs -P 1 -n 1 bash $codeDir/vcf_filter_211011.sh $outputDir/split_bed $input >>$log 2>&1
echo "extract vcf finished at $(date)" >>$log 2>&1

cat $inputgene|xargs -P 1 -n 1 python $codeDir/filter.240814.py $outputDir/split_bed/  $inputTranscript >>$log 2>&1
echo "filter snp finished at $(date)" >>$log 2>&1

cat $inputgene|xargs -P 1 -n 1 python $codeDir/trans.py $outputDir/split_bed/ >>$log 2>&1
echo "trans finished at $(date)" >>$log 2>&1

cat $inputgene|xargs -P 1 -n 1 python $codeDir/homo_deg_1110_case005.py  $outputDir1/pempirical/  $outputDir/split_bed/  10000  $samplesize  >>$log 2>&1
echo "random finished at $(date)" >>$log 2>&1

cat $inputgene|xargs -P 1 -n 1 python $codeDir/split.2.py $outputDir1/pempirical/ $inputdir2 >>$log 2>&1
echo "split_2 finished at $(date)" >>$log 2>&1

cat $inputgene|xargs -P 1 -n 1 python $codeDir/source_graph.210926.2.py 10000 $outputDir1/pempirical/  $graphdir/  $outputDir2  $outputDir/split_bed/  $input_p $input_b >>$log 2>&1
echo "graph finished at $(date)" >>$log 2>&1
