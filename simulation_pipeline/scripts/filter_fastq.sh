#!/bin/bash

#credit: http://seqanswers.com/forums/showthread.php?t=31845
# This removes reads of a below a certain length from paired read files in fastq format (e.g., R1 and R2 from the same library)
# Usage: $ bash nixshorts_PE [input fastqR1] [input fastqR2] [minimum read length to keep] 

#1. Start with inputs
R1fq=$1
R2fq=$2
minlen=$3
outdir=$4
prefix=$5

prefix1=`basename $1 | cut -f 1 -d '.'`
prefix2=`basename $2 | cut -f 1 -d '.'`

#2. Find all entries with read length less than minimum length and print line numbers, for both R1 and R2
awk -v min=$minlen '{if(NR%4==2) if(length($0)<min) print NR"\n"NR-1"\n"NR+1"\n"NR+2}' $R1fq > temp${prefix}.lines1
awk -v min=$minlen '{if(NR%4==2) if(length($0)<min) print NR"\n"NR-1"\n"NR+1"\n"NR+2}' $R2fq >> temp${prefix}.lines1

#3. Combine both line files into one, sort them numerically, and collapse redundant entries
sort -n temp${prefix}.lines1 | uniq > temp${prefix}.lines
rm temp${prefix}.lines1

#4. Remove the line numbers recorded in "lines" from both fastqs
awk 'NR==FNR{l[$0];next;} !(FNR in l)' temp${prefix}.lines $R1fq > ${outdir}/$prefix1.$minlen.fastq
awk 'NR==FNR{l[$0];next;} !(FNR in l)' temp${prefix}.lines $R2fq > ${outdir}/$prefix2.$minlen.fastq
rm temp${prefix}.lines
