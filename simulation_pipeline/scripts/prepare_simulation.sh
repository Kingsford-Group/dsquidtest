#!/bin/bash

args=("$@")
if [ ${#args[@]} -ne 2 ]; then
    echo "!!!Usage:"
    echo "         bash script/runWholeSimulation.sh <outdir> <threads>"
    exit 1
fi
OutDir=${args[0]}
threads=${args[1]}
ExecutableDir=$(pwd)

mkdir -p OutDir
mkdir -p $OutDir/GRCh38
mkdir -p $OutDir/GRCh38/Chromosomes
mkdir -p $OutDir/GRCh38/Annotation
mkdir -p $OutDir/WholeGenome

# download genome and reference
cd $OutDir/GRCh38/Chromosomes
for ((i=1; i<=5; i++)); do
    if [ ! -f ./${i}.fa ]; then
        echo `ls ./Homo_sapiens.GRCh38.dna.chromosome.${i}.fa.gz`
        wget ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${i}.fa.gz
        gunzip Homo_sapiens.GRCh38.dna.chromosome.${i}.fa.gz
        mv Homo_sapiens.GRCh38.dna.chromosome.${i}.fa ${i}.fa
        cat ${i}.fa >> ${OutDir}/genome_original.fa
    fi
done
samtools faidx ${OutDir}/genome_original.fa | cut -f 1,2 > chrom_sizes.txt
cd ../Annotation
if [ ! -f Homo_sapiens.GRCh38.89.gtf ]; then
    wget ftp://ftp.ensembl.org/pub/release-89/gtf/homo_sapiens/Homo_sapiens.GRCh38.89.gtf.gz
    gunzip Homo_sapiens.GRCh38.89.gtf.gz
fi
awk '{if(length($1)==1 && $1>=1 && $1<=5) print $0}' Homo_sapiens.GRCh38.89.gtf > genes.gtf
awk 'BEGIN{FS="\t";OFS="\t"}{split($NF,a," ");pfx="";s="";for(i=1;i<=length(a);i+=2){if(a[i]=="transcript_id"){pfx=a[i]" "a[i+1]}else{s=s" "a[i]" "a[i+1]}}if(pfx==""){}else{$NF=pfx""s;print$0} }' genes.gtf > ${OutDir}/genes_clean.gtf
