#!/bin/bash

argv=("$@")

OutputDir="" # downloaded reference and RNA-seq folder, the same as output folder
threads=4
ExeDir=$(pwd)

GenomeFasta=${OutputDir}/Ensemble75/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa
AnnotationGTF=${OutputDir}/Ensemble75/Homo_sapiens.GRCh37.75.gtf

# preparing STAR index
echo "Generating STAR index"
mkdir -p ${OutputDir}/Ensemble75/STARIndex
STAR --runThreadN ${threads} --runMode genomeGenerate --genomeDir ${OutputDir}/Ensemble75/STARIndex --genomeFastaFiles ${GenomeFasta} --sjdbGTFfile ${AnnotationGTF}

cd $OutputDir
# align reads using STAR
echo "STAR alignment"
mkdir -p ${OutputDir}/StarAlign/HCC1954
STAR  --runThreadN ${threads} --genomeDir ${OutputDir}/Ensemble75/STARIndex --readFilesIn ${OutputDir}/RNAseq/RNAHCC1954_1.fastq ${OutputDir}/RNAseq/RNAHCC1954_2.fastq --outFileNamePrefix $OutputDir/StarAlign/HCC1954/ --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --chimSegmentMin 15 --outReadsUnmapped Fastx
samtools view -Shb ${OutputDir}/StarAlign/HCC1954/Chimeric.out.sam -o ${OutputDir}/StarAlign/HCC1954/Chimeric.out.bam

mkdir -p ${OutputDir}/StarAlign/HCC1395
STAR  --runThreadN ${threads} --genomeDir ${OutputDir}/Ensemble75/STARIndex --readFilesIn ${OutputDir}/RNAseq/SRR2532336_1.fastq ${OutputDir}/RNAseq/SRR2532336_2.fastq --outFileNamePrefix $OutputDir/StarAlign/HCC1395/ --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --chimSegmentMin 15 --outReadsUnmapped Fastx
samtools view -Shb ${OutputDir}/StarAlign/HCC1395/Chimeric.out.sam -o ${OutputDir}/StarAlign/HCC1395/Chimeric.out.bam
