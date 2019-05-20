#for dir in /home/yutongq/savanna/diploid-gene-annotation/run0413/
#do
#prefix=`basename ${dir} | cut -f 1 -d '_'`


dir=$1
prefix=$2
i=$3

script_dir="/home/yutongq/savanna/SQUID-simulation/true-SV-sim/scripts/"

echo ${dir}
echo ${prefix}
echo rep${i}
# RSVSIM

# annotation

# read simulation
echo "Simulating Reads..."
mkdir -p ${dir}/reads/${prefix}_rep${i}
Rscript-3.3.2 ${script_dir}/SimulationReads.R ${dir}/${prefix}_${i}_transcripts.fa ${dir}/reads/${prefix}_rep${i}
perl ${script_dir}fasta_to_fastq.pl ${dir}/reads/${prefix}_rep${i}/sample_01_1.fasta > ${dir}/reads/${prefix}_rep${i}_1.fastq
perl ${script_dir}fasta_to_fastq.pl ${dir}/reads/${prefix}_rep${i}/sample_01_2.fasta > ${dir}/reads/${prefix}_rep${i}_2.fastq
#perl ${script_dir}fasta_to_fastq.pl ${dir}/sample_02_1.fasta > ${dir}/sample_1_02_1.fastq
#perl${script_dir}fasta_to_fastq.pl ${dir}/sample_02_2.fasta > ${dir}/sample_1_02_2.fastq

mkdir -p ${dir}/reads/${prefix}_rep$((i+1))
Rscript-3.3.2 ${script_dir}/SimulationReads.R ${dir}/${prefix}_$((i+1))_transcripts.fa ${dir}/reads/${prefix}_rep$((i+1))
perl ${script_dir}fasta_to_fastq.pl ${dir}/reads/${prefix}_rep$((i+1))/sample_01_1.fasta > ${dir}/reads/${prefix}_rep$((i+1))_1.fastq
perl ${script_dir}fasta_to_fastq.pl ${dir}/reads/${prefix}_rep$((i+1))/sample_01_2.fasta > ${dir}/reads/${prefix}_rep$((i+1))_2.fastq
#perl ${script_dir}fasta_to_fastq.pl ${dir}/sample_02_1.fasta > ${dir}/sample_2_02_1.fastq
#perl ${script_dir}fasta_to_fastq.pl ${dir}/sample_02_2.fasta > ${dir}/sample_2_02_2.fastq

bash ${script_dir}/filter_fastq.sh ${dir}/reads/${prefix}_rep${i}_1.fastq ${dir}/reads/${prefix}_rep${i}_1.fastq 50 ${dir}/reads ${prefix}_rep${i}
bash ${script_dir}/filter_fastq.sh ${dir}/reads/${prefix}_rep${i}_2.fastq ${dir}/reads/${prefix}_rep${i}_2.fastq 50 ${dir}/reads ${prefix}_rep${i}
bash ${script_dir}/filter_fastq.sh ${dir}/reads/${prefix}_rep$((i+1))_1.fastq ${dir}/reads/${prefix}_rep$((i+1))_1.fastq 50 ${dir}/reads ${prefix}_rep$((i+1))
bash ${script_dir}/filter_fastq.sh ${dir}/reads/${prefix}_rep$((i+1))_2.fastq ${dir}/reads/${prefix}_rep$((i+1))_2.fastq 50 ${dir}/reads ${prefix}_rep$((i+1))

# align
mkdir -p ${dir}/aligned
echo "Aligning..."
STAR --runThreadN 16 --genomeDir /home/yutongq/savanna/SQUID-simulation/diploid-gene-annotation/STAR_genome --readFilesIn ${dir}/reads/${prefix}_rep${i}_1.50.fastq ${dir}/reads/${prefix}_rep${i}_2.50.fastq --outFileNamePrefix ${dir}/aligned/${prefix}_${i}. --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --chimSegmentMin 20 --outSAMstrandField intronMotif --limitBAMsortRAM 29787429000
STAR --runThreadN 16 --genomeDir /home/yutongq/savanna/SQUID-simulation/diploid-gene-annotation/STAR_genome --readFilesIn ${dir}/reads/${prefix}_rep$((i+1))_1.50.fastq ${dir}/reads/${prefix}_rep$((i+1))_2.50.fastq --outFileNamePrefix ${dir}/aligned/${prefix}_$((i+1)). --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --chimSegmentMin 20 --outSAMstrandField intronMotif --limitBAMsortRAM 29787429052
samtools view -bS ${dir}/aligned/${prefix}_${i}.Chimeric.out.sam > ${dir}/aligned/${prefix}_${i}.Chimeric.out.bam
samtools view -bS ${dir}/aligned/${prefix}_$((i+1)).Chimeric.out.sam > ${dir}/aligned/${prefix}_$((i+1)).Chimeric.out.bam

# run squid on haploid
mkdir -p ${dir}/hap_predictions
echo "------ Running SQUID on haploid for ${prefix}_rep${i} ----"
/home/yutongq/SQUID/squid/bin/squid -b ${dir}/aligned/${prefix}_${i}.Aligned.sortedByCoord.out.bam  \
                                    -c ${dir}/aligned/${prefix}_${i}.Chimeric.out.bam \
                                    -o ${dir}/hap_predictions/squid_${prefix}_rep${i} -G 1 > ${dir}/hap_predictions/squid_${prefix}_rep${i}.log.txt

echo "------ Running SQUID on haploid for ${prefix}_rep$((i+1)) ----"
/home/yutongq/SQUID/diploid-squid/bin/squid -b ${dir}/aligned/${prefix}_$((i+1)).Aligned.sortedByCoord.out.bam \
                                    -c ${dir}/aligned/${prefix}_$((i+1)).Chimeric.out.bam \
                                    -o ${dir}/hap_predictions/squid_${prefix}_rep$((i+1)) -G 1 > ${dir}/hap_predictions/squid_${prefix}_rep$((i+1)).log.txt

echo "---- Running DIPLOID SQUID for ${prefix}_rep${i}----"
/home/yutongq/SQUID/diploid-squid/bin/squid -b ${dir}/aligned/${prefix}_${i}.Aligned.sortedByCoord.out.bam  \
                                    -c ${dir}/aligned/${prefix}_${i}.Chimeric.out.bam \
                                    -o ${dir}/hap_predictions/diploid_${prefix}_rep${i}  -G 1 > ${dir}/hap_predictions/diploid_${prefix}_rep${i}.log.txt

echo "---- Running DIPLOID SQUID for ${prefix}_rep$((i+1))----"
/home/yutongq/SQUID/diploid-squid/bin/squid -b ${dir}/aligned/${prefix}_$((i+1)).Aligned.sortedByCoord.out.bam  \
                                    -c ${dir}/aligned/${prefix}_$((i+1)).Chimeric.out.bam \
                                    -o ${dir}/hap_predictions/diploid_${prefix}_rep$((i+1)) -G 1> ${dir}/hap_predictions/diploid_${prefix}_rep$((i+1)).log.txt

# evaluate haploid squid
python3.6 ~/SQUID/squidtest/src/VerifySVpred.py 1 ${dir}/${prefix}_${i}_adjacency.tsv ${dir}/hap_predictions/squid_${prefix}_rep${i}_sv.txt ${dir}/hap_predictions/squid_${prefix}_rep${i}_out.txt ${dir}/hap_predictions/squid_${prefix}_rep${i}_true.txt
python3.6 ~/SQUID/squidtest/src/VerifySVpred.py 1 ${dir}/${prefix}_$((i+1))_adjacency.tsv ${dir}/hap_predictions/squid_${prefix}_rep$((i+1))_sv.txt ${dir}/hap_predictions/squid_${prefix}_rep$((i+1))_out.txt ${dir}/hap_predictions/squid_${prefix}_rep$((i+1))_true.txt
python3.6 ~/SQUID/squidtest/src/VerifySVpred.py 1 ${dir}/${prefix}_${i}_adjacency.tsv ${dir}/hap_predictions/diploid_${prefix}_rep${i}_sv.txt ${dir}/hap_predictions/diploid_${prefix}_rep${i}_out.txt ${dir}/hap_predictions/diploid_${prefix}_rep${i}_true.txt
python3.6 ~/SQUID/squidtest/src/VerifySVpred.py 1 ${dir}/${prefix}_$((i+1))_adjacency.tsv ${dir}/hap_predictions/diploid_${prefix}_rep$((i+1))_sv.txt ${dir}/hap_predictions/diploid_${prefix}_rep$((i+1))_out.txt ${dir}/hap_predictions/diploid_${prefix}_rep$((i+1))_true.txt

# merge bam files
mkdir -p ${dir}/merged_bam
samtools merge -f ${dir}/merged_bam/${prefix}_rep${i}$((i+1)).Aligned.merged.bam ${dir}/aligned/${prefix}_$((i+1)).Aligned.sortedByCoord.out.bam ${dir}/aligned/${prefix}_${i}.Aligned.sortedByCoord.out.bam
samtools merge -f ${dir}/merged_bam/${prefix}_rep${i}$((i+1)).Chimeric.merged.bam ${dir}/aligned/${prefix}_$((i+1)).Chimeric.out.bam ${dir}/aligned/${prefix}_${i}.Chimeric.out.bam
samtools view -s 0.5 -bS ${dir}/merged_bam/${prefix}_rep${i}$((i+1)).Aligned.merged.bam > ${dir}/merged_bam/${prefix}_rep${i}$((i+1)).Aligned.merged.half.bam
samtools view -s 0.5 -bS ${dir}/merged_bam/${prefix}_rep${i}$((i+1)).Chimeric.merged.bam > ${dir}/merged_bam/${prefix}_rep${i}$((i+1)).Chimeric.merged.half.bam

# run squid
mkdir -p ${dir}/dip_predictions
echo "---- Running SQUID for ${prefix}_rep${i}$((i+1)) ----"
/home/yutongq/SQUID/squid/bin/squid -b ${dir}/merged_bam/${prefix}_rep${i}$((i+1)).Aligned.merged.half.bam \
                                    -c ${dir}/merged_bam/${prefix}_rep${i}$((i+1)).Chimeric.merged.half.bam \
                                    -o ${dir}/dip_predictions/squid_${prefix}_rep${i}$((i+1)) -G 1 > ${dir}/dip_predictions/squid_${prefix}_rep${i}$((i+1)).log.txt


echo "---- Running DIPLOID SQUID for ${prefix}_rep${i}$((i+1))----"
/home/yutongq/SQUID/diploid-squid/bin/squid -b ${dir}/merged_bam/${prefix}_rep${i}$((i+1)).Aligned.merged.half.bam \
                                    -c ${dir}/merged_bam/${prefix}_rep${i}$((i+1)).Chimeric.merged.half.bam \
                                -o ${dir}/dip_predictions/diploid_${prefix}_rep${i}$((i+1)) -G 1 > ${dir}/dip_predictions/diploid_${prefix}_rep${i}$((i+1)).log.txt

# merge truth files
cat ${dir}/${prefix}_${i}_adjacency.tsv ${dir}/${prefix}_$((i+1))_adjacency.tsv >  ${dir}/${prefix}_${i}$((i+1))_adjacency.tsv

# evaluate squid
python3.6 ~/SQUID/squidtest/src/VerifySVpred.py 1 ${dir}/${prefix}_${i}$((i+1))_adjacency.tsv ${dir}/dip_predictions/diploid_${prefix}_rep${i}$((i+1))_sv.txt ${dir}/dip_predictions/diploid_${prefix}_rep${i}$((i+1))_out.txt ${dir}/dip_predictions/diploid_${prefix}_rep${i}$((i+1))_true.txt
python3.6 ~/SQUID/squidtest/src/VerifySVpred.py 1 ${dir}/${prefix}_${i}$((i+1))_adjacency.tsv ${dir}/dip_predictions/squid_${prefix}_rep${i}$((i+1))_sv.txt ${dir}/dip_predictions/squid_${prefix}_rep${i}$((i+1))_out.txt ${dir}/dip_predictions/squid_${prefix}_rep${i}$((i+1))_true.txt
