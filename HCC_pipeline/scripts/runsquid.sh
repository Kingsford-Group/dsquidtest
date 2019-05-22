OutputDir=$(pwd)
scriptDir=$(pwd)/scripts

# SQUID
echo "predicting TSV with SQUID"
mkdir -p ${OutputDir}/TSVprediction/squid_1954
~/SQUID/squid/bin/squid -b ${OutputDir}/StarAlign/HCC1954/Aligned.sortedByCoord.out.bam -c ${OutputDir}/StarAlign/HCC1954/Chimeric.out.bam -o ${OutputDir}/TSVprediction/squid_1954/squidtsv -G 1
awk 'BEGIN{FS="\t";} {if(substr($0,0,1)=="#" || !($1~"M" || $1~"G" || $4~"M" || $4~"G")) print $0;}' ${OutputDir}/TSVprediction/squid_1954/squidtsv_sv.txt > ${OutputDir}/TSVprediction/squid_1954/squidtsv_sv_final.txt
python3.6 ${scriptDir}/VerifyFusionGene.py 1 1954_ref/SV1954_gene.bedpe ${OutputDir}/TSVprediction/squid_1954/squidtsv_sv_final.txt 
python3.6 ${scriptDir}/VerifyFusionGene.py 1 1954_ref/SV1954_nongene.bedpe ${OutputDir}/TSVprediction/squid_1954/squidtsv_sv_final.txt 

mkdir -p ${OutputDir}/TSVprediction/squid_1395
~/SQUID/squid/bin/squid -b ${OutputDir}/StarAlign/HCC1395/Aligned.sortedByCoord.out.bam -c ${OutputDir}/StarAlign/HCC1395/Chimeric.out.bam -o ${OutputDir}/TSVprediction/squid_1395/squidtsv -G 1
awk 'BEGIN{FS="\t";} {if(substr($0,0,1)=="#" || !($1~"M" || $1~"G" || $4~"M" || $4~"G")) print $0;}' ${OutputDir}/TSVprediction/squid_1395/squidtsv_sv.txt > ${OutputDir}/TSVprediction/squid_1395/squidtsv_sv_final.txt
python3.6 ${scriptDir}/VerifyFusionGene.py 1 1395_ref/SV1395_gene.bedpe ${OutputDir}/TSVprediction/squid_1395/squidtsv_sv_final.txt 
python3.6 ${scriptDir}/VerifyFusionGene.py 1 1395_ref/SV1395_nongene.bedpe ${OutputDir}/TSVprediction/squid_1395/squidtsv_sv_final.txt 
# DIPLOID SQUID
echo "predicting TSV with DUPLOID-SQUID"
~/SQUID/diploid-squid/bin/squid -b ${OutputDir}/StarAlign/HCC1954/Aligned.sortedByCoord.out.bam -c ${OutputDir}/StarAlign/HCC1954/Chimeric.out.bam -o ${OutputDir}/TSVprediction/squid_1954/diploidsquidtsv -G 1
awk 'BEGIN{FS="\t";} {if(substr($0,0,1)=="#" || !($1~"M" || $1~"G" || $4~"M" || $4~"G")) print $0;}' ${OutputDir}/TSVprediction/squid_1954/diploidsquidtsv_sv.txt > ${OutputDir}/TSVprediction/squid_1954/diploidsquidtsv_sv_final.txt
python3.6 ${scriptDir}/VerifyFusionGene.py 1 1395_ref/SV1954_gene.bedpe ${OutputDir}/TSVprediction/squid_1954/diploidsquidtsv_sv_final.txt 
python3.6 ${scriptDir}/VerifyFusionGene.py 1 1395_ref/SV1954_nongene.bedpe ${OutputDir}/TSVprediction/squid_1954/diploidsquidtsv_sv_final.txt 
~/SQUID/diploid-squid/bin/squid -b ${OutputDir}/StarAlign/HCC1395/Aligned.sortedByCoord.out.bam -c ${OutputDir}/StarAlign/HCC1395/Chimeric.out.bam -o ${OutputDir}/TSVprediction/squid_1395/diploidsquidtsv -G 1
awk 'BEGIN{FS="\t";} {if(substr($0,0,1)=="#" || !($1~"M" || $1~"G" || $4~"M" || $4~"G")) print $0;}' ${OutputDir}/TSVprediction/squid_1395/diploidsquidtsv_sv.txt > ${OutputDir}/TSVprediction/squid_1395/diploidsquidtsv_sv_final.txt
python3.6 ${scriptDir}/VerifyFusionGene.py 1 1395_ref/SV1395_gene.bedpe ${OutputDir}/TSVprediction/squid_1395/diploidsquidtsv_sv_final.txt 
python3.6 ${scriptDir}/VerifyFusionGene.py 1 1395_ref/SV1395_nongene.bedpe ${OutputDir}/TSVprediction/squid_1395/diploidsquidtsv_sv_final.txt 
