dir=$1
script_dir=/home/yutongq/savanna/SQUID-simulation/true-SV-sim/scripts/
for f in ./SVRNA*
do
    prefix=`basename $f`
    echo ${prefix}
    cd ${f}
    echo "Post processing..."
    python3.6 ${script_dir}"/post_process.py" 
    cd ..
    echo "Get Transcriptome..."
    python3.6 ${script_dir}"/RearrangedTranscripts_3.py" ${prefix}/final genome_original.fa genes_clean.gtf ${dir}/${prefix}
done

