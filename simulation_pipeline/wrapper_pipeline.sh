#for dir in /home/yutongq/savanna/diploid-gene-annotation/run0413/
#do
#prefix=`basename ${dir} | cut -f 1 -d '_'`

dir=$1
script_dir="/home/yutongq/savanna/SQUID-simulation/true-SV-sim/scripts/"
rm ${dir}/args.txt

echo "getting transcript"
bash ${script_dir}/get_transcript.sh ${dir} 

for prefix in SVRNA100 SVRNA250 SVRNA400
do
    for i in 1 3 5
    do
        echo "${dir} ${prefix} ${i}">> ${dir}/args.txt
    done
done
cat ${dir}/args.txt | xargs -n 3 -P 9 bash ${script_dir}/pipeline.sh
