dir=$1
out_name=$2
echo -e "name\tobj\tfrac"> $out_name
for f in ${dir}/*log*txt
do
    obj=`grep OBJ $f |  awk -F "\t" '{SUM+=$2}END{print SUM}'`
    all=`tail -3 ${f} | head -1 | awk -F ": " '{print $2}'`
    res=`tail -3 ${f} | tail -1 | awk -F ": " '{print $2}'`
    frac=`bc -l <<< ${res}/${all}`
    prefix=`basename $f | cut -f 1 -d "."`
    echo -e "${prefix}\t${obj}\t${frac}" >> $out_name
done
