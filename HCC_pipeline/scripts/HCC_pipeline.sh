OutputDir=$1
scriptDir=$(pwd)

bash ${scriptDir}/downloadHCC.sh ${OutputDir}
bash ${scriptDir}/prepare_HCC.sh ${OutputDir}
bash ${scriptDir}/runsquid.sh ${OutputDir}
