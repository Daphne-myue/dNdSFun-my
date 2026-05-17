################################################################################################################
# Step 1 : Format the samping mutations from dNdS-Fun into LARVA and Run LARVA
# Please find some other related files of LARVA in the Benchmarked Part
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=OncodriverFML
#SBATCH --mem=20G
#SBATCH -c 8
#SBATCH -o ./log/OncodriverFML_%A-%a_out.txt
#SBATCH -e ./log/OncodriverFML_%A-%a_error.txt
#SBATCH -p amd-ep2,amd-ep2-short,intel-sc3
#SBATCH --qos=normal
#SBATCH --array=1-300

source ~/.bashrc
conda activate OncodriverFML
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

## user's own commands below
ID=${1}
seed=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
reg=${2}
numtag=${3}
genenum=${4}
drivernum=${5}

wkDir="/path/OncodriverFML/PCAWGALL/"
inputDir="/path/"
mkdir -p ${inputDir}/04_simulateInput_OncodriveFML/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}
mkdir -p ${inputDir}/05_simulateOutput_OncodriveFML/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}/PCAWGALL.${reg}.${numtag}.${genenum}.${drivernum}.${seed}

mutFil="${inputDir}/04_simulateInput/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}/PCAWGALL.${reg}.${numtag}.${genenum}.${drivernum}.${seed}.tsv"
OncodriveFMLIn="${inputDir}/04_simulateInput_OncodriveFML/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}/PCAWGALL.${reg}.${numtag}.${genenum}.${drivernum}.${seed}"
OncodriveFMLout="${inputDir}/05_simulateOutput_OncodriveFML/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}/PCAWGALL.${reg}.${numtag}.${genenum}.${drivernum}.${seed}/"

echo -e "CHROMOSOME\tPOSITION\tREF\tALT\tSAMPLE" > ${OncodriveFMLIn}.OncodriverFML.txt
cat ${mutFil} | awk 'sub(/^"chr"/,"");{printf $2"\t""%.f""\t"$4"\t"$5"\t"$1"\n", $3}' | sed 's/^chr//g' >> ${OncodriveFMLIn}.OncodriverFML.txt


region_file="/path/OncodriverFML/PCAWGALL/annotations/${reg}.tsv.gz"
signature_correlation="wg"
oncodrivefml -i ${OncodriveFMLIn}.OncodriverFML.txt --cores 8 -e ${region_file} --signature-correction ${signature_correlation} --seed 123 --force --output ${OncodriveFMLout}

################################################################################################################
# Step 2 : Use the hpc to run the scripts multiple times
################################################################################################################
rm run_05_OncodriverFML.sh
for Ref_num in 1000 10000 100000 1000000 10000000 100000000 5000 50000 500000 5000000
   do for Gene_num in 5 10 20 50 100
      do for Driver_num in 5 10 20 50 100 1000
      do 
        echo "sbatch sinfo_01_OncodriverFML.sh seed.txt cds-exon ${Ref_num} ${Gene_num} ${Driver_num}" >> run_05_OncodriverFML.sh
      done
   done
done

