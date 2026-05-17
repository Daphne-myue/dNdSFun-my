################################################################################################################
# Step 1 : Format the samping mutations from dNdS-Fun into LARVA and Run LARVA
# Please find some other related files of LARVA in the Benchmarked Part
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=LARVA
#SBATCH --mem=40G
#SBATCH -c 16
#SBATCH -o ./log/LARVA_%A-%a_out.txt
#SBATCH -e ./log/LARVA_%A-%a_error.txt
#SBATCH -p amd-ep2,amd-ep2-short
#SBATCH --qos=normal
#SBATCH --array=1-3

module load R/4.0.5
## user's own commands below
ID=${1}
seed=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
reg=${2}
numtag=${3}
genenum=${4}
drivernum=${5}

for drivernum in 5 10 20 50 100 1000
do

wkDir="/path/LARVA/PCAWGALL"
inputDir="/path/"
mkdir -p  ${inputDir}/05_simulateOutput_LARVA/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}/PCAWGALL.${reg}.${numtag}.${genenum}.${drivernum}.${seed}/
cd ${inputDir}/05_simulateOutput_LARVA/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}/PCAWGALL.${reg}.${numtag}.${genenum}.${drivernum}.${seed}/

mutFil="${inputDir}/04_simulateInput/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}/PCAWGALL.${reg}.${numtag}.${genenum}.${drivernum}.${seed}.tsv"
LARVAIn="${inputDir}/04_simulateInput_LARVA/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}/PCAWGALL.${reg}.${numtag}.${genenum}.${drivernum}.${seed}"
LARVAout="${inputDir}/05_simulateOutput_LARVA/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}/PCAWGALL.${reg}.${numtag}.${genenum}.${drivernum}.${seed}"

cat ${mutFil} | awk '{printf $2"\t""%.f""\t""%.f""\t"$1"\t"$4"\n", $3-1,$3}' > ${LARVAIn}.txt
/path/larva -vf ${LARVAIn}.txt -af ${wkDir}/annotations/gc19_pc.${reg}.bed -o ${LARVAout}.out

cd /path/01_script_others
rm -r ${inputDir}/05_simulateOutput_LARVA/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}/PCAWGALL.${reg}.${numtag}.${genenum}.${drivernum}.${seed}/

done

################################################################################################################
# Step 2 : Use the hpc to run the scripts multiple times
################################################################################################################
rm run_05_LARVA.sh
for Ref_num in 1000 10000 100000 1000000 10000000 100000000 5000 50000 500000 5000000
   do for Gene_num in 5 10 20 50 100
      do for Driver_num in 5 10 20 50 100 1000
      do 
        echo "sbatch sinfo_01_LARVA.sh seed.txt cds-exon ${Ref_num} ${Gene_num} ${Driver_num}" >> run_05_LARVA.sh
      done
   done
done
