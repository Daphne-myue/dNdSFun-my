################################################################################################################
# Step 1 : Format the samping mutations from dNdS-Fun into DriverPower and Run DriverPower
# Please find some other related files of DriverPower in the Benchmarked Part
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=Driverpower
#SBATCH --mem=20G
#SBATCH -c 16
#SBATCH -o ./log/Driverpower_%A-%a_out.txt
#SBATCH -e ./log/Driverpower_%A-%a_error.txt
#SBATCH -p amd-ep2,amd-ep2-short
#SBATCH --qos=normal
#SBATCH --array=1-3

## user's own commands below
ID=${1}
seed=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
reg=${2}
numtag=${3}
genenum=${4}
drivernum=${5}

wkDir="/path/PCAWG/"
inputDir="/path/"
mkdir -p ${inputDir}/04_simulateInput_DriverPower/PCAWG/gc19_pc.${reg}/${numtag}_${genenum}
mkdir -p ${inputDir}/05_simulateOutput_DriverPower/PCAWG/gc19_pc.${reg}/${numtag}_${genenum}

mutFil="${inputDir}/04_simulateInput/PCAWG/gc19_pc.${reg}-exon/${numtag}_${genenum}/PCAWG.${reg}-exon.${numtag}.${genenum}.${drivernum}.${seed}.tsv"
DriverPowerFil="${inputDir}/04_simulateInput_DriverPower/PCAWG/gc19_pc.${reg}/${numtag}_${genenum}/PCAWG.${reg}.${numtag}.${genenum}.${drivernum}.${seed}"

echo -e "#chr\tstart\tend\tref\talt\tdonor" > ${DriverPowerFil}".Driverpower.tsv"
cat ${mutFil} | awk '{printf $2"\t""%.f""\t""%.f""\t"$4"\t"$5"\t"$1"\n", $3-1,$3}' >> ${DriverPowerFil}".Driverpower.tsv"
cat ${mutFil} | awk 'BEGIN {OFS="\t"} {print $2, $3-1, $3, $7}' | bedtools intersect -a /path/Driverpower/PCAWG/annotations/gc19_pc.${reg}.excDupSite.bed -b stdin -wa -wb > ${DriverPowerFil}".CADD_ele.tsv"

printf "binID\tCADD\n" > ${DriverPowerFil}".CADD.tsv"
bedtools groupby -i ${DriverPowerFil}".CADD_ele.tsv" -g 4 -c 8 -o mean  >> ${DriverPowerFil}".CADD.tsv"

prepare="/path/DriverPower-1.0.2/prepare.py"
callable="/path/DriverPower-1.0.2/data/callable.bed.gz"
PCAWG_elements="/path/DriverPower-1.0.2/data/PCAWG_test_genomic_elements.bed12.gz"
train_element="/path/Test/train_elements.tsv.gz"
test_element="/path/Driverpower/PCAWG/annotations/gc19_pc.${reg}.excDupSite.bed"
mut_fil=${DriverPowerFil}".Driverpower.tsv"
train_y="/path/Test/train_y.tsv"
test_y=${DriverPowerFil}".Driverpower.test_y.tsv"

# Training responses
python3 ${prepare} ${mut_fil} ${test_element} ${callable} ${test_y}
mkdir -p ${inputDir}/05_simulateOutput_DriverPower/PCAWG/gc19_pc.${reg}/${numtag}_${genenum}/PCAWG.${reg}.${numtag}.${genenum}.${drivernum}.${seed}/
driverpower infer \
    --feature "${wkDir}/test_feature.hdf5" \
    --response ${test_y} \
    --model "/path/Driverpower/PCAWG/output/tutorial.GBM.model.pkl" \
    --name 'DriverPower_burden_function' \
    --outDir "${inputDir}/05_simulateOutput_DriverPower/PCAWG/gc19_pc.${reg}/${numtag}_${genenum}/PCAWG.${reg}.${numtag}.${genenum}.${drivernum}.${seed}/" \
    --funcScore ${DriverPowerFil}".CADD.tsv" \
    --funcScoreCut "CADD:0.01"


################################################################################################################
# Step 2 : Use the hpc to run the scripts multiple times
################################################################################################################
rm run_05_DriverPower.sh
for Ref_num in 1000 10000 100000 1000000 10000000 100000000 5000 50000 500000 5000000
   do for Gene_num in 5 10 20 50 100
      do for Driver_num in 5 10 20 50 100 1000
      do 
        echo "sbatch sinfo_01_DriverPower.sh seed.txt cds ${Ref_num} ${Gene_num} ${Driver_num}" >> run_05_DriverPower.sh
      done
   done
done
