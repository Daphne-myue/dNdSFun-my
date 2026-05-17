################################################################################################################
# step 1 : format the genomic annotations used in dNdS-Fun into the DriverPower
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=Driverpower
#SBATCH --ntasks-per-node=2
###SBATCH --mem=20G
#SBATCH -o ./log/Driverpower_%A-%a_out.txt
#SBATCH -e ./log/Driverpower_%A-%a_error.txt
#SBATCH -p amd-ep2,amd-ep2-short
#SBATCH --qos=huge
#SBATCH --array=1-1

## user's own commands below
ID=${1}
reg=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
rawDir="/path/PROJECT/05_SELECTION/01_XIWEI/my_resource/"
wkDir="/path/Driverpower/PCAWG"

cat ${rawDir}/gc19_pc.${reg}/gc19_pc.${reg}.excDupSite.bed | awk '{printf $1"\t"$2"\t"$3"\t"$4"\n"}' > ${wkDir}/annotations/gc19_pc.${reg}.excDupSite.bed

################################################################################################################
# step 2 : format the somatic mutations used in dNdS-Fun into the DriverPower
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=Driverpower
#SBATCH --ntasks-per-node=2
###SBATCH --mem=20G
#SBATCH -o ./log/Driverpower_%A-%a_out.txt
#SBATCH -e ./log/Driverpower_%A-%a_error.txt
#SBATCH -p amd-ep2,amd-ep2-short
#SBATCH --qos=huge
#SBATCH --array=1-25

## user's own commands below
ID=${1}
tissue=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
reg=${2}
rawDir="/path/PCAWG/singleproject_tissue/"
wkDir="/path/Driverpower/PCAWG/singleproject_tissue"

echo -e "#chr\tstart\tend\tref\talt\tdonor" > ${wkDir}/PCAWG_SNVs_allmutations.${tissue}.gc19_pc.${reg}.Driverpower.tsv
cat ${rawDir}/PCAWG_SNVs_allmutations.${tissue}.gc19_pc.${reg}.all.impScore.txt | awk '{printf $2"\t""%.f""\t""%.f""\t"$4"\t"$5"\t"$1"\n", $3-1,$3}' >> ${wkDir}/PCAWG_SNVs_allmutations.${tissue}.gc19_pc.${reg}.Driverpower.tsv

################################################################################################################
# step 3 : add the functional impact scores into the mutations
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=Driverpower
#SBATCH --ntasks-per-node=2
###SBATCH --mem=20G
#SBATCH -o ./log/Driverpower_%A-%a_out.txt
#SBATCH -e ./log/Driverpower_%A-%a_error.txt
#SBATCH -p amd-ep2,amd-ep2-short
#SBATCH --qos=huge
#SBATCH --array=1-25

## user's own commands below
ID=${1}
tissue=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
reg=${2}
rawDir="/path/PROJECT/05_SELECTION/07_Tri_QuantSet_excDupSite_NC/04_GeneSet/01_Resource/PCAWG/singleproject_tissue/"
wkDir="/path/Driverpower/PCAWG/singleproject_tissue"

cat ${rawDir}/PCAWG_SNVs_allmutations.${tissue}.gc19_pc.${reg}.all.impScore.txt | \
awk 'BEGIN {OFS="\t"} {print $2, $3-1, $3, $7}' | \
bedtools intersect -a ${wkDir}/../annotations/gc19_pc.${reg}.excDupSite.bed -b stdin -wa -wb > ${wkDir}/../CADD_anno/PCAWG_SNVs_allmutations.${tissue}.gc19_pc.${reg}.CADD_ele.tsv

printf "binID\tCADD\n" > ${wkDir}/../CADD_anno/PCAWG_SNVs_allmutations.${tissue}.gc19_pc.${reg}.CADD.tsv
bedtools groupby -i ${wkDir}/../CADD_anno/PCAWG_SNVs_allmutations.${tissue}.gc19_pc.${reg}.CADD_ele.tsv -g 4 -c 8 -o mean  >> ${wkDir}/../CADD_anno/PCAWG_SNVs_allmutations.${tissue}.gc19_pc.${reg}.CADD.tsv

################################################################################################################
# step 4 : Run the DriverPower
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=driverpower
#SBATCH --ntasks-per-node=2
###SBATCH --mem=200G
#SBATCH -o ./log/driverpower_%A-%a_out.txt
#SBATCH -e ./log/driverpower_%A-%a_error.txt
#SBATCH -p amd-ep2-short
#SBATCH --qos=normal
#SBATCH --array=1-1

## user's own commands below
pwd;hostname;date
#source ~/.bashrc
source /home/yangjianLab/zhengmengyue/miniconda3/bin/activate driverpower

## user's own commands below
ID=${1}
tissue=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
reg=${2}

wkDir="/path/Driverpower/PCAWG/"
prepare="/path/SOFTWARE.bak/DriverPower/DriverPower-1.0.2/prepare.py"
callable="/path/SOFTWARE.bak/DriverPower/DriverPower-1.0.2/data/callable.bed.gz"
PCAWG_elements="/path/SOFTWARE.bak/DriverPower/DriverPower-1.0.2/data/PCAWG_test_genomic_elements.bed12.gz"
train_element="/path/SOFTWARE.bak/DriverPower/Test/train_elements.tsv.gz"
test_element="${wkDir}/annotations/gc19_pc.${reg}.excDupSite.bed"
mut_fil="${wkDir}/singleproject_tissue/PCAWG_SNVs_allmutations.${tissue}.gc19_pc.${reg}.Driverpower.tsv"
train_y="/path/SOFTWARE.bak/DriverPower/Test/train_y.tsv"
test_y="${wkDir}/singleproject_tissue/PCAWG_SNVs_allmutations.${tissue}.gc19_pc.${reg}.Driverpower.test_y.tsv"

python ${prepare} ${mut_fil} ${test_element} ${callable} ${test_y}
mkdir -p ${wkDir}/output/${tissue}/${reg}
driverpower infer \
    --feature ${wkDir}/test_feature.hdf5 \
    --response ${test_y} \
    --model ./output/tutorial.GBM.model.pkl \
    --name 'DriverPower_burden_function' \
    --outDir ${wkDir}/output/${tissue}/${reg}/ \
    --funcScore ${wkDir}/CADD_anno/PCAWG_SNVs_allmutations.${tissue}.gc19_pc.${reg}.CADD.tsv \
    --funcScoreCut "CADD:0.01"


