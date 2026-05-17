################################################################################################################
# step 1 : format the genomic annotations used in dNdS-Fun into the Larva
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=LARVA
#SBATCH --ntasks-per-node=2
###SBATCH --mem=20G
#SBATCH -o ./log/LARVA_%A-%a_out.txt
#SBATCH -e ./log/LARVA_%A-%a_error.txt
#SBATCH -p amd-ep2,amd-ep2-short
#SBATCH --qos=huge
#SBATCH --array=1-11

## user's own commands below
ID=${1}
reg=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
rawDir="/path/my_resource/"
wkDir="/path/LARVA/PCAWG"

cat ${rawDir}/gc19_pc.${reg}/gc19_pc.${reg}.excDupSite.bed | awk '{printf $1"\t"$2"\t"$3"\t"$4"\n"}' > ${wkDir}/annotations/gc19_pc.${reg}.bed

################################################################################################################
# step 2 : format the somatic mutations used in dNdS-Fun into the Larva
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=larva
#SBATCH --ntasks-per-node=2
#SBATCH --mem=20G
#SBATCH -c 8
#SBATCH -o ./log/larva_%A-%a_out.txt
#SBATCH -e ./log/larva_%A-%a_error.txt
#SBATCH -p amd-ep2,amd-ep2-short
#SBATCH --qos=huge
#SBATCH --array=1-25

## user's own commands below
ID=${1}
tissue=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`

cat /path/singleproject_tissue/PCAWG_SNVs_allmutations.${tissue}.txt | awk '{printf $1"\t""%.f""\t""%.f""\t"$8"\t"$6"\n", $2-1,$2}' > /path/LARVA/PCAWG/singleproject_tissue/PCAWG_SNVs_allmutations.${tissue}.larva.txt

cat /path/PCAWG_SNVs_allmutations.txt | awk '{printf $1"\t""%.f""\t""%.f""\t"$8"\t"$6"\n", $2-1,$2}' > /path/LARVA/PCAWG/singleproject_tissue/PCAWG_SNVs_allmutations.all.larva.txt

################################################################################################################
# step 3 : Run larva
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=larva
#SBATCH --ntasks-per-node=2
#SBATCH --mem=20G
#SBATCH -c 8
#SBATCH -o ./log/larva_%A-%a_out.txt
#SBATCH -e ./log/larva_%A-%a_error.txt
#SBATCH -p amd-ep2,amd-ep2-short
#SBATCH --qos=huge
#SBATCH --array=1-25

## user's own commands below
ID=${1}
tissue=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
bed=$2

variant_file="singleproject_tissue/PCAWG_SNVs_allmutations.${tissue}.larva.txt"
annotation_prefix=${bed}
out_file=${tissue}.${annotation_prefix}.result.txt
annotation_file="/path/SOFTWARE.bak/LARVA/resource/noncoding-annotations/${annotation_prefix}.bed"
echo "larva -vf ${variant_file} -af ${annotation_file} -o singleproject_tissue_output/${out_file}


