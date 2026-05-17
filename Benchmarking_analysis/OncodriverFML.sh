################################################################################################################
# step 1 : format the genomic annotations used in dNdS-Fun into the OncodriveFML
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=OncodriverFML
#SBATCH --ntasks-per-node=2
###SBATCH --mem=20G
#SBATCH -o ./log/OncodriverFML_%A-%a_out.txt
#SBATCH -e ./log/OncodriverFML_%A-%a_error.txt
#SBATCH -p amd-ep2,amd-ep2-short
#SBATCH --qos=huge
#SBATCH --array=1-11

## user's own commands below
ID=${1}
reg=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
echo -e "CHROMOSOME\tSTART\tEND\tSTRAND\tELEMENT\tSEGMENT\tSYMBOL" > /path/OncodriverFML/PCAWG/annotations/${reg}.tsv
cat /path/my_resource/gc19_pc.${reg}/gc19_pc.${reg}.excDupSite.bed |sed 's/^chr//g'| awk '{printf $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$4"\t"$4"\n"}' >> /path/OncodriverFML/PCAWG/annotations/${reg}.tsv 
gzip /path/OncodriverFML/PCAWG/annotations/${reg}.tsv

################################################################################################################
# step 2 : format the somatic mutations used in dNdS-Fun into the OncodriveFML
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=OncodriverFML
#SBATCH --ntasks-per-node=2
###SBATCH --mem=20G
#SBATCH -o ./log/OncodriverFML_%A-%a_out.txt
#SBATCH -e ./log/OncodriverFML_%A-%a_error.txt
#SBATCH -p amd-ep2,amd-ep2-short
#SBATCH --qos=huge
#SBATCH --array=1-25

## user's own commands below
ID=${1}
tissue=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`

echo -e "CHROMOSOME\tPOSITION\tREF\tALT\tSAMPLE" > /path/OncodriverFML/PCAWG/singleproject_tissue/PCAWG_SNVs_allmutations.${tissue}.OncodriverFML.txt
cat /path/PCAWG/singleproject_tissue/PCAWG_SNVs_allmutations.${tissue}.txt | sed 's/^chr//g'| awk 'sub(/^"chr"/,"");{printf $1"\t""%.f""\t"$4"\t"$5"\t"$6"\n", $2}' >> /path/OncodriverFML/PCAWG/singleproject_tissue/PCAWG_SNVs_allmutations.${tissue}.OncodriverFML.txt

################################################################################################################
# step 3 : Run oncodriveFML
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=OncodriverFML
#SBATCH --ntasks-per-node=2
###SBATCH --mem=200G
#SBATCH -o ./log/OncodriverFML_%A-%a_out.txt
#SBATCH -e ./log/OncodriverFML_%A-%a_error.txt
#SBATCH -p amd-ep2-short
#SBATCH --qos=normal
#SBATCH --array=1-1

## user's own commands below
pwd;hostname;date
source ~/.bashrc
source /home/yangjianLab/zhengmengyue/miniconda3/bin/activate OncodriverFML
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

ID=${1}
tissue=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
reg=${2}
mkdir output/${tissue}/${reg}

variants_file="/path/OncodriverFML/PCAWG/singleproject_tissue/PCAWG_SNVs_allmutations.${tissue}.OncodriverFML.txt"
region_file="/path/OncodriverFML/PCAWG/annotations/${reg}.tsv.gz"
signature_correlation="wg"
oncodrivefml -i ${variants_file} -e ${region_file} --signature-correction ${signature_correlation} --seed 123 --force --output output/${tissue}/${reg}/