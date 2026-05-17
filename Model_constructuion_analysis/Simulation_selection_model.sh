################################################################################################################
# step 1 : Split the trinucleotide pools into 2 parts (less functiona and more functional parts)
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=Mimic
#SBATCH --ntasks-per-node=2
#SBATCH --mem=17G
#SBATCH -o ./log3/cal_globaldnds_%A-%a_out.txt
#SBATCH -e ./log3/cal_globaldnds_%A-%a_error.txt
#SBATCH -p intel-debug
#SBATCH --qos=debug
#SBATCH --array=1-50

## user's own commands below
pwd;hostname;date

dichtomy_threshold=$1
reg=$2
cat gc19_pc.${reg}_CADD_score_ALL.impScore.txt | awk '{if(let $6<=${dichtomy_threshold}){print $0}}' > gc19_pc.${reg}_CADD_score_ALL.Less.impScore.txt
cat gc19_pc.${reg}_CADD_score_ALL.impScore.txt | awk '{if(let $6>${dichtomy_threshold}){print $0}}' > gc19_pc.${reg}_CADD_score_ALL.More.impScore.txt

################################################################################################################
# step 2 : Randomly extract a half number of mutations from less functional mutations
# Randomly extract a half number of mutations from more functional mutations
# Duplicate the differennt proportion of mutations in more functional mutations (Positive selection)
# Discard the differennt proportion of mutations in more functional mutations (Negative selection)
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=Mimic
#SBATCH --ntasks-per-node=2
#SBATCH --mem=17G
#SBATCH -o ./log3/cal_globaldnds_%A-%a_out.txt
#SBATCH -e ./log3/cal_globaldnds_%A-%a_error.txt
#SBATCH -p intel-debug
#SBATCH --qos=debug
#SBATCH --array=1-50

## user's own commands below
pwd;hostname;date

ID=$1
seed_num=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
reg=${2}
total_num=${3}
variant_name=${4}

let "neutral_num=${total_num}*5/10"
let "select_num=${total_num}*5/10"
echo ${neutral_num}
echo ${select_num}

wkDir="/path/"
mkdir -p ${wkDir}/06_MimcSelection_50/01_MimcInput/${reg}/
mkdir -p ${wkDir}/06_MimcSelection_50/01_MimcOutput/${reg}/

DOWNFILE="${wkDir}/QuantSet_50_0.50/gc19_pc.${reg}/CADD_score/gc19_pc.${reg}_CADD_score_ALL.Less.impScore.txt"
UPFILE="${wkDir}/QuantSet_50_0.50/gc19_pc.${reg}/CADD_score/gc19_pc.${reg}_CADD_score_ALL.More.impScore.txt"

SamplingFile_neutral="${wkDir}/06_MimcSelection_50/01_MimcInput/${reg}/gc19_pc.${reg}_CADD_score.${variant_name}.neutral.${seed_num}.impScore.txt"
SamplingFile_select="${wkDir}/06_MimcSelection_50/01_MimcInput/${reg}/gc19_pc.${reg}_CADD_score.${variant_name}.select.${seed_num}.impScore.txt"
SamplingFile_input="${wkDir}/06_MimcSelection_50/01_MimcInput/${reg}/gc19_pc.${reg}_CADD_score.${variant_name}.NoSelect.${seed_num}.impScore.txt"
SamplingFile_select_prefix="${wkDir}/06_MimcSelection_50/01_MimcInput/${reg}/gc19_pc.${reg}_CADD_score.${variant_name}"

SamplingFile_select_tmp2="${wkDir}/06_MimcSelection_50/01_MimcInput/${reg}/gc19_pc.${reg}_CADD_score.${variant_name}.${seed_num}.tmp2"
SamplingFile_select_tmp3="${wkDir}/06_MimcSelection_50/01_MimcInput/${reg}/gc19_pc.${reg}_CADD_score.${variant_name}.${seed_num}.tmp3"

shuf -n${neutral_num} ${DOWNFILE} -o ${SamplingFile_neutral}
shuf -n${select_num} ${UPFILE} -o ${SamplingFile_select}
cat ${SamplingFile_neutral} ${SamplingFile_select} > ${SamplingFile_input}

for rate in 5 10 15 20 25 30 ## set the different proportion of mutations to duplicate or discard
do
  let "select_num2=${select_num}/100*${rate}"
  let "select_num3=${select_num}-${select_num2}"

  shuf -n${select_num2} ${UPFILE} -o ${SamplingFile_select_tmp2}
  shuf -n${select_num3} ${UPFILE} -o ${SamplingFile_select_tmp3}

  cat ${SamplingFile_input}  ${SamplingFile_select_tmp2} > ${SamplingFile_select_prefix}.PosSelect_${rate}.${seed_num}.impScore.txt
  cat ${SamplingFile_neutral}  ${SamplingFile_select_tmp3} > ${SamplingFile_select_prefix}.NegSelect_${rate}.${seed_num}.impScore.txt
done


################################################################################################################
# step 3 : Run dNdS-Fun method
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=CADD_dnds
#SBATCH --ntasks-per-node=2
#SBATCH --mem=30G
#SBATCH -o ./log3/CADD_dnds_%A-%a_out.txt
#SBATCH -e ./log3/CADD_dnds_%A-%a_error.txt
#SBATCH -p intel-e5,amd-ep2,amd-ep1
#SBATCH --qos=normal
#SBATCH --array=1-1

## user's own commands below
#conda deactivate
pwd;hostname;date

ID=$1
seed_num=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
reg=${2}
select_group=${3}
variant_name=${4}
neutral=50
quant="0.50"

wkDir="/path/"
gencode_anno="gc19_pc.${reg}"

threshold_fil="${wkDir}/QuantSet_${neutral}_${quant}/${gencode_anno}/threshold.txt"
threshold_positive=`cat ${threshold_fil} | cut -f 2 | head -n 1`
threshold_neutral=`cat ${threshold_fil} | cut -f 2 | tail -n 1`

##files
R_command="/path/src/cal_CADD_dndsWGS.NEG.R"
mutsFile="${wkDir}/06_MimcSelection_50/01_MimcInput/${reg}/gc19_pc.${reg}_CADD_score.${variant_name}.${select_group}.${seed_num}.impScore.txt"
refDb_element="${wkDir}/QuantSet_${neutral}_${quant}/${gencode_anno}/CADD_score/${gencode_anno}_CADD_score_ALL_${quant}.rda"
globaldnds_outFile="${wkDir}/06_MimcSelection_50/01_MimcOutput/${reg}/CADD_globadnds.${reg}_CADD_score.${variant_name}.${select_group}.${seed_num}.all_genes.out"
dndsloc_outFile="${wkDir}/06_MimcSelection_50/01_MimcOutput/${reg}/CADD_dndsloc.${reg}_CADD_score.${variant_name}.${select_group}.${seed_num}.single_genes.out"
dndscv_outFile="${wkDir}/06_MimcSelection_50/01_MimcOutput/${reg}/CADD_dndscv.${reg}_CADD_score.${variant_name}.${select_group}.${seed_num}.single_genes.out"
gene_group="MimcSelect"
clonal_group="MimcSelect"
project_group="MimcSelect"
mut_resource="MimcSelect"
negbeta="1"
iscv="nocv"
model="1"

Rscript ${R_command} ${mutsFile} ${refDb_element} ${reg} ${neutral} ${quant} ${threshold_positive} ${threshold_neutral} ${gene_group} ${clonal_group} ${project_group} ${mut_resource} ${globaldnds_outFile} ${dndsloc_outFile} ${dndscv_outFile} ${negbeta} ${iscv} ${model}
