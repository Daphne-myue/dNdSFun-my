################################################################################################################
# step 1 : remove the the mutations of the same number of the known gene sets to form the sampling excGene group
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=PCAWG
#SBATCH --ntasks-per-node=2
#SBATCH --mem=20G
#SBATCH -o ./log/excGene_Sampling_%A-%a_out.txt
#SBATCH -e ./log/excGene_Sampling_%A-%a_error.txt
#SBATCH -p intel-sc3,amd-ep2,amd-ep2-short
#SBATCH --qos=huge
#SBATCH --array=1-300

## user's own commands below
pwd;hostname;date

ID=$1
seed_num=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
reg=${2}
neutral=${3}
quant=${4}
gene_group=${5}
clonal_group=${6}
project_group=${7}
mut_resource=${8}
negbeta=${9}
model=${10}

module load R/4.0.5
wkDir="/path/"
gencode_anno="gc19_pc.${reg}"
mkdir -p ${wkDir}/09_excGene_Sampling_50/01_RandomIn/${mut_resource}/${clonal_group}/${gene_group}/
inputDir="${wkDir}/09_excGene_Sampling_50/01_RandomIn/${mut_resource}/${clonal_group}/${gene_group}/"
mkdir -p ${wkDir}/09_excGene_Sampling_50/02_RandomOut/${mut_resource}/${clonal_group}/${gene_group}/
outputDir="${wkDir}/09_excGene_Sampling_50/02_RandomOut/${mut_resource}/${clonal_group}/${gene_group}/"

threshold_fil="${wkDir}/QuantSet_${neutral}_${quant}/${gencode_anno}/threshold.txt"
threshold_positive=`cat ${threshold_fil} | cut -f 2 | head -n 1`
threshold_neutral=`cat ${threshold_fil} | cut -f 2 | tail -n 1`

input="${wkDir}/04_GeneSet/01_Resource/${mut_resource}/${clonal_group}/${mut_resource}_SNVs_allmutations.${project_group}.gc19_pc.${reg}_CADD_score.all.impScore.txt"
input_bed="${wkDir}/04_GeneSet/01_Resource/${mut_resource}/${clonal_group}/${mut_resource}_SNVs_allmutations.${project_group}.gc19_pc.${reg}_CADD_score.all.impScore.bed"
input_exc="${inputDir}/${mut_resource}_SNVs_allmutations.${project_group}.gc19_pc.${reg}_CADD_score.${gene_group}_exc.${seed_num}.impScore.txt"

genefil="/path/GeneList/gencode.v19.${gene_group}.bed"
allgenefil="/path/GeneList/gencode.v19.gene.noXY.bed"
genefil_sampling="${inputDir}/gencode.v19.${project_group}.${seed_num}.bed"

line_num=`cat ${genefil} | grep -v "chrX" | grep -v "chrY" | grep -v "chrM" | wc -l | cut -d " " -f 1`
shuf -n${line_num} ${allgenefil} -o ${genefil_sampling}
bedtools intersect -v -a ${input_bed} -b ${genefil_sampling} -wa | awk '{printf $4"\t"$5"\t""%.f""\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\n", $6}' | sort | uniq > ${input_exc}


################################################################################################################
# step 2 : Run dNdS-Fun of the sampling excGene profiles
# use the seed_num to generate different random mutations, and replicated about 5,000 times
################################################################################################################
##CADD method  -  geneset_exc
R_command="/path/src/cal_CADD_dndsWGS.NEG.R"
mutsFile=${input_exc}
refDb_element="${wkDir}/QuantSet_${neutral}_${quant}/${gencode_anno}/CADD_score/${gencode_anno}_CADD_score_ALL_${quant}.rda"
globaldnds_outFile="${outputDir}/CADD_globadnds.${mut_resource}.${reg}.${gene_group}.${clonal_group}.${project_group}.${gene_group}_exc.${seed_num}.all_genes.out"
dndsloc_outFile="${outputDir}/CADD_dndsloc.${mut_resource}.${reg}.${gene_group}.${clonal_group}.${project_group}.${gene_group}_exc.${seed_num}.single_genes.out"
dndscv_outFile="${outputDir}/CADD_dndscv.${mut_resource}.${reg}.${gene_group}.${clonal_group}.${project_group}.${gene_group}_exc.${seed_num}.single_genes.out"
iscv="nocv"

Rscript ${R_command} ${mutsFile} ${refDb_element} ${reg} ${neutral} ${quant} ${threshold_positive} ${threshold_neutral} ${gene_group}_exc ${clonal_group} ${project_group} ${mut_resource} ${globaldnds_outFile} ${dndsloc_outFile} ${dndscv_outFile} ${negbeta} ${iscv} ${model}










