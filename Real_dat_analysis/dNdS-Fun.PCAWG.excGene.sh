################################################################################################################
# step 1 : remove the the mutations of known gene sets to form the excGene group
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=PCAWG
#SBATCH --ntasks-per-node=4
#SBATCH --mem=30G
#SBATCH -o ./log/PickGene_%A-%a_out.txt
#SBATCH -e ./log/PickGene_%A-%a_error.txt
#SBATCH -p intel-sc3,amd-ep2,amd-ep2-short
#SBATCH --qos=huge
#SBATCH --array=1-8

## user's own commands below
pwd;hostname;date

ID=$1
reg=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
mut_resource=$2
clonal_group=$3
project_group=$4
geneset=$5

wkDir="/path/"
input="${wkDir}/04_GeneSet/01_Resource/${mut_resource}/${clonal_group}/${mut_resource}_SNVs_allmutations.${project_group}.gc19_pc.${reg}_CADD_score.all.impScore.txt"
input_bed="${wkDir}/04_GeneSet/01_Resource/${mut_resource}/${clonal_group}/${mut_resource}_SNVs_allmutations.${project_group}.gc19_pc.${reg}_CADD_score.all.impScore.bed"
output_inc="${wkDir}/04_GeneSet/01_Resource/${mut_resource}/${clonal_group}/${mut_resource}_SNVs_allmutations.${project_group}.gc19_pc.${reg}_CADD_score.${geneset}_inc.impScore.txt"
output_exc="${wkDir}/04_GeneSet/01_Resource/${mut_resource}/${clonal_group}/${mut_resource}_SNVs_allmutations.${project_group}.gc19_pc.${reg}_CADD_score.${geneset}_exc.impScore.txt"
genefil="/path/GeneList/gencode.v19.${geneset}.bed"

cat ${input} | grep -v ^"V1" | awk '{printf $2"\t""%.f""\t""%.f""\t"$0"\n", $3-1,$3}' > ${input_bed}
bedtools intersect -v -a ${input_bed} -b ${genefil} -wa | awk '{printf $4"\t"$5"\t""%.f""\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\n", $6}' | sort | uniq > ${output_exc}

################################################################################################################
# step 2 : To run this script for each gene set among all cancer subtypes
################################################################################################################
rm run_09_PickGenes.sh
for tissue in `cat PCAWGALL_tissue.meta.list`
do for geneset in `cat /path/04_GeneSet/01_Resource/GeneList/excGene.list`
do
  echo "sbatch S33_PickGeneGroup.PCAWG.bedtools.sh gc19_pc.excDupSite.les.list PCAWGALL singleproject_tissue ${tissue} ${geneset}" >> run_09_PickGenes.sh
  echo "sbatch S33_PickGeneGroup.PCAWG.bedtools.sh gc19_pc.excDupSite.les.list PCAWGALL allproject all ${geneset}" >> run_09_PickGenes.sh
done
done

################################################################################################################
# step 3 : Run the dNdS-Fun
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=CADD_t
###SBATCH --ntasks-per-node=6
#SBATCH --mem=50G
#SBATCH -c 20
#SBATCH -o ./log5/CADD_%A-%a_out.txt
#SBATCH -e ./log5/CADD_%A-%a_error.txt
#SBATCH -p intel-sc3,amd-ep2,amd-ep2-short
#SBATCH --qos=normal
#SBATCH --array=1-1

pwd;hostname;date

reg=$1
neutral=$2
quant=$3
gene_group=$4
clonal_group=$5
project_group=$6
mut_resource=$7
negbeta=$8
model=$9

module load R/4.0.5
runDir="/path/"
gencode_anno="gc19_pc.${reg}"
mkdir -p ${runDir}/04_GeneSet_${neutral}/02_dndsOut/${mut_resource}/${clonal_group}

threshold_fil="${runDir}/QuantSet_${neutral}_${quant}/${gencode_anno}/threshold.txt"
threshold_positive=`cat ${threshold_fil} | cut -f 2 | head -n 1`
threshold_neutral=`cat ${threshold_fil} | cut -f 2 | tail -n 1`

##files
R_command="/path/src/cal_CADD_dndsWGS.NEG.R"
mutsFile="${runDir}/04_GeneSet/01_Resource/${mut_resource}/${clonal_group}/${mut_resource}_SNVs_allmutations.${project_group}.${gencode_anno}_CADD_score.${gene_group}.impScore.txt"
refDb_element="${runDir}/QuantSet_${neutral}_${quant}/${gencode_anno}/CADD_score/${gencode_anno}_CADD_score_ALL_${quant}.rda"
globaldnds_outFile="${runDir}/04_GeneSet_${neutral}/02_dndsOut/${mut_resource}/${clonal_group}/CADD_globadnds.${mut_resource}.${reg}.${gene_group}.${clonal_group}.${project_group}.all_genes.negbeta_${negbeta}.out"
dndsloc_outFile="${runDir}/04_GeneSet_${neutral}/02_dndsOut/${mut_resource}/${clonal_group}/CADD_dndsloc.${mut_resource}.${reg}.${gene_group}.${clonal_group}.${project_group}.single_genes.negbeta_${negbeta}.out"
dndscv_outFile="${runDir}/04_GeneSet_${neutral}/02_dndsOut/${mut_resource}/${clonal_group}/CADD_dndscv.${mut_resource}.${reg}.${gene_group}.${clonal_group}.${project_group}.single_genes.negbeta_${negbeta}.out"
iscv="nocv"

##command
Rscript ${R_command} ${mutsFile} ${refDb_element} ${reg} ${neutral} ${quant} ${threshold_positive} ${threshold_neutral} ${gene_group} ${clonal_group} ${project_group} ${mut_resource} ${globaldnds_outFile} ${dndsloc_outFile} ${dndscv_outFile} ${negbeta} ${iscv} ${model}

date

################################################################################################################
# step 4 : Run the dNdS-Fun for each gene set among all cancer subtypes
################################################################################################################
rm run_05_CADD_dnds.PCAWGALL.excGene.multi.sh
for reg in `cat ../gc19_pc.excDupSite.les.list`
do for geneset in `cat /path/excGene.list`
do for project in `cat 04_GeneSet/01_Resource/PCAWGALL/PCAWGALL_tissue.meta.list`
do
  echo "sbatch S05_CADD_dnds_3models.NEG.neutral.excGene.sh ${reg} 50 0.50 ${geneset}_exc singleproject_tissue ${project} PCAWGALL 1 1" >> run_05_CADD_dnds.PCAWGALL.excGene.multi.sh
done
done
done











