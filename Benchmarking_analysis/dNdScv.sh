################################################################################################################
# step 1 : split the samples of TCGA and PCAWG into difference cancer subtypes
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=split
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH -o ./log/split_%A-%a_out.txt
#SBATCH -e ./log/split_%A-%a_error.txt
#SBATCH -p intel-sc3,amd-ep2,amd-ep2-short
#SBATCH --qos=huge
#SBATCH --array=1-32

ID=$1
tissue=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
grep -w ${tissue} singleproject_single/TCGA_SNVs_allmutations.bed > singleproject_single/TCGA_SNVs_allmutations.${tissue}.bed
grep -w ${tissue} singleproject_single/PCAWG_SNVs_allmutations.bed > singleproject_single/PCAWG_SNVs_allmutations.${tissue}.bed

################################################################################################################
# step 2 : Run the dNdScv
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=Sampling
#SBATCH --ntasks-per-node=2
#SBATCH --mem=10G
#SBATCH -o ./log/cal_globaldnds_%A-%a_out.txt
#SBATCH -e ./log/cal_globaldnds_%A-%a_error.txt
#SBATCH -p intel-sc3,amd-ep2,amd-ep2-short
#SBATCH --qos=huge
#SBATCH --array=1-1

## user's own commands below
##ID=$1
##id=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
reg="cds-exon"
project_group=$1
gene_group=$2
clonal_group=$3
mut_resource=$4
opt_num=$5

module load R/4.0.5
export R_LIBS_USER="/path/SOFTWARE.bak/R_LIB_4.0.5"
wkDir="/path/"
mkdir -p ${wkDir}/04_GeneSet/02_dndsOut/${mut_resource}/${clonal_group}/

R_command="${wkDir}/dndscv.R"
mutsFile="${wkDir}/singleproject_single/PCAWG_SNVs_allmutations.${tissue}.bed"
refDb_element="/path/dndscv-0.1.0/data/refcds_hg19.rda"
substmodel="/path/dndscv-0.1.0/data/submod_192r_3w.rda"
cancergene="/path/dndscv-0.1.0/data/cancergenes_cgc81.rda"
cv_hg19="/path/dndscv-0.1.0/data/covariates_hg19.rda"
globaldnds_outFile="${wkDir}/04_GeneSet/02_dndsOut/${mut_resource}/${clonal_group}/refCDS_globadnds.${mut_resource}.${reg}.${gene_group}.${clonal_group}.${project_group}.all_genes.qua.out"
dndsloc_outFile="${wkDir}/04_GeneSet/02_dndsOut/${mut_resource}/${clonal_group}/refCDS_dndsloc.${mut_resource}.${reg}.${gene_group}.${clonal_group}.${project_group}.single_genes.out"
dndscv_outFile="${wkDir}/04_GeneSet/02_dndsOut/${mut_resource}/${clonal_group}/refCDS_dndscv.${mut_resource}.${reg}.${gene_group}.${clonal_group}.${project_group}.single_genes.out"
cv_hg19="NULL"

##command
Rscript ${R_command} ${mutsFile} ${refDb_element} ${substmodel} ${cancergene} ${cv_hg19} ${project_group} ${clonal_group} ${gene_group} ${globaldnds_outFile} ${dndsloc_outFile} ${dndscv_outFile} ${mut_resource} ${opt_num}

