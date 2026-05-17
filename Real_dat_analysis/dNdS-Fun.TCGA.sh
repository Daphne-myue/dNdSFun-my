################################################################################################################
# step 1 : split the samples of TCGA data into difference cancer subtypes
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
grep -w ${tissue} singleproject_single/TCGA_SNVs_allmutations.bed > singleproject_single/TCGA_SNVs_allmutations.${tissue}.bed

################################################################################################################
# step 2 : R code to add the impact score into the mutations
################################################################################################################
rm(list = ls())
library("data.table")
args <- commandArgs(trailingOnly = TRUE)
impScore_file <- args[1]
input <- args[2]
output <- args[3]

# Load and process impScore data
impScore_dat <- as.data.frame(fread(impScore_file))
colnames(impScore_dat) <- c("CHR.1", "START.1", "END.1", "Gene", "Chrom", "Pos", "Ref", "Alt", 
                            "RawScore", "PHRED", "RawScore_rankscore")
# Create a unique key for each row based on Chromosome, Position, Reference, and Alternate alleles
impScore_dat$KEY <- paste0(impScore_dat$Chrom, impScore_dat$Pos, impScore_dat$Ref, impScore_dat$Alt)

# Load and process mutation data
mutation_dat <- as.data.frame(fread(input))
colnames(mutation_dat) <- c("CHR", "START", "END", "REF", "ALT", "SAMPLE", "FLANK", "tmp", 
                            "PROJECT", "DONOR", "Read_count", "ALT_count", "REF_count", "vaf")
# Create a unique key for each mutation entry
mutation_dat$KEY <- paste0(mutation_dat$CHR, mutation_dat$END, mutation_dat$REF, mutation_dat$ALT)

# Convert data frames to data tables and merge them on the "KEY" column
dt1 <- data.table(impScore_dat, key = "KEY")
dt2 <- data.table(mutation_dat, key = "KEY")
dt <- as.data.frame(merge(dt2, dt1))

# Select relevant columns for output
dtres <- dt[, c("CHR", "START", "END", "DONOR", "CHR", "END", "REF", "ALT", 
                "RawScore", "PHRED", "RawScore_rankscore")]

# Write the results to a file without quotes, column names, or row names
fwrite(dtres, output, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


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

## user's own commands below
#conda deactivate
pwd;hostname;date

reg=$1
##reg=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
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

## add functional impact scores into data
for chr in {1..22}
do
   Rscript S00_Add_impScore.TCGA.R /path/CADD_Score/${gencode_anno}/${gencode_anno}_CADD_score.impScore.${chr}.bed splitchr/TCGA_SNVs_allmutations.${chr}.txt splitchr/TCGA_SNVs_allmutations.txt.${gencode_anno}_CADD_score.impScore.${chr}.txt
done

cat splitchr/TCGA_SNVs_allmutations.txt.${gencode_anno}_CADD_score.impScore.${chr}.txt > splitchr/TCGA_SNVs_allmutations.txt.${gencode_anno}_CADD_score.impScore.all.txt

threshold_fil="${runDir}/QuantSet_${neutral}_${quant}/${gencode_anno}/threshold.txt"
threshold_positive=`cat ${threshold_fil} | cut -f 2 | head -n 1`
threshold_neutral=`cat ${threshold_fil} | cut -f 2 | tail -n 1`

## run the dNdS-Fun
R_command="/path/src/cal_CADD_dndsWGS.NEG.R"
mutsFile="${runDir}/04_GeneSet/01_Resource/${mut_resource}/${clonal_group}/${mut_resource}_SNVs_allmutations.${project_group}.${gencode_anno}_CADD_score.${gene_group}.impScore.txt"
refDb_element="${runDir}/QuantSet_${neutral}_${quant}/${gencode_anno}/CADD_score/${gencode_anno}_CADD_score_ALL_${quant}.rda"
globaldnds_outFile="${runDir}/04_GeneSet_${neutral}/02_dndsOut/${mut_resource}/${clonal_group}/CADD_globadnds.${mut_resource}.${reg}.${gene_group}.${clonal_group}.${project_group}.all_genes.negbeta_${negbeta}.out"
cv_hg19="/storage/yangjianLab/zhengmengyue/PROJECT/05_SELECTION/06_RefCDS/dndscv-0.1.0/data/covariates_hg19.rda"
dndsloc_outFile="${runDir}/04_GeneSet_${neutral}/02_dndsOut/${mut_resource}/${clonal_group}/CADD_dndsloc.${mut_resource}.${reg}.${gene_group}.${clonal_group}.${project_group}.single_genes.negbeta_${negbeta}.out"
dndscv_outFile="${runDir}/04_GeneSet_${neutral}/02_dndsOut/${mut_resource}/${clonal_group}/CADD_dndscv.${mut_resource}.${reg}.${gene_group}.${clonal_group}.${project_group}.single_genes.negbeta_${negbeta}.out"
iscv="nocv"

Rscript ${R_command} ${mutsFile} ${refDb_element} ${reg} ${neutral} ${quant} ${threshold_positive} ${threshold_neutral} ${gene_group} ${clonal_group} ${project_group} ${mut_resource} ${globaldnds_outFile} ${dndsloc_outFile} ${dndscv_outFile} ${negbeta} ${iscv} ${model}














