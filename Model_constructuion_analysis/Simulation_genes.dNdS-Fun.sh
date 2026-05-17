################################################################################################################
# Step 1 : download the COSMIC somatic mutations and use liftover to tranlate into GRCh37 version
# Then selection the mutations of CGC genes as the simulated pool
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=split
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH -o ./log/split_%A-%a_out.txt
#SBATCH -e ./log/split_%A-%a_error.txt
#SBATCH -p intel-sc3,amd-ep2,amd-ep2-short
#SBATCH --qos=normal
#SBATCH --array=1-22

ID=$1
chr=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`

liftOver ../02_out/CosmicMutantExportCensus.GRCh38.chr${chr}.tsv /storage/yangjianLab/zhengmengyue/PROJECT/05_SELECTION/15_FunctionalScores/00_Resource/FAVOR/hg38ToHg19.over.chain.gz ../02_out/CosmicMutantExportCensus.GRCh37.chr${chr}.tsv ../02_out/CosmicMutantExportCensus.Unmap.chr${chr}.tsv
cat ../02_out/CosmicMutantExportCensus.GRCh37.chr${chr}.tsv > ../02_out/CosmicMutantExportCensus.GRCh37.ALL.tsv

for gene in `cat CGCgene.list`
do 
  grep -w "${gene}" ../02_out/CosmicMutantExportCensus.GRCh37.ALL.tsv >> ../02_out/CosmicMutantExportCensus.GRCh37.CGC.tsv
done

################################################################################################################
# Step 2 : Randomly selected the mutations of COSMIC CGC genes into the simulated data
################################################################################################################
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
randomfil <- args[1]
prefix <- args[2]
genelog <- args[3]
numtag <- args[4]
genenum <- as.numeric(args[5])
seed <- as.numeric(args[6])
set.seed(seed)

# Load random mutation data
random_dat <- as.data.frame(fread(randomfil, header = TRUE))
colnames(random_dat) <- c("donor", "chr", "pos", "ref", "alt", "RawScore", "PHRED", "RawScore_rankscore")

# Load and select genes for driver mutation simulation
Genefil <- "/path/CGCgene.list"  # Path to candidate gene list
Genedat <- as.data.frame(fread(Genefil, header = FALSE))
totGenenum <- length(Genedat$V1)
# Randomly select specified number of genes from candidate list
GenePick <- as.data.frame(matrix(Genedat[sample(1:totGenenum, genenum, replace = FALSE), ], ncol = 1))
fwrite(GenePick, genelog, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Define driver mutation numbers to iterate through
drivernum_list <- c(5, 10, 15, 20, 25, 30, 50, 100, 200, 500, 1000)
Genedat <- as.data.frame(fread(genelog, header = FALSE))
Genelist <- Genedat$V1

# For each specified driver mutation count, generate simulated mutations for each gene
for (x2 in seq_along(drivernum_list)) {
  drivernum <- drivernum_list[x2]
  tmp_res <- data.frame()
  
  # Process each gene in selected gene list
  for (x1 in seq_along(Genelist)) {
    Gene <- Genelist[x1]
    tmpdat <- as.data.frame(fread(paste0("/path/", Gene, ".GRCh37.CADD.functional.tsv")))
    tmp_dat_num <- nrow(tmpdat)
    
    # Randomly sample one mutation from gene data and replicate it to create driver mutations
    tmp_dat_tmp <- tmpdat[sample(1:tmp_dat_num, 1, replace = TRUE), ]
    tmp_dat_pick <- tmp_dat_tmp[rep(1, drivernum), ]
    colnames(tmp_dat_pick) <- c("donor", "chr", "pos", "ref", "alt", "RawScore", "PHRED", "RawScore_rankscore", "Gene")
    
    # Finalize output by renaming donors and selecting specific columns
    tmp_dat_out <- tmp_dat_pick[, c("donor", "chr", "pos", "ref", "alt", "RawScore", "PHRED", "RawScore_rankscore")]
    tmp_dat_out$donor <- paste0("donor", seq_len(nrow(tmp_dat_out)))
    tmp_res <- rbind(tmp_dat_out, tmp_res)
  }
  
  # Combine random data with driver mutations and save to output file
  tmp_res <- as.data.frame(tmp_res)
  tmp_out <- rbind(random_dat, tmp_res)
  fwrite(tmp_out, paste0(prefix, ".", numtag, ".", genenum, ".", drivernum, ".", seed, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

################################################################################################################
# step 3 : Randomly sampling the mutations multiple times
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=excDupSite
#SBATCH --mem=25G
#SBATCH -o ./log5/make_%A-%a_out.txt
#SBATCH -e ./log5/make_%A-%a_error.txt
#SBATCH -p intel-sc3,amd-ep2,amd-ep2-short
#SBATCH --qos=normal
#SBATCH --array=1-3

# User's commands
pwd; hostname; date

module load R/4.0.5
export R_LIBS_USER="/storage/yangjianLab/zhengmengyue/SOFTWARE.bak/R_LIB_4.0.5"
ID="$1"
id=$(head -n "${SLURM_ARRAY_TASK_ID}" "$ID" | tail -n1 | awk '{print $0}')
seed="${id}"
reg="$2"
clonal_group="$3"
project_group="$4"
numtag="$5"
genenum="$6"

GeneList="/path/01_scripts/CGCgene.list"
GenePath="/path/02_out/CGCgene/TP53.GRCh37.CADD.tsv"
GenePath="/path/02_out/CGCgene/TP53.GRCh37.CADD.functional.tsv"

realpath="/path/PCAWG/"
runDir="/path/"
inputDir="/path/"
mkdir -p "${inputDir}/04_simulateInput/PCAWG/gc19_pc.${reg}/${numtag}_${genenum}/"
realfil="${realpath}/${clonal_group}/PCAWG_SNVs_allmutations.${project_group}.gc19_pc.${reg}_CADD_score.all.impScore.txt"

for chr in {1..22}; do 
  reffil="${runDir}/02_CADDPool/gc19_pc.${reg}/gc19_pc.${reg}_CADD_trinucl_chr${chr}.txt"
  output="${inputDir}/04_simulateInput/PCAWG/gc19_pc.${reg}/${numtag}_${genenum}/PCAWG.${reg}.${numtag}_${genenum}.${seed}.chr${chr}.txt"
  /soft/devtools/R/R-4.0.5_installation/bin/Rscript S03_SimData.R "${realfil}" "${reffil}" "${output}" "${numtag}" "${chr}" "${seed}"
done

randomfil="${inputDir}/04_simulateInput/PCAWG/gc19_pc.${reg}/${numtag}_${genenum}/PCAWG.${reg}.${numtag}_${genenum}.${seed}.ALL.txt"
cat "${inputDir}/04_simulateInput/PCAWG/gc19_pc.${reg}/${numtag}_${genenum}/PCAWG.${reg}.${numtag}_${genenum}.${seed}.chr"*.txt > "${randomfil}"
rm "${inputDir}/04_simulateInput/PCAWG/gc19_pc.${reg}/${numtag}_${genenum}/PCAWG.${reg}.${numtag}_${genenum}.${seed}.chr"*.txt

genelog="${inputDir}/04_simulateInput/PCAWG/gc19_pc.${reg}/${numtag}_${genenum}/CGCgene.${seed}.txt"
prefix="${inputDir}/04_simulateInput/PCAWG/gc19_pc.${reg}/${numtag}_${genenum}/PCAWG.${reg}"
/soft/devtools/R/R-4.0.5_installation/bin/Rscript S04_SimDriver.R "${randomfil}" "${prefix}" "${genelog}" "${numtag}" "${genenum}" "${seed}"

echo "/soft/devtools/R/R-4.0.5_installation/bin/Rscript S04_SimDriver.functional.recurrent.R ${randomfil} ${genelog} ${numtag} ${genenum} ${seed}"

################################################################################################################
# step 4 : Run dNdS-Fun
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=excDupSite
#SBATCH --mem=40G
#SBATCH -c 16
#SBATCH -o ./log/dNdS_%A-%a_out.txt
#SBATCH -e ./log/dNdS_%A-%a_error.txt
#SBATCH -p intel-sc3,amd-ep2,amd-ep2-short
#SBATCH --qos=normal
#SBATCH --array=2-20

## user's own commands below
pwd;hostname;date

ID=$1
id=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
seed=${id}
reg=$2
clonal_group=$3
project_group=$4
numtag=$5
genenum=$6

inputDir="/path/"

mkdir -p ${inputDir}/05_simulateOutput_functional/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}/
rda="/path/gc19_pc.${reg}/CADD_score/gc19_pc.${reg}_CADD_score_ALL_0.50.rda"

for drivernum in 5 10 20 50 100 1000
do
mutFil="${inputDir}/04_simulateInput_functional/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}/PCAWGALL.${reg}.${numtag}.${genenum}.${drivernum}.${seed}.tsv"
globalout="${inputDir}/05_simulateOutput_functional/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}/PCAWGALL.${reg}.${numtag}.${genenum}.${drivernum}.${seed}.global.out"
geneout="${inputDir}/05_simulateOutput_functional/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}/PCAWGALL.${reg}.${numtag}.${genenum}.${drivernum}.${seed}.dndscv.out"

/soft/devtools/R/R-4.0.5_installation/bin/Rscript S05_dNdSFun.R ${mutFil} ${rda} ${reg} ${globalout} ${geneout}
done

