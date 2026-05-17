################################################################################################################
# Step 1 : Generate the trinucleotide pools from the reference genomes
# This pool is prepare to sampling mutations under null
################################################################################################################
# Clear the workspace
rm(list = ls())

# Load necessary library
library("data.table")

# Parse arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
model <- args[1]
output <- args[2]

# Function to swap A/T and C/G bases and reverse the DNA sequence
swap_bases_and_reverse <- function(dna_string) {
  # Temporarily replace A and T with placeholders to prevent conflicts during replacement
  temp_string <- chartr("AT", "XY", dna_string)
  # Swap C and G
  temp_string <- chartr("CG", "GC", temp_string)
  # Replace placeholders with A and T
  final_string <- chartr("XY", "TA", temp_string)
  # Reverse the sequence
  out_string <- rev(strsplit(final_string, "")[[1]])
  return(out_string)
}

# Load the pre-built model, which contains the RefElement list of gene elements
load(model) 
genenum <- length(RefElement)          # Number of genes
partnum <- ceiling(genenum / 100)      # Number of parts to split the genes into (100 per part)

# Loop through each part and process the gene data
for (part in seq_len(partnum)) {
  file.remove(paste0(output, ".part_", part, ".txt")) # Remove existing part file

  # Define gene index range for the current part
  genestart <- ((part - 1) * 100) + 1
  geneend <- min(part * 100, genenum)

  # Loop through genes within the specified range
  for (x1 in seq(genestart, geneend)) {
    res_dat <- data.frame()               # Initialize results data frame for each gene
    genelist <- RefElement[[x1]]          # Extract gene information

    # Extract chromosome and strand information
    chromosome <- gsub("chr", "", genelist$chr)
    genename <- genelist$gene_name
    nr <- nrow(genelist$intervals_element) # Number of intervals within the gene

    # Retrieve sequence elements based on strand orientation
    if (genelist$strand == "+") {
      upseqelement <- genelist$seq_element1up
      seqelement <- genelist$seq_element
      downseqelement <- genelist$seq_element1down
    } else {
      upseqelement <- swap_bases_and_reverse(as.character(genelist$seq_element1up))
      seqelement <- swap_bases_and_reverse(as.character(genelist$seq_element))
      downseqelement <- swap_bases_and_reverse(as.character(genelist$seq_element1down))
    }

    # Generate a list of all positions within each interval
    position_list <- list()
    for (x2 in seq_len(nr)) {
      start <- genelist$intervals_element[x2, 1]
      end <- genelist$intervals_element[x2, 2]
      position_list <- c(position_list, seq(start, end))
    }

    # Generate possible alternate alleles at each position and record tri-nucleotide context
    tnum <- length(genelist$seq_element)
    for (x3 in seq_len(tnum)) {
      nt <- c("A", "C", "G", "T")
      ref <- as.character(seqelement[x3])
      alt <- nt[nt != ref]  # Exclude the reference base

      for (alt_base in alt) {
        triref <- paste0(upseqelement[x3], seqelement[x3], downseqelement[x3])
        trialt <- paste0(upseqelement[x3], alt_base, downseqelement[x3])
        res <- c(chromosome, position_list[x3], ref, alt_base, triref, trialt)
        res_dat <- rbind(res_dat, res)
      }
    }

    # Remove duplicate rows and save results for the current part
    res_dat <- unique(res_dat)
    fwrite(res_dat, paste0(output, ".part_", part, ".txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
  }
}

################################################################################################################
# step 2 : generated the pools of each chromosome for each genomic regions using above scripts
################################################################################################################
rm run_01_Model2Pool.sh
for reg in "cds-exon" "ss" "5utr" "3utr" "prom" "whole-gene"
do
	echo "sbatch sbatch_01_Model2Pool.sh chromosome.txt ${reg}" >> run_01_Model2Pool.sh
done

################################################################################################################
# Step 3 : Add the functional impact scores into the pool
# These functional impact scores doesn't work in the sampling process
# Just to avoid duplication of effort to add functional impact scores to each sampling file after sampling
################################################################################################################
library(data.table);library(R.utils)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
cadd <- args[1]      # Path to CADD data file
input <- args[2]     # Path to input data file
output <- args[3]    # Path to output file

# Load and preprocess Pool data
Pool_dat <- fread(input)
colnames(Pool_dat) <- c("chr", "pos", "ref", "alt", "triref", "trialt")
Pool_dat$key <- paste0("chr", Pool_dat$chr, Pool_dat$pos, Pool_dat$ref, Pool_dat$alt)
Pool_dat2 <- Pool_dat[, .(key, triref, trialt)]

# Load and preprocess CADD reference data
Ref_dat <- fread(cadd)
colnames(Ref_dat) <- c("gene", "chr", "pos", "ref", "alt", "RawScore", "PHRED", "RawScore_rankscore")
Ref_dat$key <- paste0(Ref_dat$chr, Ref_dat$pos, Ref_dat$ref, Ref_dat$alt)

# Convert Ref data to data.table with 'key' for merging
dt1 <- data.table(Pool_dat2, key = "key")
dt2 <- data.table(Ref_dat, key = "key")
dtres <- merge(dt1, dt2, by = "key")

# Select relevant columns for the output and write to file
dtout <- dtres[, .(gene, chr, pos, ref, alt, RawScore, PHRED, RawScore_rankscore, triref, trialt)]
fwrite(dtout, output, sep = "\t", row.names = FALSE, quote = FALSE)

################################################################################################################
# Step 4 : Add the scores for each chromosome for each genomic regions using above scripts
# Merge 22 chromosomes together
################################################################################################################
rm run_02_PoolAddScore.sh
for reg in "cds-exon" "ss" "5utr" "3utr" "prom" "whole-gene"
do
	echo "sbatch sbatch_02_PoolAddScore.sh chromosome.txt ${reg}" >> run_02_PoolAddScore.sh
done

#merge the caddscores file together
rm run_02_MergeScore.sh
for reg in "cds-exon" "ss" "5utr" "3utr" "prom" "whole-gene"
do
  echo "cat 02_CADDPool/gc19_pc.${reg}/gc19_pc.${reg}_CADD_trinucl_chr*.txt > 02_CADDPool/gc19_pc.${reg}/gc19_pc.${reg}_CADD_trinucl_ALL.txt" >> run_02_MergeScore.sh
 done

################################################################################################################
# Step 5 : Do the simulations one time
################################################################################################################
library(data.table)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
realfil <- args[1]  # File path for the real data file
input <- args[2]    # File path for the reference data
output <- args[3]   # File path for the output file
numtag <- as.numeric(args[4])  # Number of samples to draw
seed <- as.numeric(args[5])    # Random seed for reproducibility
set.seed(seed)

# Load and sample the real data
real_dat_0 <- as.data.frame(fread(realfil))
sample_size <- min(numtag, nrow(real_dat_0))  # Ensure sample size is within data bounds
real_dat <- real_dat_0[sample(1:nrow(real_dat_0), sample_size, replace = TRUE), ]

# Set column names for consistency
colnames(real_dat) <- c("donor", "chr", "pos", "ref", "alt", "RawScore", "PHRED", "RawScore_rankscore")
real_dat$key <- paste0(real_dat$chr, real_dat$pos, real_dat$ref, real_dat$alt)
dt1 <- data.table(real_dat, key = "key")

# Load and prepare the reference data
ref_dat <- as.data.frame(fread(input))
colnames(ref_dat) <- c("gene", "chr", "pos", "ref", "alt", "RawScore", "PHRED", "RawScore_rankscore", "triref", "trialt")
ref_dat$key <- paste0(ref_dat$chr, ref_dat$pos, ref_dat$ref, ref_dat$alt)
ref_dat$trikey <- paste0(ref_dat$triref, ">", ref_dat$trialt)
ref_dat_subset <- ref_dat[, .(key, triref, trialt, trikey)]
reftri_counts <- as.data.frame(table(ref_dat_subset$trikey))

# Merge sampled data with reference data by 'key'
dt2 <- data.table(ref_dat_subset, key = "key")
dtres <- merge(dt1, dt2, all.x = TRUE)
dtres_na <- dtres[is.na(dtres$trialt), ]  # Subset rows with missing 'trialt' values

# Count trinucleotide occurrences in merged data
tri_dat <- as.data.frame(table(dtres$trikey))
colnames(tri_dat) <- c("trinucleotide", "count")

# Initialize results data frame
res_dat <- data.frame()

# Sample and prepare data for each trinucleotide type
for (x1 in 1:nrow(tri_dat)) {
  trinucleotide <- tri_dat$trinucleotide[x1]
  count <- as.numeric(tri_dat$count[x1])
  
  # Filter reference data for current trinucleotide and sample rows
  ref_dat_tmp <- ref_dat[ref_dat$trikey == trinucleotide, ]
  sample_indices <- sample(1:nrow(ref_dat_tmp), size = count, replace = TRUE)
  ref_dat_sampled <- ref_dat_tmp[sample_indices, ]
  
  # Assign donor labels and select relevant columns
  ref_dat_sampled$donor <- paste0("donor", seq_len(nrow(ref_dat_sampled)))
  sampled_subset <- ref_dat_sampled[, c("donor", "chr", "pos", "ref", "alt", "RawScore", "PHRED", "RawScore_rankscore", "triref", "trialt")]
  
  # Append sampled data to results
  res_dat <- rbind(res_dat, sampled_subset)
}

# Replace donor labels with original donor list
donor_list <- real_dat$donor
sample_indices_donors <- sample(1:length(donor_list), nrow(res_dat), replace = TRUE)
res_dat$donor <- donor_list[sample_indices_donors]

# Write results to output file
fwrite(res_dat, output, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

################################################################################################################
# Step 6 : Do the simulation multiple times 
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=excDupSite
#SBATCH --mem=25G
#SBATCH -o ./log5/make_%A-%a_out.txt
#SBATCH -e ./log5/make_%A-%a_error.txt
#SBATCH -p intel-sc3,amd-ep2,amd-ep2-short
#SBATCH --qos=normal
#SBATCH --array=1-300

## user's own commands below
pwd;hostname;date

ID=$1
id=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
seed=${id}
reg=$2 #set the reg in "cds-exon" "ss" "5utr" "3utr" "prom" "whole-gene" "intron" "intergenic"
clonal_group=$3
project_group=$4
numtag=$5 ## set the toal number of mutations with 500，1K，10K，100K，1000K

realpath="/path/"
runDir="/path/simulation/"
mkdir -p ${runDir}/03_SimData/Input/PCAWG/gc19_pc.${reg}/${numtag}/
realfil="${realpath}/${clonal_group}/PCAWG_SNVs_allmutations.${project_group}.gc19_pc.${reg}_CADD_score.all.impScore.txt"
input="${runDir}/02_CADDPool/gc19_pc.${reg}/gc19_pc.${reg}_CADD_trinucl_ALL.txt"
output="${runDir}/03_SimData/Input/PCAWG/gc19_pc.${reg}/${numtag}/PCAWG.${reg}.${numtag}.${seed}.txt"

/soft/devtools/R/R-4.0.5_installation/bin/Rscript S03_SimData.R ${realfil} ${input} ${output} ${numtag} ${seed}


################################################################################################################
# Step 7: Run dNdS-Fun for each sampled file
# The dNdS-Fun method integrates global selection assessment and identification of selected genes.
# As a result, this step produces both global selection metrics and selected gene lists under the null model.
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=CADD_t
#SBATCH --ntasks-per-node=3
#SBATCH --mem=40G
#SBATCH -o ./log3/CADD_%A-%a_out.txt
#SBATCH -e ./log3/CADD_%A-%a_error.txt
#SBATCH -p amd-ep2,amd-ep2-short,intel-sc3
#SBATCH --qos=huge
#SBATCH --array=1-300

## user's own commands below
#conda deactivate
pwd;hostname;date

ID=$1
seed=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
reg=$2
neutral=$3
quant=$4
clonal_group=$5
project_group=$6
mut_resource=$7
negbeta=$8
model=$9
numtag=$10

modelDir="/path/"
gencode_anno="gc19_pc.${reg}"
threshold_fil="${modelDir}/QuantSet_${neutral}_${quant}/${gencode_anno}/threshold.txt"
threshold_positive=`cat ${threshold_fil} | cut -f 2 | head -n 1`
threshold_neutral=`cat ${threshold_fil} | cut -f 2 | tail -n 1`

runDir="/path/"
mkdir -p ${runDir}/gc19_pc.${reg}/${numtag}/

##files
R_command="/path/src/cal_CADD_dndsWGS.NEG.R"
mutsFile="${runDir}/Input/${mut_resource}/gc19_pc.${reg}/${numtag}/${mut_resource}.${reg}.${numtag}.${seed}.txt"
refDb_element="${modelDir}/QuantSet_${neutral}_${quant}/${gencode_anno}/CADD_score/${gencode_anno}_CADD_score_ALL_${quant}.rda"

globaldnds_outFile="${runDir}/gc19_pc.${reg}/${numtag}/CADD_globadnds.${mut_resource}.${reg}.${project_group}.${seed}.out"
dndsloc_outFile="${runDir}/gc19_pc.${reg}/${numtag}/CADD_dndsloc.${mut_resource}.${reg}.${project_group}.${seed}.out"
dndscv_outFile="${runDir}/gc19_pc.${reg}/${numtag}/CADD_dndscv.${mut_resource}.${reg}.${project_group}.${seed}.out"
iscv="nocv"
gene_group="all"

##command
/soft/devtools/R/R-4.0.5_installation/bin/Rscript ${R_command} ${mutsFile} ${refDb_element} ${reg} ${neutral} ${quant} ${threshold_positive} ${threshold_neutral} ${gene_group} ${clonal_group} ${project_group} ${mut_resource} ${globaldnds_outFile} ${dndsloc_outFile} ${dndscv_outFile} ${negbeta} ${iscv} ${model}
