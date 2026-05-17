################################################################################################################
# step 1 : Calculate the number of overlaping genes in observed data
################################################################################################################
# 加载必要的库
library(data.table)
library(stringr)

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 初始化结果数据框
res_dat <- data.frame()
sig_dat <- data.frame()

# 读取基因分组列表和显著选择基因列表
gene_group_dat <- fread("GeneSet/Gene.list", header = FALSE)
gene_group_list <- unique(gene_group_dat$V1)
reg_list <- c("cds", "5utr", "3utr", "whole-gene", "ss", "prom")

# 读取显著选择基因数据
siggene_dat <- fread("Negatively_selected genes.txt")

# 遍历每个区域类型
for (reg in reg_list) {
  # 获取特定区域的显著选择基因列表
  siggene_dat_reg <- subset(siggene_dat, REGION == reg)
  siggene_list <- unique(siggene_dat_reg$GENE)
  line_num <- length(siggene_list)
  
  # 将区域和显著基因数量存入 sig_dat
  sig_res <- data.frame(region = reg, gene_count = line_num)
  sig_dat <- rbind(sig_dat, sig_res)
  
  # 遍历每个基因组
  for (gene_group in gene_group_list) {
    # 读取基因组数据
    Gene_dat <- fread(paste0("GeneSet/gencode.v19.", gene_group, ".bed"))
    CGC_gene <- unique(Gene_dat$V7)
    
    # 计算交集基因数目
    overlap_gene <- intersect(siggene_list, CGC_gene)
    overlap_num <- length(overlap_gene)
    
    # 存储结果
    res <- data.frame(
      gene_group = gene_group,
      Region = reg,
      samp_num = line_num,
      overlap_num = overlap_num,
      seed = "observed"
    )
    res_dat <- rbind(res_dat, res)
  }
}

# 输出结果到文件
fwrite(res_dat, "S01_OverlapGene.RealDat.result.txt", sep = "\t", quote = FALSE, row.names = FALSE)

################################################################################################################
# step 2 : calculate the p-values for one gene set
################################################################################################################
# 加载必要的库
library(data.table)
library(stringr)

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)
line_num <- as.numeric(args[1])
reg <- args[2]
seed_num <- as.numeric(args[3])
output <- args[4]
gene_group <- args[5]
output2 <- args[6]

# 读取实际数据集并过滤
OverlapGene_dat <- fread("10_GeneOverlap/01_RESOURCE/S01_OverlapGene.RealDat.result.txt")
OverlapGene_dat_filtered <- subset(OverlapGene_dat, Region == reg & gene_group == gene_group, 
                                   select = c("gene_group", "Region", "samp_num", "window", "overlap_num", "seed"))

# 进行采样
res_dat <- data.frame()
for (seed in 1:seed_num) {
    # 设置随机数种子并读取基因池
    pool_dat <- fread(paste0("10_GeneOverlap/01_RESOURCE/gc19_pc.", reg, ".probe.list"), header = FALSE)
    sample_dat <- pool_dat[sample(1:nrow(pool_dat), line_num), ]
    colnames(sample_dat) <- "name"
    sample_dat$gene <- str_split_fixed(sample_dat$name, "::", 4)[, 3]
    
    # 读取基因组数据并计算重叠基因数目
    Gene_dat <- fread(paste0("GeneSet/gencode.v19.", gene_group, ".bed"))
    CGC_gene <- unique(Gene_dat$V7)
    overlap_num <- length(intersect(sample_dat$gene, CGC_gene))
    
    # 将结果存储到 res_dat
    res <- data.frame(gene_group = gene_group, Region = reg, samp_num = line_num, 
                      overlap_num = overlap_num, seed = seed)
    res_dat <- rbind(res_dat, res)
}

# 保存采样结果到文件
fwrite(res_dat, output, sep = "\t", quote = FALSE, row.names = FALSE)

# 计算统计值和 p 值
res_dat$overlap_num <- as.numeric(res_dat$overlap_num)
observed_val <- as.numeric(OverlapGene_dat_filtered$overlap_num)
samp_median <- median(res_dat$overlap_num)
samp_mean <- mean(res_dat$overlap_num)
enrich_median <- observed_val / samp_median
enrich_mean <- observed_val / samp_mean
sd_val <- sd(res_dat$overlap_num)
zscore <- (observed_val - samp_mean) / sd_val

# 计算 p 值
res_dat_combined <- rbind(res_dat, OverlapGene_dat_filtered)
res_dat_combined <- res_dat_combined[order(-res_dat_combined$overlap_num), ]
rank <- which(res_dat_combined$seed == "observed")
pval <- rank / seed_num

# 将统计结果存储到 p_res
p_res <- data.frame(
    gene_group = gene_group, Region = reg, sig = OverlapGene_dat_filtered$samp_num, 
    overlap = observed_val, window = OverlapGene_dat_filtered$window, rank = rank, 
    pval = pval, samp_median = samp_median, samp_mean = samp_mean, 
    enrich_median = enrich_median, enrich_mean = enrich_mean, sd = sd_val, zscore = zscore
)
# 保存统计结果到文件
fwrite(p_res, output2, sep = "\t", quote = FALSE, row.names = FALSE)

################################################################################################################
# step 3 : run the sampling test for different gene sets
################################################################################################################
 #!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=OverlapGene
#SBATCH --ntasks-per-node=2
#SBATCH --mem=30G
#SBATCH -o ./log/OverlapGene_%A-%a_out.txt
#SBATCH -e ./log/OverlapGene_%A-%a_error.txt
#SBATCH -p intel-sc3,amd-ep2,amd-ep2-short
#SBATCH --qos=huge
#SBATCH --array=1-1

# 打印工作目录、主机名和日期
pwd;hostname;date

# 接收命令行参数
line_num=${1}
reg=${2}
seed_num=${3}
mut_resource=${4}
gene_group=${5}

# 工作目录和文件路径
wkDir="/path/"
gencode_anno="gc19_pc.${reg}"
output_dir="${wkDir}/10_GeneOverlap/03_OUT/${mut_resource}"
mkdir -p "${output_dir}"

output="${output_dir}/OverlapGene.${reg}.${gene_group}.${seed_num}.result.txt"
output2="${output_dir}/OverlapGene.${reg}.${gene_group}.${seed_num}.p.txt"

# 运行 R 脚本
Rscript S01_OverlapGene.gene.R "${line_num}" "${reg}" "${seed_num}" "${output}" "${gene_group}" "${output2}"
