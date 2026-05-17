################################################################################################################
# step 1 : GO and KEGG enrichment analysis for Positively selected genes
################################################################################################################
rm(list=ls())
setwd("/Users/lizi/Documents/01_Westlake/01-Yang\'s-lab/18_CancerSomaticGermline/GeneComparison/Germline/")
args <- commandArgs(trailingOnly = TRUE)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(biomaRt)
library(org.Hs.eg.db)
library(data.table)
library(ggplot2)
library(cowplot)
library(stringr)
library(GOSemSim)
library(rrvgo)
library(rutils)

############ GermlineLocal
############ 
germline_local_list <- as.data.frame(fread("./02_out/Positive_selected_genes.txt",header=F))
length(germline_local_list$V1)
germline_local_list$external_gene_name <- germline_local_list$V1
germline_local_list <- as.data.frame(germline_local_list)

ensembl <- useEnsembl(biomart = "ensembl",
                      dataset = "hsapiens_gene_ensembl",
                      GRCh = "37")
mart <- useMart("ensembl","hsapiens_gene_ensembl")
gene <- germline_local_list$external_gene_name
gene_entrez <- bitr(gene,fromType="SYMBOL",
                    toType=c("ENSEMBL", "ENTREZID"),
                    OrgDb="org.Hs.eg.db")
write.csv(gene_entrez,paste0("./02_out/Positive_selected_genes.biomaRt.csv"),row.names=FALSE)
############

data(geneList, package="DOSE") #富集分析的背景基因集
gene <- names(geneList)

############ GermlineLocal
############ 
gene_entrez <- as.data.frame(fread(paste0("./02_out/Positive_selected_genes.biomaRt.csv")))

go_ALL <- enrichGO(gene = unique(gene_entrez$ENTREZID),
                   universe = names(geneList), #背景基因集
                   OrgDb = org.Hs.eg.db,
                   ont = "ALL", #也可以是 CC  BP  MF中的一种
                   pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                   pvalueCutoff = 1, #P值会过滤掉很多，可以全部输出
                   qvalueCutoff = 0.05,
                   readable = FALSE) #Gene ID 转成gene Symbol ，易读
go_sig <- summary(go_ALL)
fwrite(go_sig,paste0("02_out/Positive_selected_genes.GO-enrich.csv"),sep=",")

kk <- enrichKEGG(gene=unique(gene_entrez$ENTREZID),
                 organism='hsa',
                 keyType = "ncbi-geneid",
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "none")

kk_sig <- as.data.frame(kk)
fwrite(as.data.frame(kk_sig),paste0("02_out/Positive_selected_genes.KEGG-enrich.csv"),sep=",")
############



