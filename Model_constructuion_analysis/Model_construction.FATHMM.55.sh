################################################################################################################
# step 1 : Construct models of FATHMM scores
# this step is to construct the model for different genomic regions
# we extract all possible mutations with functional impact scores for each genomic regions
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=excDupSite
#SBATCH --mem=10G
#SBATCH -o ./log/make_%A-%a_out.txt
#SBATCH -e ./log/make_%A-%a_error.txt
#SBATCH -p amd-ep2-short,amd-ep2,intel-sc3
#SBATCH --qos=normal
#SBATCH --array=1-22

## user's own commands below
# Print the working directory, hostname, and date for logging purposes
pwd; hostname; date

# Set task ID based on the line from the provided ID file corresponding to the SLURM task ID
ID=$1
id=$(head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1)

# Define working directories
runDir="/path/"
wkDir="${runDir}/QuantSet_${neutral}_${quant}"
mkdir -p ${wkDir}

# Define annotation and output directories
gencode_anno="gc19_pc.${4}"
mkdir -p ${wkDir}/${gencode_anno}/FATHMM_score

# Set chromosome variable
chr="chr${id}"
outfile="${wkDir}/${gencode_anno}/FATHMM_score/${gencode_anno}_FATHMM_score_${chr}.excDupSite.txt"
bed="/path/my_resource/${gencode_anno}/${gencode_anno}.excDupSite.bed"

# Initialize output file with headers
rm -f ${outfile}
echo -ne "gene.id\tgene.name\telement.id\tchr\tpos\tref\talt\tRawScore\tPHRED\tRawScore_rankscore\tinterval.start\tinterval.end\tstrand\n" >${outfile}
cat ${bed} | awk -v var=$chr '{if($1==var)print $1":"$2+1"-"$3}' | sort -k 1,1 -T /path/03_Tri_QuantSet_FATHMM/tmp/${4}_${id}/ | uniq | xargs tabix /path/${chr}.hg19.FIS.anno.FATHMM.txt.gz | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3,$4,$5}' | bedtools intersect -a stdin -b ${bed} -wo | awk 'BEGIN{OFS="\t"}{split($10,a,"::");print a[4],a[3],$10,$1,$3,$4,$5,$6,$6,$6,$8+1,$9,$12}' | sort -T /path/tmp/${4}_${id}/ | uniq >>${outfile}

rm ${wkDir}/${gencode_anno}/FATHMM_score/${chr}_o
cat ${wkDir}/${gencode_anno}/FATHMM_score/${gencode_anno}_FATHMM_score_${chr}.excDupSite.txt | cut -f8 | grep -v "RawScore" >>${wkDir}/${gencode_anno}/FATHMM_score/${chr}_o
cat ${wkDir}/${gencode_anno}/FATHMM_score/${gencode_anno}_FATHMM_score_${chr}.excDupSite.txt | cut -f 2,4,5,6,7,8,9,10 | grep --color=auto -v ^"gene.name" > ${wkDir}/${gencode_anno}/FATHMM_score/${gencode_anno}_FATHMM_score_${chr}.impScore.txt

# Print the working directory, hostname, and date again for logging
pwd; hostname; date


################################################################################################################
# step 2 : Get the dichotomy threshold of all mutations using FATHMM scores
# here, we used 5:5 ratio for more functional and less functional mutations 
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=excDupSite
#SBATCH --mem=20G
#SBATCH -o ./log/cal_globaldnds_%A-%a_out.txt
#SBATCH -e ./log/cal_globaldnds_%A-%a_error.txt
#SBATCH -p amd-ep2,amd-ep2-short,intel-sc3
#SBATCH --qos=normal
#SBATCH --array=1-1

## user's own commands below
pwd;hostname;date

module load R/4.0.5
reg=${1} #including cds,ss,prom,5utr,3utr,whole-gene,intergenic,intron
neutral=50
quant=0.50
runDir="/path/"
wkDir="/path/QuantSet_${neutral}_${quant}"
mkdir -p ${wkDir}

gencode_anno="gc19_pc.${reg}"
mkdir -p ${wkDir}/${gencode_anno}
mkdir -p ${wkDir}/${gencode_anno}/FATHMM_score

cat ${runDir}/QuantSet_${neutral}_${quant}/${gencode_anno}/FATHMM_score/chr*_o > ${runDir}/QuantSet_${neutral}_${quant}/${gencode_anno}/FATHMM_score/total_o
threshold_positive=`Rscript ./quantile.R ./QuantSet_50_0.50/gc19_pc.${reg}/FATHMM_score/total_o ${quant}`
threshold_neutral=${threshold_positive}

rm ${wkDir}/${gencode_anno}/threshold.txt
echo -ne "${quant}\t${threshold_positive}\tpositive\n" >>${wkDir}/${gencode_anno}/threshold.txt
echo -ne "${quant}\t${threshold_neutral}\tneutral\n" >>${wkDir}/${gencode_anno}/threshold.txt

################################################################################################################
# step 3 : Construct the dNdS-Fun model for each chromosome of FATHMM scores
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=buildrda
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200G
#SBATCH -o ./log/cal_globaldnds_%A-%a_out.txt
#SBATCH -e ./log/cal_globaldnds_%A-%a_error.txt
#SBATCH -p intel-sc3,amd-ep2,amd-ep2-short
#SBATCH --qos=huge
#SBATCH --array=1-22

## user's own commands below
ID=${1}
id=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
chr=chr${id}
neutral=50
quant=0.50
gencode_anno="gc19_pc.${4}"
wkDir="/path/QuantSet_${neutral}_${quant}"

threshold_fil="${wkDir}/${gencode_anno}/threshold.txt"
threshold_positive=`cat ${threshold_fil} | cut -f 2 | head -n 1`
threshold_neutral=`cat ${threshold_fil} | cut -f 2 | tail -n 1`

elementfile="./QuantSet_${neutral}_0.${neutral}/${gencode_anno}/FATHMM_score/${gencode_anno}_FATHMM_score_${chr}.excDupSite.txt"
genomefile="/path/GRCh37.primary_assembly.genome.fasta"
outfile="${wkDir}/${gencode_anno}/FATHMM_score/${gencode_anno}_FATHMM_score_${chr}.rda"

Rscript /storage/yangjianLab/zhengmengyue/PROJECT/05_SELECTION/01_XIWEI/my_scripts/buildRefElement_FATHMM_score_rda.R ${elementfile} ${genomefile} ${outfile} ${chr} RawScore ${threshold_positive} ${threshold_neutral}
echo "run done"

################################################################################################################
# step 4 : Merge 22 chromosomes together
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=buildrda
#SBATCH --ntasks-per-node=2
#SBATCH --mem=50G
#SBATCH -o ./log/cal_globaldnds_%A-%a_out.txt
#SBATCH -e ./log/cal_globaldnds_%A-%a_error.txt
#SBATCH -p amd-ep2,amd-ep2-short,intel-sc3
#SBATCH --qos=normal
#SBATCH --array=1-1

## user's own commands below
ID=${1}
reg=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
neutral=50
quant=0.50
gencode_anno="gc19_pc.${reg}"
wkDir="/path/QuantSet_${neutral}_${quant}"

threshold_fil="${wkDir}/${gencode_anno}/threshold.txt"
threshold_positive=`cat ${threshold_fil} | cut -f 2 | head -n 1`
threshold_neutral=`cat ${threshold_fil} | cut -f 2 | tail -n 1`

prefix="${wkDir}/${gencode_anno}/FATHMM_score/${gencode_anno}_FATHMM_score_"
outfile="${wkDir}/${gencode_anno}/FATHMM_score/${gencode_anno}_FATHMM_score_ALL_${quant}.rda"
Rscript ./combine_rda.R ${prefix} ${outfile}

