################################################################################################################
# Step 1 : Format the samping mutations from dNdS-Fun into DIGDriver and Run DIGDriver
# Please find some other related files of DIGDriver in the Benchmarked Part
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=DIGDriver
#SBATCH --mem=10G
#SBATCH -c 4
#SBATCH -o ./log/Driverpower_%A-%a_out.txt
#SBATCH -e ./log/Driverpower_%A-%a_error.txt
#SBATCH -p amd-ep2,amd-ep2-short,intel-sc3
#SBATCH --qos=normal
#SBATCH --array=1-3

## user's own commands below
ID=${1}
seed=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
reg=${2}
numtag=${3}
genenum=${4}
drivernum=${5}

conda activate DIGDriver
wkDir="/path/"
inputDir="/path/"
DigPreprocess="/path/DigPreprocess.py"
HG19="/path/hg19.fasta"
mutFil="${inputDir}/04_simulateInput_functional/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}/PCAWGALL.${reg}.${numtag}.${genenum}.${drivernum}.${seed}.tsv"
mkdir -p ${inputDir}/04_simulateInput_DIGDriver/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}/
DIGDriverFil="${inputDir}/04_simulateInput_DIGDriver/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}/PCAWGALL.${reg}.${numtag}.${genenum}.${drivernum}.${seed}"

cat ${mutFil}  | awk '{printf $2"\t""%.f""\t""%.f""\t"$1"\t"$4"\t"$5"\n", $3,$3}' | sed 's/^chr//g' > ${DIGDriverFil}".tsv"
cat ${DIGDriverFil}".tsv" | sort -V -k1,1 -k2,2 | gzip -c > ${DIGDriverFil}".tmp.txt.gz"
rm ${DIGDriverFil}".annot.txt"
touch ${DIGDriverFil}".annot.txt"
${DigPreprocess} annotMutationFile ${DIGDriverFil}".tmp.txt.gz" ${HG19} 
                ${DIGDriverFil}".annot.txt" --n-procs 5
gzip -f ${DIGDriverFil}".annot.txt"
rm ${DIGDriverFil}".tmp.txt.gz"

DigPretrain="/path/DigPretrain.py"
DigDriver="/path/DigDriver.py"
Pretrained_model="/path/Pancan_SNV_MNV_INDEL.Pretrained.h5"
rawDir="/path/"
hg19="/path/hg19.fasta"
REGION="/path/annotions/grch37.PCAWG_noncoding.bed"
NAME="PCAWG_all_elts"
SAMPLE="PCAWGALL.${reg}.${numtag}.${genenum}.${drivernum}.${seed}"
outdir="${inputDir}/05_simulateOutput_DIGDriver/PCAWGALL/gc19_pc.${reg}/${numtag}_${genenum}/"
mkdir -p ${outdir}

${DigDriver} elementDriver ${DIGDriverFil}".annot.txt.gz" \
    ${Pretrained_model} $NAME \
    --f-bed $REGION --outpfx ${SAMPLE} --outdir ${outdir}/

################################################################################################################
# Step 1 : Use the hpc to run the scripts multiple times
################################################################################################################
rm run_05_DIGDriver.sh
for Ref_num in 1000 10000 100000 1000000 10000000 100000000 5000 50000 500000 5000000
   do for Gene_num in 5 10 20 50 100
      do for Driver_num in 5 10 20 50 100 1000
      do 
        echo "sbatch sinfo_01_DIGDriver.sh seed.txt cds-exon ${Ref_num} ${Gene_num} ${Driver_num}" >> run_05_DIGDriver.sh
      done
   done
done
