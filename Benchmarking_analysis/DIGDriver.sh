################################################################################################################
# step 1 : Format mutations as a zero-indexed bed file with the following six columns
################################################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=DIGDriver
#SBATCH --ntasks-per-node=2
###SBATCH --mem=20G
#SBATCH -o ./log/DIGDriver_%A-%a_out.txt
#SBATCH -e ./log/DIGDriver_%A-%a_error.txt
#SBATCH -p amd-ep2,amd-ep2-short
#SBATCH --qos=huge
#SBATCH --array=1-25

## user's own commands below
ID=${1}
tissue=`head -n ${SLURM_ARRAY_TASK_ID} $ID | tail -n1 | awk '{print $0}'`
tissue="Uterus-AdenoCA"
rawDir="/path/PCAWG/singleproject_tissue/"
wkDir="/path/DIGDriver/PCAWG/"
DigPreprocess="/path/SOFTWARE.bak/DIGDriver/scripts/DigPreprocess.py"
HG19="/path/SOFTWARE.bak/DIGDriver/Test/dig_data_files/hg19.fasta"

cat ${rawDir}/PCAWG_SNVs_allmutations.${tissue}.txt | sed 's/^chr//g' | awk '{printf $1"\t""%.f""\t""%.f""\t"$4"\t"$5"\t"$6"\n", $2,$3}' > ${wkDir}/singleproject_tissue/PCAWG_SNVs_allmutations.${tissue}.DIGDriver.tsv
cat ${wkDir}/singleproject_tissue/PCAWG_SNVs_allmutations.${tissue}.DIGDriver.tsv | sort -V -k1,1 -k2,2 | gzip -c > ${wkDir}/singleproject_tissue/PCAWG_SNVs_allmutations.${tissue}.DIGDriver.tmp.txt.gz

################################################################################################################
# step 2 : Use Dig to annotate the mutation bed file created in step 1
################################################################################################################
rm ${wkDir}/singleproject_tissue/PCAWG_SNVs_allmutations.${tissue}.DIGDriver.annot.txt
${DigPreprocess} annotMutationFile \
        ${wkDir}/singleproject_tissue/PCAWG_SNVs_allmutations.${tissue}.DIGDriver.tmp.txt.gz \
        ${HG19} \
        ${wkDir}/singleproject_tissue/PCAWG_SNVs_allmutations.${tissue}.DIGDriver.annot.txt \
        --n-procs 5
gzip ${wkDir}/singleproject_tissue/PCAWG_SNVs_allmutations.${tissue}.DIGDriver.annot.txt
rm ${wkDir}/singleproject_tissue/PCAWG_SNVs_allmutations.${tissue}.DIGDriver.tmp.txt.gz


rawDir="/path/PCAWG/allproject/"
wkDir="/path/DIGDriver/PCAWG/singleproject_tissue"
tissue="all"
echo -e "#chr\tstart\tend\tref\talt\tid" > ${wkDir}/PCAWG_SNVs_allmutations.${tissue}.gc19_pc.${reg}_CADD_score.DIGDriver.tsv
cat ${rawDir}/PCAWG_SNVs_allmutations.${tissue}.gc19_pc.${reg}_CADD_score.all.impScore.txt | grep -v ^"V1"| awk '{printf $2"\t""%.f""\t""%.f""\t"$4"\t"$5"\t"$1"\n", $3-1,$3}' >> ${wkDir}/PCAWG_SNVs_allmutations.${tissue}.gc19_pc.${reg}_CADD_score.DIGDriver.tsv
gzip ${wkDir}/singleproject_tissue/PCAWGALL_SNVs_allmutations.${tissue}.DIGDriver.tsv

################################################################################################################
# step 3 : Run DIGDriver to find the driver elements
################################################################################################################
${DigDriver} elementDriver \
    ${wkDir}/singleproject_tissue/PCAWG_SNVs_allmutations.${tissue}.DIGDriver.tsv.gz \
    /path/SOFTWARE.bak/DIGDriver/Test/mutation_maps/${Pretrained_model} \
    $NAME \
    --f-bed $REGION \
    --outpfx ${tissue}.$NAME \
    --outdir output/${tissue}/${NAME}

################################################################################################################
# step 4 : The matched DIGDriver models for data
################################################################################################################
# Pretrain Model for PCAWG data
Billiary Biliary-AdenoCA_SNV_MNV_INDEL.Pretrained.h5
Bladder Bladder-TCC_SNV_MNV_INDEL.Pretrained.h5
ColoRect-AdenoCA ColoRect-AdenoCA_SNV_MNV_INDEL_msi_low.Pretrained.h5
Eso-AdenoCA Eso-AdenoCa_SNV_MNV_INDEL.Pretrained.h5
Head-SCC Head-SCC_SNV_MNV_INDEL.Pretrained.h5
Liver Liver-HCC_SNV_MNV_INDEL.Pretrained.h5
Lung-AdenoCA Lung-AdenoCA_SNV_MNV_INDEL.Pretrained.h5
Lung-SCC Lung-SCC_SNV_MNV_INDEL.Pretrained.h5
Lymph-BNHL Lymph-BNHL_SNV_MNV_INDEL.Pretrained.h5
Lymph-CLL Lymph-CLL_SNV_MNV_INDEL.Pretrained.h5
meta_Bone Bone-Osteosarc_SNV_MNV_INDEL.Pretrained.h5
meta_Breast Breast_tumors_SNV_MNV_INDEL.Pretrained.h5
meta_Cervix Female_reproductive_system_tumors_SNV_MNV_INDEL_msi_low.Pretrained.h5
meta_CNS CNS_tumors_SNV_MNV_INDEL.Pretrained.h5
meta_Kidney Kidney_tumors_SNV_MNV_INDEL.Pretrained.h5
meta_Myeloid Hematopoietic_tumors_SNV_MNV_INDEL.Pretrained.h5
meta_SoftTissue Pancan_SNV_MNV_INDEL.Pretrained.h5
Ovary Ovary-AdenoCA_SNV_MNV_INDEL.Pretrained.h5
Panc-AdenoCA Panc-AdenoCA_SNV_MNV_INDEL.Pretrained.h5
Panc-Endocrine Panc-Endocrine_SNV_MNV_INDEL.Pretrained.h5
Prost-AdenoCA Prost-AdenoCA_SNV_MNV_INDEL.Pretrained.h5
Skin-Melanoma Skin-Melanoma_SNV_MNV_INDEL.Pretrained.h5
Stomach-AdenoCA Stomach-AdenoCA_SNV_MNV_INDEL_msi_low.Pretrained.h5
Thy-AdenoCA Pancan_SNV_MNV_INDEL.Pretrained.h5
Uterus-AdenoCA Uterus-AdenoCA_SNV_MNV_INDEL_msi_low.Pretrained.h5

# Pretrain Model for Genomics England 100kGP data
ADULT_GLIOMA: Glioma_tumors_SNV_MNV_INDEL.Pretrained.h5
BLADDER: Bladder-TCC_SNV_MNV_INDEL.Pretrained.h5
BREAST: Breast-AdenoCa_SNV_MNV_INDEL.Pretrained.h5
CHILDHOOD:  Pancan_SNV_MNV_INDEL.Pretrained.h5
COLORECTAL: ColoRect-AdenoCA_SNV_MNV_INDEL_msi_low.Pretrained.h5
ENDOMETRIAL_CARCINOMA: Female_reproductive_system_tumors_SNV_MNV_INDEL_msi_low.Pretrained.h5
HAEMONC: Hematopoietic_tumors_SNV_MNV_INDEL.Pretrained.h5
HEPATOPANCREATOBILIARY: Panc-Endocrine_SNV_MNV_INDEL.Pretrained.h5
LUNG: Lung-AdenoCA_SNV_MNV_INDEL.Pretrained.h5
MALIGNANT_MELANOMA: Skin-Melanoma_SNV_MNV_INDEL.Pretrained.h5
ORAL_OROPHARYNGEAL: Digestive_tract_tumors_SNV_MNV_INDEL_msi_low.Pretrained.h5
OVARIAN: Female_reproductive_system_tumors_SNV_MNV_INDEL_msi_low.Pretrained.h5
PROSTATE: Prost-AdenoCA_SNV_MNV_INDEL.Pretrained.h5
RENAL: Pancan_SNV_MNV_INDEL.Pretrained.h5
SARCOMA: Sarcoma_tumors_SNV_MNV_INDEL.Pretrained.h5
TESTICULAR_GERM_CELL_TUMOURS: Digestive_tract_tumors_SNV_MNV_INDEL_msi_low.Pretrained.h5
UPPER_GASTROINTESTINAL: Digestive_tract_tumors_SNV_MNV_INDEL_msi_low.Pretrained.h5








