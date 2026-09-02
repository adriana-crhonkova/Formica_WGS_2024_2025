##### FILTERING SAMPLES
# 36_vcffilt.sh
# ------------- START OF THE BASH SCRIPT  ---------------------------------------------------------------------------------------------
#!/bin/bash
#SBATCH --job-name=vcffilt
#SBATCH --output=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/pca/vcffilt_%j.out
#SBATCH --error=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/pca/vcffilt_%j.err
#SBATCH --time=01:00:00
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --mail-user=ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# Load environment
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/bcftools.env

echo "Setting Paths"
VCFPATH=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/filt
VCFIN=DP8.AN10.noScaff0003.mac2.thin20kb.vcf.gz

echo "Filtering Samples aq, lu, ru, pol, pra + references"
vcftools --gzvcf $VCFPATH/$VCFIN \
  --keep $VCFPATH/ind_list_D6 \
  --recode --recode-INFO-all \
  --stdout | bgzip > $VCFPATH/DP8.AN10.noScaff0003.mac2.thin20kb.LMUF_Ref_aqlurupolpra.vcf.gz
#  ind_list_D6 contains all LMUF except F. exsecta + reference samples from Finland and Switzerland


echo "Filtering Samples aq, lu, ru, pol + references"
vcftools --gzvcf $VCFPATH/$VCFIN \
  --keep $VCFPATH/ind_list_RefLMUF_aqlupolru.txt \
  --recode --recode-INFO-all \
  --stdout | bgzip > $VCFPATH/DP8.AN10.noScaff0003.mac2.thin20kb.LMUF_Ref_aqlurupol.vcf.gz
#  ind_list_RefLMUF_aqlupolru.txt contains all LMUF and References except F. exsecta and F. pratensis

echo "Filtering Samples aq, lu, ru, pol, pra"
vcftools --gzvcf $VCFPATH/$VCFIN \
  --keep $VCFPATH/ind.list.onlyLMUF.txt \
  --recode --recode-INFO-all \
  --stdout | bgzip > $VCFPATH/DP8.AN10.noScaff0003.mac2.thin20kb.LMUFaqlurupolpra.vcf.gz
# ind.list.onlyLMUF.txt contains all LMUF except F. exsecta
  

echo "Filtering Samples aq, lu, ru, pol"
vcftools --gzvcf $VCFPATH/$VCFIN \
  --keep $VCFPATH/ind.list.LMUF.aqlupolru.txt \
  --recode --recode-INFO-all \
  --stdout | bgzip > $VCFPATH/DP8.AN10.noScaff0003.mac2.thin20kb.LMUFaqlurupol.vcf.gz
# ind.list.LMUF.aqlupolru.txt contains all LMUF except F. exsecta and F. pratensis

echo "DONE" 

# ------------- END OF THE BASH SCRIPT  --------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 4156866 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State
sacct -M biohpc_gen -j 4157220 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State



15_PCA.sh
# ------------- START OF THE BASH SCRIPT  ---------------------------------------------------------------------------------------------

#!/bin/bash
#SBATCH --job-name=PCA
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --output=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/pca/pca_snprelate_%j.out
#SBATCH --error=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/pca/pca_snprelate_%j.err

# Load environment
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/pca.env

# Run the R script
Rscript /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/run_pca_snprelate.R

# ------------- END OF THE BASH SCRIPT  --------------------------------------------------------------------------------------------------


# LMUF and references aquilonia, lugubris, polyctena, rufa, pratensis
# ----------------------------------------------- R script --------------------------------------------------------
# run_pca_snprelate_LMUF_Ref.R
library(SNPRelate)

# Define input and output paths
vcf.fn <- "/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/filt/DP8.AN10.noScaff0003.mac2.thin20kb.LMUF_Ref_aqlurupolpra.vcf.gz"
gds.fn <- "/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/pca/DP8_LMUFRefaqlurupolpra.gds"
out.prefix <- "/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/pca/DP8_LMUFRefaqlurupolpra_PCA"

cat("=== Converting VCF to GDS ===\n")
snpgdsVCF2GDS(vcf.fn, gds.fn, method="copy.num.of.ref")

cat("=== Opening GDS file ===\n")
genofile <- snpgdsOpen(gds.fn)

cat("=== Checking SNP summary ===\n")
print(snpgdsSummary(genofile))

cat("=== Running PCA ===\n")
pca <- snpgdsPCA(genofile, num.thread=1, autosome.only=FALSE)
pc.percent <- pca$varprop * 100

# Save PCA results
pca.df <- data.frame(
  sample = pca$sample.id,
  PC1 = pca$eigenvect[,1],
  PC2 = pca$eigenvect[,2],
  PC3 = pca$eigenvect[,3],
  PC4 = pca$eigenvect[,4],
  stringsAsFactors = FALSE
)

write.csv(pca.df, paste0(out.prefix, "_scores.csv"), row.names = FALSE)
write.csv(round(pc.percent, 2), paste0(out.prefix, "_variance.csv"), row.names = FALSE)

cat("=== PCA completed successfully ===\n")
snpgdsClose(genofile)


# -------------------------------------------------------------- END of R script -----------------------------------------------------------------
sacct -M biohpc_gen -j 4156938 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State

# LMUF and Ref aquilonia, lugubris, polyctena, rufa
# ----------------------------------------------- R script --------------------------------------------------------
# run_pca_snprelate_RefLMUF_aqlurupol.R
library(SNPRelate)

# Define input and output paths
vcf.fn <- "/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/filt/DP8.AN10.noScaff0003.mac2.thin20kb.LMUF_Ref_aqlurupol.vcf.gz"
gds.fn <- "/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/pca/DP8_LMUFRefaqlurupol.gds"
out.prefix <- "/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/pca/DP8_LMUFRefaqlurupol_PCA"

cat("=== Converting VCF to GDS ===\n")
snpgdsVCF2GDS(vcf.fn, gds.fn, method="copy.num.of.ref")

cat("=== Opening GDS file ===\n")
genofile <- snpgdsOpen(gds.fn)

cat("=== Checking SNP summary ===\n")
print(snpgdsSummary(genofile))

cat("=== Running PCA ===\n")
pca <- snpgdsPCA(genofile, num.thread=1, autosome.only=FALSE)
pc.percent <- pca$varprop * 100

# Save PCA results
pca.df <- data.frame(
  sample = pca$sample.id,
  PC1 = pca$eigenvect[,1],
  PC2 = pca$eigenvect[,2],
  PC3 = pca$eigenvect[,3],
  PC4 = pca$eigenvect[,4],
  stringsAsFactors = FALSE
)

write.csv(pca.df, paste0(out.prefix, "_scores.csv"), row.names = FALSE)
write.csv(round(pc.percent, 2), paste0(out.prefix, "_variance.csv"), row.names = FALSE)

cat("=== PCA completed successfully ===\n")
snpgdsClose(genofile)


# -------------------------------------------------------------- END of R script -----------------------------------------------------------------
sacct -M biohpc_gen -j 4157223 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State





# LMUF aquilonia, lugubris, polyctena, rufa, pratensis
# ----------------------------------------------- R script --------------------------------------------------------
# run_pca_snprelate_LMUF_aqlurupolpra.R
library(SNPRelate)

# Define input and output paths
vcf.fn <- "/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/filt/DP8.AN10.noScaff0003.mac2.thin20kb.LMUFaqlurupolpra.vcf.gz"
gds.fn <- "/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/pca/DP8_LMUFaqlurupolpra.gds"
out.prefix <- "/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/pca/DP8_LMUFaqlurupolpra_PCA"

cat("=== Converting VCF to GDS ===\n")
snpgdsVCF2GDS(vcf.fn, gds.fn, method="copy.num.of.ref")

cat("=== Opening GDS file ===\n")
genofile <- snpgdsOpen(gds.fn)

cat("=== Checking SNP summary ===\n")
print(snpgdsSummary(genofile))

cat("=== Running PCA ===\n")
pca <- snpgdsPCA(genofile, num.thread=1, autosome.only=FALSE)
pc.percent <- pca$varprop * 100

# Save PCA results
pca.df <- data.frame(
  sample = pca$sample.id,
  PC1 = pca$eigenvect[,1],
  PC2 = pca$eigenvect[,2],
  PC3 = pca$eigenvect[,3],
  PC4 = pca$eigenvect[,4],
  stringsAsFactors = FALSE
)

write.csv(pca.df, paste0(out.prefix, "_scores.csv"), row.names = FALSE)
write.csv(round(pc.percent, 2), paste0(out.prefix, "_variance.csv"), row.names = FALSE)

cat("=== PCA completed successfully ===\n")
snpgdsClose(genofile)


# -------------------------------------------------------------- END of R script -----------------------------------------------------------------
sacct -M biohpc_gen -j 4156942 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State

# LMUF aquilonia, lugubris, polyctena, rufa
# ----------------------------------------------- R script --------------------------------------------------------
# run_pca_snprelate_LMUF_aqlurupol.R
library(SNPRelate)

# Define input and output paths
vcf.fn <- "/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/filt/DP8.AN10.noScaff0003.mac2.thin20kb.LMUFaqlurupol.vcf.gz"
gds.fn <- "/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/pca/DP8_LMUFaqlurupol.gds"
out.prefix <- "/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/pca/DP8_LMUFaqlurupol_PCA"

cat("=== Converting VCF to GDS ===\n")
snpgdsVCF2GDS(vcf.fn, gds.fn, method="copy.num.of.ref")

cat("=== Opening GDS file ===\n")
genofile <- snpgdsOpen(gds.fn)

cat("=== Checking SNP summary ===\n")
print(snpgdsSummary(genofile))

cat("=== Running PCA ===\n")
pca <- snpgdsPCA(genofile, num.thread=1, autosome.only=FALSE)
pc.percent <- pca$varprop * 100

# Save PCA results
pca.df <- data.frame(
  sample = pca$sample.id,
  PC1 = pca$eigenvect[,1],
  PC2 = pca$eigenvect[,2],
  PC3 = pca$eigenvect[,3],
  PC4 = pca$eigenvect[,4],
  stringsAsFactors = FALSE
)

write.csv(pca.df, paste0(out.prefix, "_scores.csv"), row.names = FALSE)
write.csv(round(pc.percent, 2), paste0(out.prefix, "_variance.csv"), row.names = FALSE)

cat("=== PCA completed successfully ===\n")
snpgdsClose(genofile)


# -------------------------------------------------------------- END of R script -----------------------------------------------------------------
sacct -M biohpc_gen -j 4156943 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State


scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/pca .
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/pca/DP8_LMUFRefaqlurupol_PCA_*.csv .
