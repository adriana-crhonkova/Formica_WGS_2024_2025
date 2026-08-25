# 35_PCA_mtDNA.sh
------------- START OF THE BASH SCRIPT  ---------------------------------------------------------------------------------------------

#!/bin/bash
#SBATCH --job-name=PCA_mtDNA
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --output=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/pca_snprelate_%j.out
#SBATCH --error=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/pca_snprelate_%j.err

# Load environment
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/pca.env

# Run the R script
Rscript /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/run_pca_snprelate_mtDNA.R

------------- END OF THE BASH SCRIPT  --------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 4110572 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,NodeList


# run_pca_snprelate_mtDNA.R
----------------------------------------------- R script --------------------------------------------------------
library(SNPRelate)

# Define input and output paths
vcf.fn <- "/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz"
gds.fn <- "/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/pca/LMUF_mtDNA.gds"
out.prefix <- "/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/pca/LMUF_mtDNA"

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


-------------------------------------------------------------- END of R script -----------------------------------------------------------------
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/pca .
