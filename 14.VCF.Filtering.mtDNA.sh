#### Scipt Edited from https://zenodo.org/records/7941711?preview_file=Raw_reads_to_SNPs_scripts.zip

# 27_vcf_filtering_mtDNA.sh
------------- START OF THE SCRIPT ---------------------------------------------------------------------------------------------------------------------------------------------

#!/bin/bash -l
#SBATCH -J filt_mtDNA
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/vcf_filt.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/vcf_filt.err
#SBATCH -t 12:00:00
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=6G
#SBATCH --mail-user=ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# set the script to stop immediately if any command fails or any part of a pipeline fails, preventing silent errors and corrupted results
set -eo pipefail

# Load environment
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/bcftools.env


#####
##### 0. Pre-filtering -----------
#####
echo "Pre-Filtering"

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA
bcftools norm -f /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/Reference_Genome/Formica_hybrid_v1_wFhyb_Sapis.fa -m -both -O z -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/filtered_samples_mtDNA_2_norm.vcf.gz /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/filtered_samples_mtDNA_2_raw.vcf.gz

bcftools filter --threads 1 -Oz -s+ --SnpGap 2 filtered_samples_mtDNA_2_norm.vcf.gz > filtered_samples_mtDNA_2_norm.SnpGap_2.vcf.gz && \

bcftools filter --threads 1 -Oz -e 'TYPE!="snp"' -s NonSnp -m+ filtered_samples_mtDNA_2_norm.SnpGap_2.vcf.gz > filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.vcf.gz && \

bcftools filter --threads 1 -Oz -s Balance -m+ -i 'RPL>=1 && RPR>=1 && SAF>=1 && SAR>=1' filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.vcf.gz > filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.vcf.gz && \

bcftools view --threads 1 -O z -f PASS filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.vcf.gz > filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.vcf.gz && \

bcftools view --threads 1 filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.vcf.gz | vcfallelicprimitives --keep-info --keep-geno -t decomposed | sed '/^##/! s/|/\//g' | sed 's/\.:\.:\.:\.:\.:\.:\.:\./\.\/\.:\.:\.:\.:\.:\.:\.:\./g' | bcftools sort --temp-dir $TMPDIR --max-mem 4G -O z > filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz && \


bcftools index -t filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz
echo "bcftools index -n filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz"
bcftools index -n filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz

#####
##### 1. SNP QUAL >= 30, biallelic -------------------------------------------------------------------
#####
echo "SNP QUAL >= 30, biallelic"

bcftools filter --threads 1 --include 'QUAL >= 30 && TYPE="snp"' -Oz filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz > filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz
echo "gunzip -c filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz | grep -vc '#'"
gunzip -c filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz | grep -vc '#'

bcftools view --threads 1 --min-alleles 2 --max-alleles 2 filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz -Oz > filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz
echo "gunzip -c filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz | grep -vc '#'"
gunzip -c filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz | grep -vc '#'
bcftools index -t filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz


#####
##### 2. Correcting the header --------------------------------------------------
#####
echo "Correcting the header"

# Extract header from VCF
bcftools view -h filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz > header.vcf

# Fix fields
perl -npe 's/<ID=AO,Number=A/<ID=AO,Number=\./' header.vcf | perl -npe 's/<ID=AD,Number=R/<ID=AD,Number=\./' | perl -npe 's/<ID=QA,Number=A/<ID=QA,Number=\./' | perl -npe 's/<ID=GL,Number=G/<ID=GL,Number=\./' > header_AO_AD_QA_GL.vcf

# Replace corrected header
bcftools reheader -h header_AO_AD_QA_GL.vcf -o filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz

bcftools index -t filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz
bcftools index -n filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz

echo "DONE"

------------- END OF THE SCRIPT ---------------------------------------------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 4032903 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode

       JobID    Elapsed     MaxRSS  AllocCPUS      State ExitCode
------------ ---------- ---------- ---------- ---------- --------
4032903        00:00:11                     8  COMPLETED      0:0
4032903.bat+   00:00:11                     8  COMPLETED      0:0 

# Ended up with 722 SNPs

------------------- START OF INTERACTIVE JOB -----------------------------------------------------------------------------------------------------------------------------------
salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 00:30:00 --mem=2G
mamba activate bcftools.env

awk '{print $1}' filt_samples_mtDNA_2_qc_stats.txt > samples.list.filt

# Run the loop using samtools faidx to isolate ONLY the mtDNA region
while IFS= read -r s; do
    samtools faidx /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/Reference_Genome/Formica_hybrid_v1_wFhyb_Sapis.fa mtDNA | \
    bcftools consensus -s "$s" -I filtered_samples_mtDNA_2_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz > "${s}_mtDNA.fa"
    
    # Clean up the fasta header so it is just the sample name
    sed -i "s/>.*/>${s}/" "${s}_mtDNA.fa"
done < samples.list.filt

mkdir filt_samples_mtDNA_fa_2
mv *fa filt_samples_mtDNA_fa_2


------------------- END OF INTERACTIVE JOB --------------------------------------------------------------------------------------------------------------------------
28_filt.alignment_mtDNA_2.sh
-------------------- START OF BASH JOB -----------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash 
#SBATCH -J filt.alignment_mtDNA_2
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/filt.alignment_mtDNA_2.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/filt.alignment_mtDNA_2.err
#SBATCH -t 00:30:00
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --mail-user ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# set the script to stop immediately if any command fails or any part of a pipeline fails, preventing silent errors and corrupted results
set -eo pipefail

# Load environment
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/mafft.env

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/filt_samples_mtDNA_fa_2

echo "Combining FASTA files"
# Concatenate all sample FASTAs into one file
cat *.fa > filt_mtDNAsamples.fa

echo "Alignment"
# alignment
## --thread -1 tells MAFT to use all allocated CPUs
## --anysymbol accept any valid text character which could occur because of missing data
mafft --anysymbol --thread -1 filt_mtDNAsamples.fa > filt_mtDNAsamples_aligned_2.fa

echo "DONE" 

-------------------- END OF BASH JOB ----------------------------------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 4032913 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode


salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 00:15:00 --mem=2G
mamba activate mafft.env

perl Fasta2NEXUS.pl /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/filt_samples_mtDNA_fa_2/filt_mtDNAsamples_aligned_2.fa /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/filt_samples_mtDNA_fa_2/filt_mtDNAsamples_2_aligned.nex


# use in popart to build the network
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/filt_samples_mtDNA_fa_2/filt_mtDNAsamples_2_aligned.nex .



