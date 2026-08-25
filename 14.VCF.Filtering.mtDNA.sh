#### Scipt Edited from https://zenodo.org/records/7941711?preview_file=Raw_reads_to_SNPs_scripts.zip

# 30_vcf_filtering_mtDNA.sh
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
bcftools norm -f /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/Reference_Genome/Formica_hybrid_v1_wFhyb_Sapis.fa -m -both -O z -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/filtered_samples_mtDNA_norm.vcf.gz /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/filtered_samples_mtDNA_raw.vcf.gz

bcftools filter --threads 1 -Oz -s+ --SnpGap 2 filtered_samples_mtDNA_norm.vcf.gz > filtered_samples_mtDNA_norm.SnpGap_2.vcf.gz && \

bcftools filter --threads 1 -Oz -e 'TYPE!="snp"' -s NonSnp -m+ filtered_samples_mtDNA_norm.SnpGap_2.vcf.gz > filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.vcf.gz && \

bcftools filter --threads 1 -Oz -s Balance -m+ -i 'RPL>=1 && RPR>=1 && SAF>=1 && SAR>=1' filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.vcf.gz > filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.vcf.gz && \

bcftools view --threads 1 -O z -f PASS filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.vcf.gz > filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.vcf.gz && \

bcftools view --threads 1 filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.vcf.gz | vcfallelicprimitives --keep-info --keep-geno -t decomposed | sed '/^##/! s/|/\//g' | sed 's/\.:\.:\.:\.:\.:\.:\.:\./\.\/\.:\.:\.:\.:\.:\.:\.:\./g' | bcftools sort --temp-dir $TMPDIR --max-mem 4G -O z > filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz && \


bcftools index -t filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz
echo "bcftools index -n filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz"
bcftools index -n filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz

#####
##### 1. SNP QUAL >= 30, biallelic -------------------------------------------------------------------
#####
echo "SNP QUAL >= 30, biallelic"

bcftools filter --threads 1 --include 'QUAL >= 30 && TYPE="snp"' -Oz filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz > filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz
echo "gunzip -c filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz | grep -vc '#'"
gunzip -c filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz | grep -vc '#'

bcftools view --threads 1 --min-alleles 2 --max-alleles 2 filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz -Oz > filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz
echo "gunzip -c filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz | grep -vc '#'"
gunzip -c filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz | grep -vc '#'
bcftools index -t filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz


#####
##### 2. Correcting the header --------------------------------------------------
#####
echo "Correcting the header"

# Extract header from VCF
bcftools view -h filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz > header.vcf

# Fix fields
perl -npe 's/<ID=AO,Number=A/<ID=AO,Number=\./' header.vcf | perl -npe 's/<ID=AD,Number=R/<ID=AD,Number=\./' | perl -npe 's/<ID=QA,Number=A/<ID=QA,Number=\./' | perl -npe 's/<ID=GL,Number=G/<ID=GL,Number=\./' > header_AO_AD_QA_GL.vcf

# Replace corrected header
bcftools reheader -h header_AO_AD_QA_GL.vcf -o filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz

bcftools index -t filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz
bcftools index -n filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz

echo "DONE"

------------- END OF THE SCRIPT ---------------------------------------------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 4046307 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode
# less wall time 

____________________________________________________________________________________________________________________________________________________________________


------------------- START OF INTERACTIVE JOB -----------------------------------------------------------------------------------------------------------------------------------
salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 00:30:00 --mem=2G
mamba activate bcftools.env

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/stats
awk '{print $1}' filt_samples_mtDNA_qc_stats.txt > samples.list.filt # remove the first row "sampleID" by hand
mv samples.list.filt ../samples.list.filt


cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA

# Run the loop using samtools faidx to isolate ONLY the mtDNA region
while IFS= read -r s; do
    samtools faidx /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/Reference_Genome/Formica_hybrid_v1_wFhyb_Sapis.fa mtDNA | \
    bcftools consensus -s "$s" -I filtered_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz > "${s}_mtDNA.fa"
    
    # Clean up the fasta header so it is just the sample name
    sed -i "s/>.*/>${s}/" "${s}_mtDNA.fa"
done < samples.list.filt

mkdir indiv_samples_mtDNA
mv *fa indiv_samples_mtDNA

find /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/indiv_samples_mtDNA -maxdepth 1 -type f | wc -l
# 140 files
------------------- END OF INTERACTIVE JOB --------------------------------------------------------------------------------------------------------------------------

31_filt.alignment_mtDNA.sh
-------------------- START OF BASH JOB -----------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash 
#SBATCH -J filt.alignment_mtDNA
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/filt.alignment_mtDNA.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/filt.alignment_mtDNA.err
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

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/indiv_samples_mtDNA

echo "Combining FASTA files"
# Concatenate all sample FASTAs into one file
cat *.fa > filt_mtDNAsamples.fa

echo "Alignment"
# alignment
## --thread -1 tells MAFT to use all allocated CPUs
## --anysymbol accept any valid text character which could occur because of missing data
mafft --anysymbol --thread -1 filt_mtDNAsamples.fa > filt_mtDNAsamples_aligned.fa

echo "DONE" 

-------------------- END OF BASH JOB ----------------------------------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 4046373 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode

________________________________________________________________________________________________________________________________________________________________

# https://github.com/josephhughes/Sequence-manipulation/blob/master/Fasta2Nexus.pl
# copying the script to transform fasta file into nexus file

touch Fasta2NEXUS.pl
vi Fasta2NEXUS.pl # copy the code

salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 00:15:00 --mem=2G
mamba activate mafft.env

perl Fasta2NEXUS.pl /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/indiv_samples_mtDNA/filt_mtDNAsamples_aligned.fa /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/indiv_samples_mtDNA/filt_mtDNAsamples_aligned.nex


# use in popart to build the network
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/indiv_samples_mtDNA/filt_mtDNAsamples_aligned.nex .
# use as a list of samples to filter the trait file
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/samples.list.filt .

# Popart was downloaded from https://popart.maths.otago.ac.nz/download/
## Popart workflow
# File > Open > all_mtDNAsamples_aligned.nex
### Epsilon value = 0
# Network > Minimum spanning network
# File > Export Graphics
# Network > Median joining network
# File > Export Graphics
# Statistics > Identical Sequences > Log to file? > Yes
# Statistics > All stats > Log to file? > Yes


_________________________________________________________________________________________________________________________________________________________________
## Dataset 2

# 33_vcf_filtering_mtDNA.sh
------------- START OF THE SCRIPT ---------------------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash -l
#SBATCH -J LMUF_mtDNA
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/vcf_LMUF.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/vcf_LMUF.err
#SBATCH -t 00:15:00
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
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
bcftools norm -f /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/Reference_Genome/Formica_hybrid_v1_wFhyb_Sapis.fa -m -both -O z -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/LMUF_samples_mtDNA_norm.vcf.gz /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/LMUF_samples_mtDNA_raw.vcf.gz

bcftools filter --threads 1 -Oz -s+ --SnpGap 2 LMUF_samples_mtDNA_norm.vcf.gz > LMUF_samples_mtDNA_norm.SnpGap_2.vcf.gz && \

bcftools filter --threads 1 -Oz -e 'TYPE!="snp"' -s NonSnp -m+ LMUF_samples_mtDNA_norm.SnpGap_2.vcf.gz > LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.vcf.gz && \

bcftools filter --threads 1 -Oz -s Balance -m+ -i 'RPL>=1 && RPR>=1 && SAF>=1 && SAR>=1' LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.vcf.gz > LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.vcf.gz && \

bcftools view --threads 1 -O z -f PASS LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.vcf.gz > LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.vcf.gz && \

bcftools view --threads 1 LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.vcf.gz | vcfallelicprimitives --keep-info --keep-geno -t decomposed | sed '/^##/! s/|/\//g' | sed 's/\.:\.:\.:\.:\.:\.:\.:\./\.\/\.:\.:\.:\.:\.:\.:\.:\./g' | bcftools sort --temp-dir $TMPDIR --max-mem 4G -O z > LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz && \


bcftools index -t LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz
echo "bcftools index -n LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz"
bcftools index -n LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz

#####
##### 1. SNP QUAL >= 30, biallelic -------------------------------------------------------------------
#####
echo "SNP QUAL >= 30, biallelic"

bcftools filter --threads 1 --include 'QUAL >= 30 && TYPE="snp"' -Oz LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz > LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz
echo "gunzip -c LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz | grep -vc '#'"
gunzip -c LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz | grep -vc '#'

bcftools view --threads 1 --min-alleles 2 --max-alleles 2 LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz -Oz > LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz
echo "gunzip -c LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz | grep -vc '#'"
gunzip -c LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz | grep -vc '#'
bcftools index -t LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz


#####
##### 2. Correcting the header --------------------------------------------------
#####
echo "Correcting the header"

# Extract header from VCF
bcftools view -h LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz > header.vcf

# Fix fields
perl -npe 's/<ID=AO,Number=A/<ID=AO,Number=\./' header.vcf | perl -npe 's/<ID=AD,Number=R/<ID=AD,Number=\./' | perl -npe 's/<ID=QA,Number=A/<ID=QA,Number=\./' | perl -npe 's/<ID=GL,Number=G/<ID=GL,Number=\./' > header_AO_AD_QA_GL.vcf

# Replace corrected header
bcftools reheader -h header_AO_AD_QA_GL.vcf -o LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz

bcftools index -t LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz
bcftools index -n LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz

echo "DONE"

------------- END OF THE SCRIPT ---------------------------------------------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 4052686 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode


____________________________________________________________________________________________________________________________________________________________________


------------------- START OF INTERACTIVE JOB -----------------------------------------------------------------------------------------------------------------------------------
salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 00:10:00 --mem=2G
mamba activate bcftools.env

## Check the number of SNPs across all samples
bcftools view -v snps LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz | grep -vc "^#"
# output: 247

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/stats
awk '{print $1}' LMUF_samples_mtDNA_qc_stats.txt > samples.list.LMUF # remove the first row "sampleID" by hand
mv samples.list.LMUF ../samples.list.LMUF


cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA

# Run the loop using samtools faidx to isolate ONLY the mtDNA region
while IFS= read -r s; do
    samtools faidx /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/Reference_Genome/Formica_hybrid_v1_wFhyb_Sapis.fa mtDNA | \
    bcftools consensus -s "$s" -I LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz > "${s}_mtDNA.fa"
    
    # Clean up the fasta header so it is just the sample name
    sed -i "s/>.*/>${s}/" "${s}_mtDNA.fa"
done < samples.list.LMUF


mv *fa indiv_samples_mtDNA

34_LMUF.alignment_mtDNA.sh
-------------------- START OF BASH JOB -----------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash 
#SBATCH -J LMUF.alignment_mtDNA
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/LMUF.alignment_mtDNA.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/LMUF.alignment_mtDNA.err
#SBATCH -t 00:30:00
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-user ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# set the script to stop immediately if any command fails or any part of a pipeline fails, preventing silent errors and corrupted results
set -eo pipefail

# Load environment
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/mafft.env

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/indiv_samples_mtDNA

echo "Combining FASTA files"
# Concatenate all sample FASTAs into one file
cat /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/indiv_samples_mtDNA/*.fa > LMUF_mtDNAsamples.fa

echo "Alignment"
# alignment
## --thread -1 tells MAFFT to use all allocated CPUs
## --anysymbol accept any valid text character which could occur because of missing data
mafft --anysymbol --thread -1 LMUF_mtDNAsamples.fa > LMUF_mtDNAsamples_aligned.fa

echo "DONE" 

-------------------- END OF BASH JOB ----------------------------------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 4052713 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode

________________________________________________________________________________________________________________________________________________________________

salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 00:10:00 --mem=2G
mamba activate mafft.env

perl Fasta2NEXUS.pl /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/indiv_samples_mtDNA/LMUF_mtDNAsamples_aligned.fa /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/indiv_samples_mtDNA/LMUF_mtDNAsamples_aligned.nex


# use in popart to build the network
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/indiv_samples_mtDNA/LMUF_mtDNAsamples_aligned.nex .
# use as a list of samples to filter the trait file
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/samples.list.LMUF .

________________________________________________________________________________________________________________________________________________________________________

# Checking missingness
salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 00:10:00 --mem=2G
mamba activate bcftools.env

vcftools --gzvcf LMUF_samples_mtDNA_norm.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz --missing-indv --out missing









