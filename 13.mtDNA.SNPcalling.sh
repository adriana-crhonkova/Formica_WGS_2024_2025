# Below, we will first call SNPs from the mitochonrial part of the WGS only
# For this, we will use again freebayes and need a "regions_mtDNA.txt" file with only the mtDNA region in it for the SNP calling


________________________________________________________________________________________________________________________________________________________________
# Creating regions_mtDNA.txt for mtDNA SNP calling 
------------------- START OF INTERACTIVE JOB -------------------------------------------------------------------------------------------------------------------
salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 04:00:00 --mem=2G

#modify regions file so that Freebayes accepts it
cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/Reference_Genome
awk 'BEGIN{OFS=""} {print $1,":",$2,"-",$3}' Formica_hybrid_v1_wFhyb_Sapis_5e6_data_regions.fa.fai > Formica_hybrid_v1_5e6_data_regions.tmp

#remove regions that require different SNP calling parameters if they are needed: Wolbachia (wFhyb*), and Spiroplasma (Spiroplasma*), as well as Scaffold00, but KEEP mtDNA
grep -v -e wFhyb -e Spiroplasma -e Scaffold00 Formica_hybrid_v1_5e6_data_regions.tmp > Formica_hybrid_v1_5e6_data_regions_wmtDNA.txt ; rm Formica_hybrid_v1_5e6_data_regions.tmp
grep -e mtDNA Formica_hybrid_v1_5e6_data_regions_wmtDNA.txt > regions_mtDNA.txt

------------------- END OF INTERACTIVE JOB ---------------------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 3649566 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode

__________________________________________________________________________________________________________________________________________________________________
# Quality control of mtDNA and filtering of samples with low coverage

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/
mkdir 04.VCF.mtDNA
cd 04.VCF.mtDNA
mkdir logs

## Quality Control for mtDNA
# going back to mosdepth statistics
cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/stats/coverage
# find lines with mtDNA in mosdepth summary files and prints them to one file (non-clipped overlaps)
grep -Hw "mtDNA" *merged.mosdepth.summary.txt > mosdepth_all_mtDNA_summary.txt

scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/stats/coverage/mosdepth_all_mtDNA_summary.txt .

### Running samtools coverage on all bam files (also the reference samples)
# 22_bam_mtDNA_qc.sh
---------------------------- START OF BASH JOB -----------------------------------------------------------------------------------------------------
#!/bin/bash 
#SBATCH -J bam_mtDNA_qc
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/logs/bam_mtDNA_qc.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/logs/bam_mtDNA_qc.err
#SBATCH -t 02:00:00
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-user ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# Load environment
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/samtools.env

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/bam_final

echo "START"

# Print the header first
echo -e "Sample\trname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq" > samtools_mtdna_coverage_summary.txt

echo "Running samtools coverage"

while IFS= read -r bam_path; do
    [ -z "$bam_path" ] && continue
    sample_name=$(basename "$bam_path" | cut -d'.' -f1)
    stats=$(samtools coverage -r mtDNA "$bam_path" | grep -v '^#')
    if [ ! -z "$stats" ]; then
        echo -e "${sample_name}\t${stats}" >> samtools_mtdna_coverage_summary.txt
    else
        echo -e "${sample_name}\tmtDNA_not_found_or_empty" >> samtools_mtdna_coverage_summary.txt
    fi
done < all.bam.list

mv samtools_mtdna_coverage_summary.txt /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/stats/coverage/samtools_mtdna_coverage_summary.txt

echo "FINISHED"

---------------------------- END OF BASH JOB -----------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 3990612 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode

       JobID    Elapsed     MaxRSS  AllocCPUS      State ExitCode
------------ ---------- ---------- ---------- ---------- --------
3990612        00:01:42                     1  COMPLETED      0:0
3990612.bat+   00:01:42     10880K          1  COMPLETED      0:0

# the wall time can be reduced next time

scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/bam_final/samtools_mtdna_coverage_summary.txt .


# Based on the samtools coverage statistics ("samtools_mtdna_coverage_summary.txt") I set a threshold of coverage breadth >70, depth >15X. 
# I filtered out 36 samples, all of them were sequenced in the batch Formica_2025 (only sequenced once). 
# Only sample from this batch which passed the filters was LMUF_00118c with coverage breadth 71.1283, depth 178.853.
# no sample from the merged batches of 2024 samples were filtered out, neither any of the reference samples. Seems there was some issue with library prep?


# directory with bam files
cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/bam_final
touch low_coverage_samples.txt
vi low_coverage_samples.txt # paste the names of samples to be filtered
# remove the samples and save into "filtered.bam.list"
grep -v -F -f low_coverage_samples.txt all.bam.list > mtDNA.coverage.PASS.bam.list
# Sanity check: 193 (lines in all.bam.list) - 157 (lines in mtDNA.coverage.PASS.bam.list) = 36 (number of samples which should be filtered)
touch unknown.samples.txt
vi unknown.samples.txt # paste RNxxx and F. exsecta samples

grep -v -F -f unknown.samples.txt mtDNA.coverage.PASS.bam.list > mtDNA.coverage.PASS.rm.unknown.bam.list
# Sanity check: 157 (lines in mtDNA.coverage.PASS.bam.list) - 143 (lines in mtDNA.coverage.PASS.rm.unknown.bam.list) = 14 (number of samples which should be filtered)

touch exsecta.samples.txt
vi exsecta.samples.txt
grep -v -F -f exsecta.samples.txt mtDNA.coverage.PASS.rm.unknown.bam.list > mtDNA.final.bam.list
# Sanity check: 143 (lines in mtDNA.coverage.PASS.rm.unknown.bam.list) - 140 (lines in mtDNA.final.bam.list) = 3 (number of samples which should be filtered)
# there are 4 F. exsecta samples, but LMUF_00427b was already filtered for low coverage of mtDNA

_______________________________________________________________________________________________________________________________________________________________
# SNP CALLING


# For nuclear DNA SNP calling I used following parameters
# --skip-coverage was calculated as: 193 samples, max coverage across all samples 46.4X -> 46.4*193=8955.2 -> --skip-coverage 10000 will only skip super high-depth artifacts
# using the same logic: 140 samples, max coverage across mtDNA of all samples 5874X -> 5874*140=822 360 -> --skip-coverage 850000

# Checked scripts of Ina's paper "Semipermeable species boundaries ..." (Satokangas et. al, 2023) 

# /appl/soft/bio/bioconda/miniconda3/envs/freebayes/bin/freebayes-puhti \
#  -time 72 \
#  -regions ref/mtDNA_50kb_regions.txt \
#  -f ref/mtDNA_wFhyb_Sapis.fa \
#  -L vcf/bam.list \
#  -k --genotype-qualities --pooled-continuous -F 0.05 --ploidy 1 --min-coverage 1000 \ 
#  -out vcf/mt_F005_C1000_all_samples_raw.vcf


---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Here is updated SNP calling based on Ina's scripts for paper "Semipermeable species boundaries create opportunities for gene flow and adaptive potential"
# https://zenodo.org/records/7941711?preview_file=Raw_reads_to_SNPs_scripts.zip


# 29.SNPcalling.filt.mtDNA.sh
-------------------- START OF BASH JOB -----------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash 
#SBATCH -J SNPcalling.filt.mtDNA.2
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/SNPcalling.filt.mtDNA.2.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/SNPcalling.filt.mtDNA.2.err
#SBATCH -t 48:00:00
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-user ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# set the script to stop immediately if any command fails or any part of a pipeline fails, preventing silent errors and corrupted results
set -eo pipefail

# Load environment
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/freebayes_updated.env

REF=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/Reference_Genome
RES=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA
BAM=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/bam_final

echo "Started SNP calling"

freebayes-parallel \
  $REF/regions_mtDNA.txt \
  $SLURM_CPUS_PER_TASK \
  -f $REF/Formica_hybrid_v1_wFhyb_Sapis.fa \
  -L $BAM/mtDNA.final.bam.list \
  -k \
  --genotype-qualities \
  --pooled-continuous \
  -F 0.05 \
  --ploidy 1 \
  --skip-coverage 850000 \
  --use-best-n-alleles 3 \
  > "$RES/filtered_samples_mtDNA_raw.vcf"

# removed --min-coverage 1000 \ flag, because 5 samples (even after filtering) have mean coverage depth <1000 (eg. LMUF_00118c with coverage breadth 71.1283, mean depth 178.853X, max depth 1263X).

echo "Finished SNP calling"

-------------------- END OF BASH JOB -----------------------------------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 4046250 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode
# use less wall time

___________________________________________________________________________________________________________________________________________________________________
# Quality Control of called SNPs

-------------------- START OF INTERACTIVE JOB -----------------------------------------------------------------------------------------------------------------------------------
salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 00:30:00 --mem=2G

mamba activate bcftools.env

####### Quality Check #######

# Check the number of variants
grep -vc "^#" filtered_samples_mtDNA_raw.vcf
# output: 1095

## Check the number of SNPs across all samples
bcftools view -v snps filtered_samples_mtDNA_raw.vcf | grep -vc "^#"
# output: 1016

## Check the quality (Quality 30 = 1 in 1,000 chance that the variant call is wrong; Quality 20 = 1 in 100, Quality 10 = 1 in 10)
bcftools view -v snps -i 'QUAL<30' filtered_samples_mtDNA_raw.vcf | grep -vc "^#"
# output: 590

## Check the depth per site
bcftools query -f '%CHROM\t%POS[\t%DP]\n' filtered_samples_mtDNA_raw.vcf \
| awk '{
    sum=0; n=0;
    for(i=3;i<=NF;i++){
        if($i!="." && $i>0){sum+=$i; n++}
    }
    if(n>0) print $1,$2,sum/n
}' > average_depth_per_position_filt.txt
# output is contig name, position, Depth
mv average_depth_per_position_filt.txt stats/average_depth_per_position_filt.txt


## bcftools stats
# Create a header
echo -e "Sample_ID\tnHapRef\tnHapAlt\tnMissing" > filt_samples_mtDNA_qc_stats.txt

# Get the haploid columns ($12, $13, and $14)
bcftools stats -s - filtered_samples_mtDNA_raw.vcf | grep "^PSC" | awk '{print $3"\t"$12"\t"$13"\t"$14}' >> filt_samples_mtDNA_qc_stats.txt
mv filt_samples_mtDNA_qc_stats.txt stats/filt_samples_mtDNA_qc_stats.txt
####### SORTING, COMPRESSING, AND INDEXING #####

# In the next step we "Sort" and "compress" the VCF file in one step.
bcftools sort -m 1G -Oz -o filtered_samples_mtDNA_raw.vcf.gz -T ./tmp_sort filtered_samples_mtDNA_raw.vcf

# Index
tabix -p vcf filtered_samples_mtDNA_raw.vcf.gz
bcftools index -n filtered_samples_mtDNA_raw.vcf.gz
# 1095 


_____________________________________________________________________________________________________________________________________________________________
## DATASET 2
# Running all LMUF samples without references

# prepare list of LMUF samples
cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/bam_final
grep -v "Ref" all.bam.list > LMUF.list

# 89 samples, max coverage across mtDNA of all samples 5874X -> 5874*89=522 786 -> --skip-coverage 550000
# 32.SNPcalling.mtDNA.LMUF.sh
-------------------- START OF BASH JOB -----------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash 
#SBATCH -J SNPcalling.mtDNA.LMUF
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/SNPcalling.mtDNA.LMUF.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/SNPcalling.mtDNA.LMUF.out
#SBATCH -t 00:45:00
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-user ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# set the script to stop immediately if any command fails or any part of a pipeline fails, preventing silent errors and corrupted results
set -eo pipefail

# Load environment
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/freebayes_updated.env

REF=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/Reference_Genome
RES=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA
BAM=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/bam_final

echo "Started SNP calling"

freebayes-parallel \
  $REF/regions_mtDNA.txt \
  $SLURM_CPUS_PER_TASK \
  -f $REF/Formica_hybrid_v1_wFhyb_Sapis.fa \
  -L $BAM/LMUF.list \
  -k \
  --genotype-qualities \
  --pooled-continuous \
  -F 0.05 \
  --ploidy 1 \
  --skip-coverage 550000 \
  --use-best-n-alleles 3 \
  > "$RES/LMUF_samples_mtDNA_raw.vcf"

# removed --min-coverage 1000 \ flag, because 5 samples (even after filtering) have mean coverage depth <1000 (eg. LMUF_00118c with coverage breadth 71.1283, mean depth 178.853X, max depth 1263X).

echo "Finished SNP calling"

-------------------- END OF BASH JOB -----------------------------------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 4052192 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode

___________________________________________________________________________________________________________________________________________________________________
# Quality Control of called SNPs

-------------------- START OF INTERACTIVE JOB -----------------------------------------------------------------------------------------------------------------------------------
salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 00:30:00 --mem=2G

mamba activate bcftools.env

####### Quality Check #######

# Check the number of variants
grep -vc "^#" LMUF_samples_mtDNA_raw.vcf
# output: 957

## Check the number of SNPs across all samples
bcftools view -v snps LMUF_samples_mtDNA_raw.vcf | grep -vc "^#"
# output: 798

## Check the quality (Quality 30 = 1 in 1,000 chance that the variant call is wrong; Quality 20 = 1 in 100, Quality 10 = 1 in 10)
bcftools view -v snps -i 'QUAL<30' LMUF_samples_mtDNA_raw.vcf | grep -vc "^#"
# output: 273

## Check the depth per site
bcftools query -f '%CHROM\t%POS[\t%DP]\n' LMUF_samples_mtDNA_raw.vcf \
| awk '{
    sum=0; n=0;
    for(i=3;i<=NF;i++){
        if($i!="." && $i>0){sum+=$i; n++}
    }
    if(n>0) print $1,$2,sum/n
}' > average_depth_per_position_LMUF.txt
# output is contig name, position, Depth
mv average_depth_per_position_LMUF.txt stats/average_depth_per_position_LMUF.txt


## bcftools stats
# Create a header
echo -e "Sample_ID\tnHapRef\tnHapAlt\tnMissing" > LMUF_samples_mtDNA_qc_stats.txt

# Get the haploid columns ($12, $13, and $14)
bcftools stats -s - LMUF_samples_mtDNA_raw.vcf | grep "^PSC" | awk '{print $3"\t"$12"\t"$13"\t"$14}' >> LMUF_samples_mtDNA_qc_stats.txt
mv LMUF_samples_mtDNA_qc_stats.txt stats/LMUF_samples_mtDNA_qc_stats.txt
####### SORTING, COMPRESSING, AND INDEXING #####

# In the next step we "Sort" and "compress" the VCF file in one step.
bcftools sort -m 1G -Oz -o LMUF_samples_mtDNA_raw.vcf.gz -T ./tmp_sort LMUF_samples_mtDNA_raw.vcf

# Index
tabix -p vcf LMUF_samples_mtDNA_raw.vcf.gz
bcftools index -n LMUF_samples_mtDNA_raw.vcf.gz
# 957




