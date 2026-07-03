### This script is edited from the one made by Patrick Krapf (based on Ina Satokangas?)

# Below, we will first call SNPs from the mitochonrial part of the WGS only
# For this, we nwill use again freebayes and need a "regions_mtDNA.txt" file with only the mtDNA region in it for the SNP calling

------------------- START OF INTERACTIVE JOB -----------------------------------------------------------------------------------------------------------------------------------
salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 04:00:00 --mem=2G

#modify regions file so that Freebayes accepts it
cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/Reference_Genome
awk 'BEGIN{OFS=""} {print $1,":",$2,"-",$3}' Formica_hybrid_v1_wFhyb_Sapis_5e6_data_regions.fa.fai > Formica_hybrid_v1_5e6_data_regions.tmp

#remove regions that require different SNP calling parameters if they are needed: Wolbachia (wFhyb*), and Spiroplasma (Spiroplasma*), as well as Scaffold00, but KEEP mtDNA
grep -v -e wFhyb -e Spiroplasma -e Scaffold00 Formica_hybrid_v1_5e6_data_regions.tmp > Formica_hybrid_v1_5e6_data_regions_wmtDNA.txt ; rm Formica_hybrid_v1_5e6_data_regions.tmp
grep -e mtDNA Formica_hybrid_v1_5e6_data_regions_wmtDNA.txt > regions_mtDNA.txt

------------------- END OF INTERACTIVE JOB -----------------------------------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 3649566 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode

##### Call SNPs #####
cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/
mkdir 04.VCF.mtDNA
cd 04.VCF.mtDNA
mkdir logs


# 20.SNPcalling.mtDNA.sh
-------------------- START OF BASH JOB -----------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash 
#SBATCH -J SNPcalling.mtDNA
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/SNPcalling.mtDNA.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/SNPcalling.mtDNA.err
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
  -L $BAM/all.bam.list \
  -k \
  --genotype-qualities \
  --skip-coverage 15200 \
  --limit-coverage 100 \
  --use-best-n-alleles 3 \
  > "$RES/all_samples_mtDNA_raw.vcf"

echo "Finished SNP calling"

-------------------- END OF BASH JOB -----------------------------------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 3649568 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode

------------------- START OF INTERACTIVE JOB -----------------------------------------------------------------------------------------------------------------------------------
salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 00:30:00 --mem=2G

mamba activate bcftools.env

## Check the number of SNPs across all samples
bcftools view -v snps all_samples_mtDNA_raw.vcf | grep -vc "^#"
# output: 171

## Check the quality (Quality 30 = 1 in 1,000 chance that the variant call is wrong; Quality 20 = 1 in 100, Quality 10 = 1 in 10)
bcftools view -v snps -i 'QUAL<30' all_samples_mtDNA_raw.vcf | grep -vc "^#"
# output: 123
bcftools view -v snps -i 'QUAL<20' all_samples_mtDNA_raw.vcf | grep -vc "^#"
# output: 121
bcftools view -v snps -i 'QUAL<10' all_samples_mtDNA_raw.vcf | grep -vc "^#"
# output: 117

## Check the tepth per site
bcftools query -f '%CHROM\t%POS[\t%DP]\n' all_samples_mtDNA_raw.vcf \
| awk '{
    sum=0; n=0;
    for(i=3;i<=NF;i++){
        if($i!="." && $i>0){sum+=$i; n++}
    }
    if(n>0) print $1,$2,sum/n
}' > average_depth_per_position.txt

# output is contig name, position, Depth

## Check the depth per site with only SNP
bcftools query -v snps -f '%CHROM\t%POS[\t%DP]\n' all_samples_mtDNA_raw.vcf \
| awk 'BEGIN{OFS="\t"} {
    sum=0; n=0;
    for(i=3;i<=NF;i++){
        if($i!="." && $i>0){sum+=$i; n++}
    }
    if(n>0) print $1,$2,sum/n
}' > average_snp_depth_per_position.tsv
# output: 154 lines
# the output have less lines then the output of bcftools view -v snps all_samples_mtDNA_raw.vcf | grep -vc "^#" (output: 171)

## bcftools stats
# Add a header row first, then append the parsed columns
echo -e "Sample_ID\tnRefHom\tnNonRefHom\tnHets\tnTransitions\tnTransversions\tnMissing" > sample_mtDNA_qc_stats_nofilt.txt
# run the stats
bcftools stats -s - all_samples_mtDNA_raw.vcf.gz | grep "^PSC" | awk '{print $3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$14}' >> sample_mtDNA_qc_stats_nofilt.txt

scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/average_depth_per_position.txt .
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/average_snp_depth_per_position.tsv .
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/sample_mtDNA_qc_stats_nofilt.txt .


# A)
################################ FILTERING OF mtDNA SNPs by Quality ############################################################################################

# Filter for QUAL >= 30 and keep only SNPs
bcftools view -v snps -i 'QUAL>=30' all_samples_mtDNA_raw.vcf -O z -o all_samples_mtDNA_filtered.vcf.gz

# Index the new file 
bcftools index all_samples_mtDNA_filtered.vcf.gz

## Check the number of SNPs across all samples
bcftools view -v snps all_samples_mtDNA_filtered.vcf.gz | grep -vc "^#"
# output 48

## Check the number of samples
bcftools query -l all_samples_mtDNA_filtered.vcf.gz | wc -l
# output: 193

# Checking the stats of retained SNPs
bcftools stats -s - all_samples_mtDNA_filtered.vcf.gz | grep "^PSC"

## saving the stats in a file
# Add a header row first, then append the parsed columns
echo -e "Sample_ID\tnRefHom\tnNonRefHom\tnHets\tnTransitions\tnTransversions\tnMissing" > sample_mtDNA_qc_stats.txt
# run the stats
bcftools stats -s - all_samples_mtDNA_filtered.vcf.gz | grep "^PSC" | awk '{print $3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$14}' >> sample_mtDNA_qc_stats.txt

scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/sample_mtDNA_qc_stats.txt .


## As the perl script did not work for me, I switched to bcftools consensus

# Run the loop using samtools faidx to isolate ONLY the mtDNA region
while IFS= read -r s; do
    samtools faidx /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/Reference_Genome/Formica_hybrid_v1_wFhyb_Sapis.fa mtDNA | \
    bcftools consensus -s "$s" -I all_samples_mtDNA_filtered.vcf.gz > "${s}_mtDNA.fa"
    
    # Clean up the fasta header so it is just the sample name
    sed -i "s/>.*/>${s}/" "${s}_mtDNA.fa"
done < samples.list

mkdir mtDNA_fa_filt
mv *fa mtDNA_fa_filt

------------------- END OF INTERACTIVE JOB --------------------------------------------------------------------------------------------------------------------------
21.alignment_mtDNA.sh
-------------------- START OF BASH JOB -----------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash 
#SBATCH -J alignment_mtDNA
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/alignment_mtDNA.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/alignment_mtDNA.err
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

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/mtDNA_fa_filt

echo "Combining FASTA files"
# Concatenate all sample FASTAs into one file
cat *.fa > all_mtDNAsamples.fa

echo "Alignment"
# alignment
## --thread -1 tells MAFT to use all allocated CPUs
## --anysymbol accept any valid text character which could occur because of missing data
mafft --anysymbol --thread -1 all_mtDNAsamples.fa > all_mtDNAsamples_aligned.fa

echo "DONE" 

-------------------- END OF BASH JOB ----------------------------------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 3948711 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode

       JobID    Elapsed     MaxRSS  AllocCPUS      State ExitCode
------------ ---------- ---------- ---------- ---------- --------
3948711        00:00:24                     4  COMPLETED      0:0
3948711.bat+   00:00:24                     4  COMPLETED      0:0

scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/mtDNA_fa/all_mtDNAsamples_aligned.fa .
# popart needs different format

# https://github.com/josephhughes/Sequence-manipulation/blob/master/Fasta2Nexus.pl
# copying the script to transform fasta file into nexus file

touch Fasta2NEXUS.pl
vi Fasta2NEXUS.pl # copy the code


salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 00:15:00 --mem=2G
mamba activate mafft.env

perl Fasta2NEXUS.pl /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/mtDNA_fa/all_mtDNAsamples_aligned.fa /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/mtDNA_fa/all_mtDNAsamples_aligned.nex

scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/mtDNA_fa/all_mtDNAsamples_aligned.nex .


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

# there should be a way how to add matedata to colour the nodes

############################# FROM PATRICK #################################################################################################

##WORK IN PROGRESS BELOW
#So, we need to change the format manually using the perl script "Fasta2Phylip.pl" from https://indra.mullins.microbiol.washington.edu/perlscript/docs/Sequence.html
#After downloading the script (do not forget to "chmod +x script), we run the script locally under 
#Linux - Root - Genomic data - Fasta2Phylip
#The code is below:
Fasta2Phylip.pl all38_mtDNAsamples_aligned.fa all38_mtDNAsamples_aligned.phylip

#The file "all38_mtDNAsamples_aligned.phylip" is now in a phylip format
#https://scikit.bio/docs/latest/generated/skbio.io.format.phylip.html
#and can be used as it is, but the suffix needs to be changed from "phylip" to "nex" (For more details, check the tutorial: https://dcsoto.github.io/2024/09/30/haplotype-network
all257_samples_aligned.nex
#This can be done locally. You can also use PGDspider to transform it! Then, you can use it in PopArt.

#### PopArt ####
# In PopArt, you load the "nex" file via the "File-Open"Function. Then you create a "Median Joining spanning network". Here you are asked to give an "Epsilon"-value. "Epsilon defines a threshold for allowing connections that aren't the shortest (most parsimonious) path between nodes (haplotypes). By setting a higher epsilon, you permit more connections, allowing the network to capture more complexity, especially when evolutionary relationships are unclear or reticulate (branching and reconnecting). (-> from CHATGPT). 

# You can also add a "trait.txt" file, with which you can color-code the nodes based on species identity. In this file, each sample is in one line and starts with the species ID followed by a numeric code of 0 and 1 indicaitng to which species it belongs.
#Note, the species IDs and the IDs in the phylip file need to be identical for this to work
#For example 0,1,0,0,0,0,0,0 means it belongs to F aquilonia. For simplicity, I've added part of such a file below:
AQU,LUG,RUF,POL,PRA,HYB1,HYB2,HYB3
102-Frufa_mtDNA.fa,0,0,1,0,0,0,0,0
108-Flug_mtDNA.fa,0,1,0,0,0,0,0,0
115-Flug_mtDNA.fa,0,1,0,0,0,0,0,0
117-Fprat_mtDNA.fa,0,0,0,0,1,0,0,0
119-Flug_mtDNA.fa,0,0,0,0,0,0,1,0
120-Fprat_mtDNA.fa,0,0,0,0,1,0,0,0
126-Flug_mtDNA.fa,0,0,0,0,0,0,1,0
52-Fpol_mtDNA.fa,0,0,0,0,0,1,0,0
72-Frufa_mtDNA.fa,0,0,1,0,0,0,0,0

#We can load both files into PopArt. First the phylip-alignment, then the trait file. Do NOT delete the alignment when loading the file.
#First run the minimum and then the median spanning network

#####################################################################################################################################################

### Troubleshooting

# After noticing that some samples in "sample_mtDNA_qc_stats_nofilt.txt" and "sample_mtDNA_qc_stats.txt" (filtered Q>30) have 0 variants.
# As a follow-up, I checked a representative of a bad sample (LMUF_00341c) and of a good sample (LMUF_00206a) final bam file in IGV (Interactive genome viewer).
# There was clearly visible that the coverage breath for mtDNA scaffold in LMUF_00341c is very poor in comparison with LMUF_00206a. 
  # note: LMUF_00206a is from the merged batch of sequencing
# Therefore, I am going to do quality control again of bam files with focus on mtDNA. Before, I used mosdepth summary for coverage calculation accross all samples. 

# going back to mosdepth statistics
cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/stats/coverage
# find lines with mtDNA in mosdepth summary files and prints them to one file (non-clipped overlaps)
grep -Hw "mtDNA" *merged.mosdepth.summary.txt > combined_mtDNA_summary.txt

scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/stats/coverage/combined_mtDNA_summary.txt .

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
echo -e "Sample\trname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq" > mtdna_coverage_summary.txt

echo "Running samtools coverage"

while IFS= read -r bam_path; do
    [ -z "$bam_path" ] && continue
    sample_name=$(basename "$bam_path" | cut -d'.' -f1)
    stats=$(samtools coverage -r mtDNA "$bam_path" | grep -v '^#')
    if [ ! -z "$stats" ]; then
        echo -e "${sample_name}\t${stats}" >> mtdna_coverage_summary.txt
    else
        echo -e "${sample_name}\tmtDNA_not_found_or_empty" >> mtdna_coverage_summary.txt
    fi
done < all.bam.list

echo "FINISHED"

---------------------------- END OF BASH JOB -----------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 3990612 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode

       JobID    Elapsed     MaxRSS  AllocCPUS      State ExitCode
------------ ---------- ---------- ---------- ---------- --------
3990612        00:01:42                     1  COMPLETED      0:0
3990612.bat+   00:01:42     10880K          1  COMPLETED      0:0

# the wall time can be reduced next time

scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/bam_final/mtdna_coverage_summary.txt .


# Based on the samtools coverage statistics ("mtdna_coverage_summary.txt") I set a threshold of coverage breadth >70, depth >15X. 
# I filtered out 36 samples, all of them were sequenced in the batch Formica_2025 (only sequenced once). 
# Only sample from this batch which passed the filters was LMUF_00118c with coverage breadth 71.1283, depth 178.853.
# no sample from the merged batches of 2024 samples were filtered out, neither any of the reference samples. Seems there was some issue with library prep?

##### SNP calling with filtered dataset
# directory with bam files
cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/bam_final
touch remove_samples.txt
vi remove_samples.txt # paste the names of samples to be filtered
# remove the samples and save into "filtered.bam.list"
grep -v -F -f remove_samples.txt all.bam.list > filtered.bam.list
# Sanity check: 193 (lines in all.bam.list) - 157 (lines in filtered.bam.list) = 36 (number of samples which should be filtered)


# 23.SNPcalling.filt.mtDNA.sh
-------------------- START OF BASH JOB -----------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash 
#SBATCH -J SNPcalling.filt.mtDNA
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/SNPcalling.filt.mtDNA.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/SNPcalling.filt.mtDNA.err
#SBATCH -t 00:30:00
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G
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
  -L $BAM/filtered.bam.list \
  -k \
  --ploidy 1 \
  --genotype-qualities \
  --skip-coverage 15200 \
  --limit-coverage 100 \
  --use-best-n-alleles 3 \
  > "$RES/filt_samples_mtDNA_raw.vcf"

echo "Finished SNP calling"

-------------------- END OF BASH JOB -----------------------------------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 3994915 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode

       JobID    Elapsed     MaxRSS  AllocCPUS      State ExitCode
------------ ---------- ---------- ---------- ---------- --------
3994915        00:01:19                    16  COMPLETED      0:0
3994915.bat+   00:01:19   6184456K         16  COMPLETED      0:0



------------------- START OF INTERACTIVE JOB -----------------------------------------------------------------------------------------------------------------------------------
salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 00:30:00 --mem=2G

mamba activate bcftools.env

# Check the number of variants
grep -vc "^#" filt_samples_mtDNA_raw.vcf
# output: 187

## Check the number of SNPs across all samples
bcftools view -v snps filt_samples_mtDNA_raw.vcf | grep -vc "^#"
# output: 171

## Check the quality (Quality 30 = 1 in 1,000 chance that the variant call is wrong; Quality 20 = 1 in 100, Quality 10 = 1 in 10)
bcftools view -v snps -i 'QUAL<30' filt_samples_mtDNA_raw.vcf | grep -vc "^#"
# output: 124
bcftools view -v snps -i 'QUAL<20' filt_samples_mtDNA_raw.vcf | grep -vc "^#"
# output: 123
bcftools view -v snps -i 'QUAL<10' filt_samples_mtDNA_raw.vcf | grep -vc "^#"
# output: 123

## Check the depth per site
bcftools query -f '%CHROM\t%POS[\t%DP]\n' filt_samples_mtDNA_raw.vcf \
| awk '{
    sum=0; n=0;
    for(i=3;i<=NF;i++){
        if($i!="." && $i>0){sum+=$i; n++}
    }
    if(n>0) print $1,$2,sum/n
}' > average_depth_per_position_filt.txt

# output is contig name, position, Depth

## bcftools stats
# Create a header
echo -e "Sample_ID\tnHapRef\tnHapAlt\tnMissing" > filt_samples_mtDNA_qc_stats.txt

# Get the haploid columns ($12, $13, and $14)
bcftools stats -s - filt_samples_mtDNA_raw.vcf | grep "^PSC" | awk '{print $3"\t"$12"\t"$13"\t"$14}' >> filt_samples_mtDNA_qc_stats.txt
------------------- END OF INTERACTIVE JOB -----------------------------------------------------------------------------------------------------------------------------------
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/filt_samples_mtDNA_qc_stats.txt .

# Some samples have high amount of missing sites after SNP filtering (RN417 have 187 missing sites out of 187 sites - 100% missingness)

------------------- START OF INTERACTIVE JOB -----------------------------------------------------------------------------------------------------------------------------------
salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 00:30:00 --mem=2G
mamba activate bcftools.env

bcftools query -f '[%SAMPLE\t%F_MISSING\n]' filt_samples_mtDNA_raw.vcf | sort -u > sample_missingness.txt

# check the number of samples passing various thresholds
for threshold in 0.10 0.20 0.30 0.50; do
    remaining=$(awk -v t="$threshold" '$3 <= t {count++} END {print count+0}' sample_missingness.txt)
    # Use awk to calculate the clean percentage (e.g., 0.10 -> 10%)
    percent=$(awk -v t="$threshold" 'BEGIN {print t*100}')
    echo -e "<= ${percent}%\t\t\t${remaining}"
done

## Thresholds of maximal missingness
# 10%                  18
# 20%                  54
# 30%                  74
# 50%                  140

# at 50% threshld, only 17 samples will be filtered out

# make a list of samples to keep
awk '$3 < 0.50 {print $1}' sample_missingness.txt > samples_to_keep.txt
# sanity check: 140 lines

# filter vcf file
bcftools view -S samples_to_keep.txt filt_samples_mtDNA_raw.vcf > filt_samples_mtDNA_VCFfilt.vcf


####### SORTING, COMPRESSING, AND INDEXING #####

# In the next step we "Sort" and "compress" the VCF file in one step.
bcftools sort -m 1G -Oz -o filt_samples_mtDNA_VCFfilt.vcf.gz -T ./tmp_sort filt_samples_mtDNA_VCFfilt.vcf

# Index
tabix -p vcf filt_samples_mtDNA_VCFfilt.vcf.gz
bcftools index -n filt_samples_mtDNA_VCFfilt.vcf.gz


# Run the loop using samtools faidx to isolate ONLY the mtDNA region
while IFS= read -r s; do
    samtools faidx /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/Reference_Genome/Formica_hybrid_v1_wFhyb_Sapis.fa mtDNA | \
    bcftools consensus -s "$s" -I filt_samples_mtDNA_VCFfilt.vcf.gz > "${s}_mtDNA.fa"
    
    # Clean up the fasta header so it is just the sample name
    sed -i "s/>.*/>${s}/" "${s}_mtDNA.fa"
done < samples_to_keep.txt

mkdir filt_samples_mtDNA_fa
mv *fa filt_samples_mtDNA_fa

------------------- END OF INTERACTIVE JOB --------------------------------------------------------------------------------------------------------------------------
24_filt.alignment_mtDNA.sh
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

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/filt_samples_mtDNA_fa

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
sacct -M biohpc_gen -j 3995659 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode


salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 00:15:00 --mem=2G
mamba activate mafft.env

perl Fasta2NEXUS.pl /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/filt_samples_mtDNA_fa/filt_mtDNAsamples_aligned.fa /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/filt_samples_mtDNA_fa/filt_mtDNAsamples_aligned.nex

scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/filt_samples_mtDNA_fa/filt_mtDNAsamples_aligned.nex .

