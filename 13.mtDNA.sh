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

mamba activate mafft.env

## Check the number of SNPs across all samples
bcftools view -v snps all_samples_mtDNA_raw.vcf | grep -vc "^#"
# output: 171

## Check the quality
bcftools view -v snps -i 'QUAL<30' all_samples_mtDNA_raw.vcf | grep -vc "^#"
# output: 123
bcftools view -v snps -i 'QUAL<20' all_samples_mtDNA_raw.vcf | grep -vc "^#"
# output: 121
bcftools view -v snps -i 'QUAL<10' all_samples_mtDNA_raw.vcf | grep -vc "^#"
# output: 117

## Check the depth
bcftools query -f '%CHROM\t%POS[\t%DP]\n' all_samples_mtDNA_raw.vcf \
| awk '{
    sum=0; n=0;
    for(i=3;i<=NF;i++){
        if($i!="." && $i>0){sum+=$i; n++}
    }
    if(n>0) print $1,$2,sum/n
}'

# output is contig name, position, Depth

# Might try also without filtering and compare the results
##### FILTERING OF mtDNA ######

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



## As the perl script did not work for me, I switched to bcftools consensus

# Run the loop using samtools faidx to isolate ONLY the mtDNA region
while IFS= read -r s; do
    samtools faidx /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/Reference_Genome/Formica_hybrid_v1_wFhyb_Sapis.fa mtDNA | \
    bcftools consensus -s "$s" -I all_samples_mtDNA_filtered.vcf.gz > "${s}_mtDNA.fa"
    
    # Clean up the fasta header so it is just the sample name
    sed -i "s/>.*/>${s}/" "${s}_mtDNA.fa"
done < samples.list

mkdir mtDNA_fa
mv *fa mtDNA_fa

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

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/mtDNA_fa

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










#####################  NO FILTERING STEPS #########################################################
####### SORTING, COMPRESSING, AND INDEXING #####

------------------- START OF INTERACTIVE JOB -----------------------------------------------------------------------------------------------------------------------------------
salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 00:30:00 --mem=2G

mamba activate bcftools.env

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA

# Sort and compress
bcftools sort -m 1G -Oz -o all_samples_mtDNA_raw.vcf.gz -T ./tmp_sort all_samples_mtDNA_raw.vcf

# Index
tabix -p vcf all_samples_mtDNA_raw.vcf.gz
bcftools index -n all_samples_mtDNA_raw.vcf.gz
# output: 187

#Extract the number of sites per mtDNA 
gunzip -c all_samples_mtDNA_raw.vcf.gz | grep -v "^#" | cut -f 1 | uniq -c > site_count_per_mtDNA.tab
# output: 187 mtDNA (means that there is 187 variant sites, some might be indels)

######### Create file for PopArt #########
## Next, we create a consensus sequence by by applying VCF variants to a reference fasta file using vcftools vcf-consensus scripts from
## https://github.com/vcftools/vcftools/blob/master/src/perl/vcf-consensus#L111
#Note, this is a pearl script. You have to copy-paste it to a file "vcf-consensus_code.pl" and make it executable with
chmod +x vcf-consensus_code.pl

# create a sample list 
bcftools query -l all_samples_mtDNA_raw.vcf.gz > samples.list


### to be edited

# forward loop, which iterates through the samples.list file and creates the consensu file
#NOTE. -H keeps the first haplotype, -s are the samples
#In one line:
while IFS= read -r s; do samtools faidx /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/Reference_Genome/Formica_hybrid_v1_wFhyb_Sapis.fa mtDNA | ./vcf-consensus_code.pl -H 1 -s "$s" all_samples_mtDNA_raw.vcf.gz > "${s}_mtDNA.fa"; done < samples.list
# doesn't work because of missing genotypes


###
#After this, you will fa.files for each sample. Move these files into a folder "mtDNA_fa", which you will create using mkdir
mkdir mtDNA_fa
mv *fa mtDNA_fa/

#Inside the folder, we now add all fa-files into a sample list
ls *fa > list.samples

------------------- END OF INTERACTIVE JOB -----------------------------------------------------------------------------------------------------------------------------------

# Popart file preparation

##############
#Now, we download the file locally to plot it in PopArt. PopArt needs a phylip format, while the data format is mafft.
scp krapfpat@puhti.csc.fi:/scratch/project_2009316/DutchSamples/X204SC23115958-Z01-F008/07.VCF_mtDNA/mtDNA_fa/all38_mtDNAsamples_aligned.fa ./ 


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

##END
