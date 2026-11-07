Adjusted from Ina Satokagas https://github.com/isatok/Wood_ant_local_genomic_patterns/blob/main/Pixy/1_allsites_vcf.sh

#Create needed directories:

mkdir /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF
mkdir /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/raw/
mkdir /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/logs/

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF

# get bam list
cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/bam_final/
#find . -type f -name "*.bam" -exec realpath {} \; > gvcf.bam.all.list 
## gvcf.bam.all.list  would be identical to all.bam.list used for SNP calling with freebayes
cp all.bam.list /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/

# get scaffold list
cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/Reference_Genome
grep -i '^>.*scaffold' Formica_hybrid_v1_wFhyb_Sapis.fa | grep -viE 'scaffold(00|03)' | sed 's/^>//' > scaffold_no0300.list


# 39.allsites_vcf_highcovsamples_mac2.sh
# ------- START OF BASH SCRIPT ----------------------------------------------------------------------------------------

#!/bin/bash -l
#SBATCH -J allsites_vcf_highcovsamples_mac2
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/logs/allsites_vcf_highcovsamples_mac2_%j_%a.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/logs/allsites_vcf_highcovsamples_mac2_%j_%a.err
#SBATCH -t 24:00:00
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --array=1-26%10
#SBATCH --ntasks 1
#SBATCH --mem=2G
#SBATCH --mail-user ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# set the script to stop immediately if any command fails or any part of a pipeline fails, preventing silent errors and corrupted results
set -eo pipefail

# Load environment
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/bcftools.env

# go to directory where raw gvcf for each scaffold will be produced
cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/raw/

# Get scaffold ID
REF=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/Reference_Genome
scaffold=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${REF}/scaffold_no0300.list)

# write gvcf
echo "Writing gvcf for ${scaffold} ..."


# bcftools call -m retains also monomorphic sites, which is needed for generation of all sites vcf
bcftools mpileup -f $REF/Formica_hybrid_v1_wFhyb_Sapis.fa \
  -b /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/all.bam.list \
  -r ${scaffold} | bcftools call -m -Oz -f GQ -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/raw/${scaffold}_allsamples.vcf.gz

# Set DP thresholds & filter missing data
#Max depth threshold: Calculate average of the mean depths and multiply by two -> follow the logic from the variant site vcf filtering where 
#sites with >2x mean ind depth are set as missing.
# IDEPTHPATH=/scratch/project_2001443/barriers_introgr_formica/vcf/filt
# cat $IDEPTHPATH/all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.idepth \
# | grep -v "INDV" | awk '{ total += $3; count++ } END { print total/count }' #15.6196
# 15.6196*2= ca. 31,24; set max avg depth to 31x.

mindp=386 #2*193
maxdp=11580 #60*193

echo "Filtering gvcf ..."
bcftools view -Ou ${scaffold}_allsamples.vcf.gz | bcftools filter -e "DP < ${mindp} | DP > ${maxdp}" --set-GTs . -Oz > ${scaffold}_allsamples_meanDP${mindp}-${maxdp}.vcf.gz &&
bcftools index -t ${scaffold}_allsamples_meanDP${mindp}-${maxdp}.vcf.gz

# Removing individuals which had in freebayes SNP calling more than 50% missingness
echo "Removing low coverage individuals ..."
vcftools --gzvcf ${scaffold}_allsamples_meanDP${mindp}-${maxdp}.vcf.gz  --remove-indv 121-Fprat --remove-indv RN417 --remove-indv 105-FaquH \
--remove-indv s353 --remove-indv 54-Frufa --remove-indv RN418 --remove-indv s354 \
--recode --recode-INFO-all --stdout | bgzip > ${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}.vcf.gz
bcftools index -t ${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}.vcf.gz

echo "Filtering missing data per site ..."
bcftools view -Ou ${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}.vcf.gz | bcftools filter -e "F_MISSING > 0.5" -Oz > \
${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}_maxNA50perc.vcf.gz


# two vcfs will be made, one with only invariant sites and one with only variants\
    # these two are then combine to make a vcf will all sites (gvcf)
# to combine they must both be indexed then combined with bcftools concat

echo "Creating vcf for invariant sites ..."
# Create a filtered VCF containing only invariant sites
vcftools --gzvcf ${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}_maxNA50perc.vcf.gz \
--max-maf 0 --recode --recode-INFO-all --stdout | bgzip -c > ${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}_maxNA50perc_invar.vcf.gz

echo "Creating vcf for variant sites ..."
# Create a filtered VCF containing only variant sites, keep SNPqual >= 30 and filter out HW excess and singletons
vcftools --gzvcf ${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}_maxNA50perc.vcf.gz --mac 2 --hwe 0.001 --minQ 30 --recode --recode-INFO-all --stdout | bgzip -c > \
${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}_maxNA50perc_var.vcf.gz

# Index both vcfs using tabix
tabix ${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}_maxNA50perc_invar.vcf.gz
tabix ${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}_maxNA50perc_var.vcf.gz

echo "Combining vcfs for ${scaffold} ..."
# Combine the two VCFs using bcftools concat
bcftools concat \
--allow-overlaps \
${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}_maxNA50perc_invar.vcf.gz \
${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}_maxNA50perc_var.vcf.gz \
-Oz > ${scaffold}_highcovsamples_filtered.vcf.gz

# Index for concat afterwards
tabix ${scaffold}_highcovsamples_filtered.vcf.gz


# Remove temp files
echo "Removing temp files ..."

# Remove intermediate files after DP filtering
rm -f ${scaffold}_allsamples_meanDP${mindp}-${maxdp}*
# rm -f Scaffold*_allsamples_meanDP*

# Remove intermediate files after missing data filtering
rm -f ${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}*
# rm -f Scaffold*_highcovsamples_meanDP*

echo "${scaffold} DONE"

# ------- END OF BASH SCRIPT ----------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 4233970 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode


########## NOT EDITED YET ###################

### Interactive script to combine the per scaffold gvcf files ### ----------------------------


# Combine VCFs
cd /scratch/project_2001443/barriers_introgr_formica/gvcf/raw/

#ls *allsamples_filtered.vcf.gz > all_filtered_vcfs.list
ls *_highcovsamples_filtered.vcf.gz > highcovsamples_filtered_vcfs.list

sinteractive...
module load biokit

bcftools concat \
--allow-overlaps \
-f highcovsamples_filtered_vcfs.list \
-Oz > ../highcovsamples_filtered.vcf.gz

tabix highcovsamples_filtered.vcf.gz


### Filter out duplicate sites from the all-sites vcf highcovsamples_filtered.vcf.gz ### ----------------------------


### remove_duplicates_highcov.sh ### 
#!/bin/bash -l
#SBATCH -J remove_duplicates_highcov
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/gvcf/logs/remove_duplicates_highcov.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/gvcf/logs/remove_duplicates_highcov.err
#SBATCH --account=project_2001443
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --ntasks 1
#SBATCH --mem=8G

module load biokit
cd /scratch/project_2001443/barriers_introgr_formica/gvcf

VCFIN=/scratch/project_2001443/barriers_introgr_formica/gvcf/highcovsamples_filtered.vcf.gz

echo "counting the initial number of sites in highcovsamples_filtered.vcf.gz ..."
bcftools index -n $VCFIN

#Filter out duplicate sites (achieve this by filtering out all that are mnp's and not homozygous "ref" or "snp" type)

echo "counting the amount of duplicates in highcovsamples_filtered.vcf.gz (divide this number by two)..."
echo $VCFIN
gunzip -c $VCFIN | cut -f2 | uniq -D | wc -l

#Keep only sites hom for reference allele or snps
bcftools filter --threads 8 -Oz -e 'TYPE!="ref" && TYPE!="snp"' -m+ $VCFIN > highcovsamples_filtered_noindels.vcf.gz
bcftools index -t highcovsamples_filtered_noindels.vcf.gz

echo "counting how many sites remain after TYPEref and TYPEsnp filtering..."
bcftools index -n highcovsamples_filtered_noindels.vcf.gz

echo "counting how many duplicates are left (divide the number by two)..."
gunzip -c highcovsamples_filtered_noindels.vcf.gz | cut -f2 | uniq -D | wc -l

echo "checking examples of what kind of variants they are - mnps? look at these coordinates in the highcovsamples_filtered_noindels.vcf.gz file..."
gunzip -c highcovsamples_filtered_noindels.vcf.gz | grep "Scaffold01" | cut -f2 | uniq -D | head     # e.g. xx, yy, zz
#bcftools view -H -r Scaffold01:5047119 highcovsamples_filtered_noindels.vcf.gz #They are mnps

#Always the first one is the snp and the second one mnp. Remove all second instances
echo "removing all second instances of duplicate sites (expecting the first ones to be the real snps)..."
bcftools norm --rm-dup all --threads 8 -Oz highcovsamples_filtered_noindels.vcf.gz > highcovsamples_filtered_noindels_rmdup.vcf.gz
bcftools index -t highcovsamples_filtered_noindels_rmdup.vcf.gz

echo "counting how many sites remain after duplicate filtering (add to this the initial #of duplicates - does it match with the input vcf site count?) ..."
bcftools index -n highcovsamples_filtered_noindels_rmdup.vcf.gz 

echo "counting how many duplicates remain (divide by two)..."
gunzip -c highcovsamples_filtered_noindels_rmdup.vcf.gz | cut -f2 | uniq -D | wc -l

### END.
