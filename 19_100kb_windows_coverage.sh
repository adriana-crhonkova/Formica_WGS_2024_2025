# Adjusted from https://github.com/isatok/Wood_ant_local_genomic_patterns/blob/main/Pixy/allsites_depth.sh
cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM
mkdir cov_highcovInds_5perSp
mkdir cov_highcovInds_5perSp/logs

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/groupfiles
cp Fin_5ind_species_fixed.tab /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/Fin_5ind_species_fixed.tab

# generating the nodupl_highcovinds.list of .bam files used in pixy
## final bam files have "Ref" prefix for references + some Finnish references have species suffix - the following code deals with that and generates the list for the following script.

awk '{
    id = $1
    if (id ~ /^[0-9]+-/) {
        sub(/-.*/, "", id)
    }
    print "Ref_" id
}' Fin_5ind_species_fixed.tab > nodupl_highcovinds.list

## The list
Ref_104
Ref_44
Ref_107
Ref_109
Ref_106
Ref_108
Ref_124
Ref_113
Ref_11
Ref_56
Ref_CBCH1_1w
Ref_CBCH2_2w
Ref_VDa_1w
Ref_CBCH3_1w
Ref_NAZa_1w
Ref_16
Ref_37
Ref_102
Ref_26
Ref_72


### Compute sequencing depth per sample (after duplicate removal)

# 49.mosdepth_highcovinds.sh
# ------- START OF BASH SCRIPT ----------------------------------------------------------------------------------------

#!/bin/bash -l
#SBATCH -J mosdepth_highcovinds
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/logs/mosdepth_highcovinds_%j.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/logs/mosdepth_highcovinds_%j.err
#SBATCH -t 04:00:00
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --array=1-20
#SBATCH --ntasks 1
#SBATCH --mem-per-cpu=4GB
#SBATCH --mail-user ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# set the script to stop immediately if any command fails or any part of a pipeline fails, preventing silent errors and corrupted results
set -eo pipefail

# Load environment mosdepth
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/mosdepth


BAMDIR=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/bam_final/
FINALDIR=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp
NAMELIST=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/nodupl_highcovinds.list

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp

# Get sample ID
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p $NAMELIST)

echo "Processing sample: ${sample}"

mosdepth -t 1 -b 100000 -n $FINALDIR/${sample}_overlap_correction $BAMDIR/$sample"_nodupl_wRG_clip.bam"

echo "Finished sample: ${sample}"

# ------- END OF BASH SCRIPT ----------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 4287683 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode

# Checking which samples were problematic
for file in mosdepth_highcovinds_*.out; do
    echo -e "$(basename "$file")\t$(head -n 1 "$file")"
done > job_sample_mapping.txt

mosdepth_highcovinds_4287694.out        Processing sample: CBCH1_1w
mosdepth_highcovinds_4287695.out        Processing sample: CBCH2_2w
mosdepth_highcovinds_4287696.out        Processing sample: VDa_1w
mosdepth_highcovinds_4287697.out        Processing sample: CBCH3_1w
mosdepth_highcovinds_4287698.out        Processing sample: NAZa_1w
# These samples failed because they don't match the name pattern "_nodupl_wRG_clip.bam" 
# These samples are from Portinha et at., 2022 - the clip overlaps step was skipped due to low fraction of read overlap (~4%)



## run the failed samples

# 50.mosdepth_highcovinds.sh
# ------- START OF BASH SCRIPT ----------------------------------------------------------------------------------------

#!/bin/bash -l
#SBATCH -J mosdepth_highcovinds
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/logs/mosdepth_highcovinds_%j.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/logs/mosdepth_highcovinds_%j.err
#SBATCH -t 04:00:00
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --array=1-5
#SBATCH --ntasks 1
#SBATCH --mem-per-cpu=4GB
#SBATCH --mail-user ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# set the script to stop immediately if any command fails or any part of a pipeline fails, preventing silent errors and corrupted results
set -eo pipefail

# Load environment mosdepth
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/mosdepth

BAMDIR=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/bam_final/
FINALDIR=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp
NAMELIST=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/nodupl_highcovinds_rest.list

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp

# Get sample ID
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p $NAMELIST)

# Select BAM filename based on sample name
if [[ "$sample" =~ ^Ref_[0-9]+$ ]]; then
    bam="$BAMDIR/${sample}_nodupl_wRG_clip.bam"
else
    bam="$BAMDIR/${sample}_nodupl_wRG.bam"
fi

echo "Processing sample: $sample"

# Run mosdepth
mosdepth -t 1 -b 100000 -n \
    "$FINALDIR/${sample}_overlap_correction" \
    "$bam"

echo "Finished sample: ${sample}"

# ------- END OF BASH SCRIPT ----------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 4288844 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode

for file in mosdepth_highcovinds_*.out; do
    echo -e "$(basename "$file")\t$(head -n 1 "$file")"
done > job_sample_mapping_2.txt




### Get average coverage per each window ----------------

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp

#paste per region coverage data for all pixy samples
touch output.txt
touch ind_cov.list
for i in *regions.bed.gz; do
    zcat "$i" | awk -v OFS='\t' '{print $4}' | paste - output.txt > temp_output.txt
    mv temp_output.txt output.txt
    echo "Processed $i"
    echo $i >> ind_cov.list
done

#paste region info as the first column
zcat Ref_102_overlap_correction.regions.bed.gz | awk -v OFS='\t' '{print $1, $2, $3}' > 1stcol.txt #extract region info
paste 1stcol.txt output.txt > mosdepth_highcovinds_100kb_regions_depth.txt
rm 1stcol.txt output.txt

#remove regions that are not used: mitchondria (mtDNA), Wolbachia (wFhyb*), and Spiroplasma (Spiroplasma*), as well as Scaffold00 and Scaffold03
grep -v -e mtDNA -e wFhyb -e Spiroplasma -e Scaffold00 -e Scaffold03 mosdepth_highcovinds_100kb_regions_depth.txt > mosdepth_highcovinds_100kb_regions_depth_coreregions.txt
# rm mosdepth_highcovinds_100kb_regions_depth.txt


sort ind_cov.list >  ind_cov_sorted.list
sort Fin_5ind_species_fixed.tab > Fin_5ind_species_fixed_sorted.tab

paste ind_cov_sorted.list <(cut -f2 Fin_5ind_species_fixed_sorted.tab) > ind_cov_species.list
# !!! ind_cov_species.list is not in the same order as samples in the mosdepth_highcovinds_100kb_regions_depth_coreregions.txt !!!
# !!! correct sample order is at ind_cov.list


##### NOTES FROM INA ######
#make the ind_cov.list to have species vcf_id as 1st col and species as 2nd col
#modify it; resulting in ind_cov_species.list that has 1st row: chr, 2nd row: start, 3rd row: end, followed by ind/species info


### PROBABLY NOT USE THIS
# Maybe change the names to short versions first
{
    printf "chr\tstart\tend"
    cut -f1 ind_cov_species.list | while read sample; do
        printf "\t%s" "$sample"
    done
    printf "\n"
    cat mosdepth_highcovinds_100kb_regions_depth_coreregions.txt
} > mosdepth_highcovinds_100kb_regions_depth_named.txt



#get the files  mosdepth_highcovinds_100kb_regions_depth_coreregions.txt, ind_cov_species.list
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/mosdepth_highcovinds_100kb_regions_depth_coreregions.txt .
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/ind_cov_species.list .
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/ind_cov.list .

# everything later moved to /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/Fin_5ind_species



# _________________________________________________________________________________________________________________________________________________________________________________________________
# Exact same samples as Ina used


cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/groupfiles
cp Fin_5ind_species_Ina_fixed.tab /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/Fin_5ind_species_Ina/Fin_5ind_species_Ina_fixed.tab


# generating the nodupl_highcovinds.list of .bam files used in pixy
## final bam files have "Ref" prefix for references + some Finnish references have species suffix - the following code deals with that and generates the list for the following script.

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/Fin_5ind_species_Ina/

awk '{
    id = $1
    if (id ~ /^[0-9]+-/) {
        sub(/-.*/, "", id)
    }
    print "Ref_" id
}' Fin_5ind_species_Ina_fixed.tab > nodupl_highcovinds_Ina.list



# 52.mosdepth_highcovinds.sh
# ------- START OF BASH SCRIPT ----------------------------------------------------------------------------------------

#!/bin/bash -l
#SBATCH -J mosdepth_highcovinds
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/logs/mosdepth_highcovinds_Ina_%j.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/logs/mosdepth_highcovinds_Ina_%j.err
#SBATCH -t 01:00:00
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --array=1-25
#SBATCH --ntasks 1
#SBATCH --mem-per-cpu=4GB
#SBATCH --mail-user ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# set the script to stop immediately if any command fails or any part of a pipeline fails, preventing silent errors and corrupted results
set -eo pipefail

# Load environment mosdepth
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/mosdepth

BAMDIR=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/bam_final/
FINALDIR=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/Fin_5ind_species_Ina/
NAMELIST=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/Fin_5ind_species_Ina/nodupl_highcovinds_Ina.list

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/Fin_5ind_species_Ina

# Get sample ID
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p $NAMELIST)

# Select BAM filename based on sample name
if [[ "$sample" =~ ^Ref_[0-9]+$ ]]; then
    bam="$BAMDIR/${sample}_nodupl_wRG_clip.bam"
else
    bam="$BAMDIR/${sample}_nodupl_wRG.bam"
fi

echo "Processing sample: $sample"

# Run mosdepth
mosdepth -t 1 -b 100000 -n \
    "$FINALDIR/${sample}_overlap_correction" \
    "$bam"

echo "Finished sample: ${sample}"

# ------- END OF BASH SCRIPT ----------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 4297480 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode


### Get average coverage per each window ----------------

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/Fin_5ind_species_Ina

#paste per region coverage data for all pixy samples
touch output.txt
touch ind_cov.list
for i in *regions.bed.gz; do
    zcat "$i" | awk -v OFS='\t' '{print $4}' | paste - output.txt > temp_output.txt
    mv temp_output.txt output.txt
    echo "Processed $i"
    echo $i >> ind_cov.list
done

#paste region info as the first column
zcat Ref_102_overlap_correction.regions.bed.gz | awk -v OFS='\t' '{print $1, $2, $3}' > 1stcol.txt #extract region info
paste 1stcol.txt output.txt > mosdepth_highcovinds_100kb_regions_depth.txt
rm 1stcol.txt output.txt

#remove regions that are not used: mitchondria (mtDNA), Wolbachia (wFhyb*), and Spiroplasma (Spiroplasma*), as well as Scaffold00 and Scaffold03
grep -v -e mtDNA -e wFhyb -e Spiroplasma -e Scaffold00 -e Scaffold03 mosdepth_highcovinds_100kb_regions_depth.txt > mosdepth_highcovinds_100kb_regions_depth_coreregions.txt
# rm mosdepth_highcovinds_100kb_regions_depth.txt


sort ind_cov.list >  ind_cov_sorted.list
sort Fin_5ind_species_Ina_fixed.tab > Fin_5ind_species_Ina_fixed_sorted.tab

paste ind_cov_sorted.list <(cut -f2 Fin_5ind_species_Ina_fixed_sorted.tab) > ind_cov_species.list
# !!! ind_cov_species.list is not in the same order as samples in the mosdepth_highcovinds_100kb_regions_depth_coreregions.txt !!!
# !!! correct sample order is at ind_cov.list

scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/Fin_5ind_species_Ina/mosdepth_highcovinds_100kb_regions_depth_coreregions.txt .
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/Fin_5ind_species_Ina/ind_cov_species.list .
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/02.BAM/cov_highcovInds_5perSp/Fin_5ind_species_Ina/ind_cov.list .






