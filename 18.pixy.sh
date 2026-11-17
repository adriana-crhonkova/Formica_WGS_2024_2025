#Adjusted from Ina Satokagas https://github.com/isatok/Wood_ant_local_genomic_patterns/blob/main/Pixy/1_allsites_vcf.sh

#### This script is to compute local Pi, Dxy and Fst estimates with PIXY ####

mkdir /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy
mkdir /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/groupfiles
mkdir /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/logs

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy

### 1. Prepare group files for PIXY ------

# DATASET 1 - Bavaria between species comparisons
cd groupfiles
vi Bav_5ind_species.tab

# INSERT FOLLOWING
LMUF_00005d.merged      aq
LMUF_00018a.merged      aq
LMUF_00031a.merged      aq
LMUF_00036a.merged      aq
LMUF_00109b.merged      aq
LMUF_00132b.merged      pol
LMUF_00312a.merged      pol
LMUF_00324c.merged      pol
LMUF_00341c.merged      pol
LMUF_00400c.merged      pol
LMUF_00010a.merged      ru
LMUF_00025a.merged      ru
LMUF_00041a.merged      ru
LMUF_00101a.merged      ru
LMUF_00118c.merged      ru
LMUF_00012a.merged      lu
LMUF_00026a.merged      lu
LMUF_00035a.merged      lu
LMUF_00202a.merged      lu
LMUF_00237a.merged      lu

# Check whether the columns are really separated by tab
cat -A Bav_5ind_species.tab
## Expected fortmat LMUF_00005d.merged^Iaq$

# If some of the output is in format "LMUF_00025a.merged      ru$", run following:
awk 'BEGIN {OFS="\t"} {$1=$1; print}' \
Bav_5ind_species.tab \
> Bav_5ind_species_fixed.tab

# sanity check 
cat -A Bav_5ind_species_fixed.tab




### 2. Compute PIXY ------
# 43.Bav_5inds_pixy.sh
# ------- START OF BASH SCRIPT ----------------------------------------------------------------------------------------

#!/bin/bash -l
#SBATCH -J Bav_5inds_pixy
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/logs/Bav_5inds_pixy_%j.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/logs/Bav_5inds_pixy_%j.err
#SBATCH -t 24:00:00
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks 8
#SBATCH --mem=4G
#SBATCH --mail-user ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# set the script to stop immediately if any command fails or any part of a pipeline fails, preventing silent errors and corrupted results
set -eo pipefail

# Load environment
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/pixy.env

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy

echo "Starting pixy"

pixy --stats pi dxy fst \
--vcf /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/highcovsamples_filtered_noindels_rmdup.vcf.gz \
--populations groupfiles/Bav_5ind_species_fixed.tab \
--window_size 100000 \
--output_folder /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy \
--output_prefix Bav_5ind_species \
--n_cores 8

echo "Done"

# ------- END OF BASH SCRIPT ----------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 4254687 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode

# the job failed, but because I was missing "\" after --output_folder path and therefore flags under it were not part of the pixy flags.
# I renamed the output files with the --output_prefix Bav_5ind_species
## Run time  02:44:31 -> 5 hours should be enough
## Used memory 275664K = 275,664 MB -> 2GB should be enough
