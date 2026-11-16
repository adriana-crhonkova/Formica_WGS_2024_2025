#Adjusted from Ina Satokagas https://github.com/isatok/Wood_ant_local_genomic_patterns/blob/main/Pixy/1_allsites_vcf.sh

#### This script is to compute local Pi, Dxy and Fst estimates with PIXY ####

mkdir /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy
mkdir /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/groupfiles
mkdir /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/logs

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy

### 1. Prepare group files for PIXY ------

# DATASET 1 - Bavaria between species comparisons




### 2. Compute PIXY ------
# Bav_5inds_pixy.sh
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
--populations groupfiles/Bav_5inds_pixy.tab \
--window_size 100000 \
--output_prefix Bav_5inds_pixy \
--n_cores 8

echo "Done"

# ------- END OF BASH SCRIPT ----------------------------------------------------------------------------------------
sacct -M biohpc_gen -j XXX --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode
