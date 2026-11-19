#Adjusted from Ina Satokagas https://github.com/isatok/Wood_ant_local_genomic_patterns/blob/main/Pixy/1_allsites_vcf.sh

#### This script is to compute local Pi, Dxy and Fst estimates with PIXY ####

mkdir /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy
mkdir /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/groupfiles
mkdir /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/logs

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy

### 1. Prepare group files for PIXY ------

# DATASET 1 - Bavaria between species comparisons
## The samples were chosen based on the least missing data + each from different location

cd groupfiles
vi Bav_5ind_species.tab

# INSERT FOLLOWING
LMUF_00548a.merged	aq
LMUF_00226a.merged	aq
LMUF_00523a.merged	aq
LMUF_00510a.merged	aq
LMUF_00020a.merged	aq
LMUF_00035a.merged	lu
LMUF_00027a.merged	lu
LMUF_00237a.merged	lu
LMUF_00202a.merged	lu
LMUF_00012a.merged	lu
LMUF_00525a.merged	pol
LMUF_00400c.merged	pol
LMUF_00324c.merged	pol
LMUF_00132b.merged	pol
LMUF_00430c.merged	pol
LMUF_00025a.merged	ru
LMUF_00010a.merged	ru
LMUF_00041a.merged	ru
LMUF_00118c.merged	ru
LMUF_00407b.merged	ru

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
#SBATCH -t 5:00:00
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks 8
#SBATCH --mem=2G
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

# new run with correct samples (the least missing data), the wrong (randomly chosen) moved to run1_wrong_samples
sacct -M biohpc_gen -j 4258695 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode


scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/Bav_5ind_species_dxy.txt .
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/Bav_5ind_species_fst.txt .
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/Bav_5ind_species_pi.txt .


# _________________________________________________________________________________________________________________________________________________
## Comparisons between species from different locations (Finland X Bavaria)

# FINNISH SAMPLES (F. polyctena from Switzerland)
104-Faqu  aq
44-FaquH	aq
107-Faqu	aq
109-Faqu	aq
106-Faqu	aq
108-Flug	lu
124-Flug	lu
113-Flug	lu
11-Flug	lu
56-Flug	lu
CBCH1_1w	pol
CBCH2_2w	pol
VDa_1w	pol
CBCH3_1w	pol
NAZa_1w	pol
16-Frufa	ru
37-Frufa	ru
102-Frufa	ru
26-Frufa	ru
72-Frufa	ru




# DATASET 2 - Aquilonia FinXBav


cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/groupfiles

vi BavXFin_5ind_aquilonia.tab
# insert following
LMUF_00548a.merged	aq_Bav
LMUF_00226a.merged	aq_Bav
LMUF_00523a.merged	aq_Bav
LMUF_00510a.merged	aq_Bav
LMUF_00020a.merged	aq_Bav
104-Faqu  aq_Fin
44-FaquH	aq_Fin
107-Faqu	aq_Fin
109-Faqu	aq_Fin
106-Faqu	aq_Fin

# Check whether the columns are really separated by tab
cat -A BavXFin_5ind_aquilonia.tab
## Expected fortmat LMUF_00005d.merged^Iaq$

# If some of the output is in format "LMUF_00025a.merged      ru$", run following:
awk 'BEGIN {OFS="\t"} {$1=$1; print}' \
BavXFin_5ind_aquilonia.tab \
> BavXFin_5ind_aquilonia_fixed.tab

# sanity check 
cat -A BavXFin_5ind_aquilonia_fixed.tab


### 2. Compute PIXY ------
# 44.BavXFin_5ind_aquilonia_pixy.sh
# ------- START OF BASH SCRIPT ----------------------------------------------------------------------------------------

#!/bin/bash -l
#SBATCH -J BavXFin_5ind_aquilonia_pixy
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/logs/BavXFin_5ind_aquilonia_pixy_%j.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/logs/BavXFin_5ind_aquilonia_pixy_%j.err
#SBATCH -t 5:00:00
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks 8
#SBATCH --mem=2G
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
--populations groupfiles/BavXFin_5ind_aquilonia_fixed.tab \
--window_size 100000 \
--output_folder /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy \
--output_prefix BavXFin_5ind_aquilonia \
--n_cores 8

echo "Done"

# ------- END OF BASH SCRIPT ----------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 4259796 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode



# DATASET 3 - Rufa FinXBav

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/groupfiles

vi BavXFin_5ind_rufa.tab
# insert following
LMUF_00025a.merged	ru_bav
LMUF_00010a.merged	ru_bav
LMUF_00041a.merged	ru_bav
LMUF_00118c.merged	ru_bav
LMUF_00407b.merged	ru_bav
16-Frufa	ru_Fin
37-Frufa	ru_Fin
102-Frufa	ru_Fin
26-Frufa	ru_Fin
72-Frufa	ru_Fin

# Check whether the columns are really separated by tab
cat -A BavXFin_5ind_rufa.tab
## Expected fortmat LMUF_00005d.merged^Iaq$

### 2. Compute PIXY ------
# 45.BavXFin_5ind_rufa_pixy.sh
# ------- START OF BASH SCRIPT ----------------------------------------------------------------------------------------

#!/bin/bash -l
#SBATCH -J BavXFin_5ind_rufa_pixy
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/logs/BavXFin_5ind_rufa_pixy_%j.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy/logs/BavXFin_5ind_rufa_pixy_%j.err
#SBATCH -t 5:00:00
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks 8
#SBATCH --mem=2G
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
--populations groupfiles/BavXFin_5ind_rufa.tab \
--window_size 100000 \
--output_folder /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/05.gVCF/pixy \
--output_prefix BavXFin_5ind_rufa \
--n_cores 8

echo "Done"

# ------- END OF BASH SCRIPT ----------------------------------------------------------------------------------------
sacct -M biohpc_gen -j XXX --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode








