#!/bin/bash 
#SBATCH -J fasta2phylip
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/fasta2phylip.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/logs/fasta2phylip.err
#SBATCH -t 00:10:00
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --mail-user ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# Load environment
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/mafft.env


# https://github.com/chinchliff/physcripts/blob/master/fasta2phylip.pl
perl Fasta2Phylip.pl /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/indiv_samples_mtDNA/LMUF_mtDNAsamples_aligned.fa /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/indiv_samples_mtDNA/LMUF_mtDNAsamples_aligned.phy

echo "Done"

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 4192027 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State,ExitCode



# use in HaploThread to build the network
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/indiv_samples_mtDNA/LMUF_mtDNAsamples_aligned_short_names.phy .
# use as a list of samples to filter the trait file
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/04.VCF.mtDNA/samples.list.LMUF .
