salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 01:00:00 --mem=2G --cpus-per-task=4
--------- INTERACTIVE JOB -----------------------------------------------------------------------------------

mamba activate admixure.env

VCF=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/filt/DP8.AN10.noScaff0003.mac2.thin20kb.vcf.gz

# Convert the VCF into PLINK format for ADMIXTURE analysis and keep only the desired individuals (in "LIST")
grep -v -F -f remove_samples.txt ind.list > /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/sample_list.txt
awk '{print $1, $1}' /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/sample_list.txt > sample_list.plink.txt

LIST=sample_list.plink.txt
OUTPATH=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture
OUTFILE=DP8.AN10.noScaff0003.mac2.thin20kb

plink --threads 4 --vcf $VCF --make-bed \
--double-id --allow-extra-chr --set-missing-var-ids @:# \
--keep $OUTPATH/$LIST \
--out $OUTPATH/$OUTFILE

# ADMIXTURE doesn't accept non-human chromosome names. Replace the first column by 0.

cd $OUTPATH
awk '{$1=0;print $0}' $OUTFILE.bim > $OUTFILE.bim.tmp
mv $OUTFILE.bim.tmp $OUTFILE.bim

------------- END OF INTERACTIVE SESSION -------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 3590242 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State

# 17_admixture.sh
------------- START OF THE SCRIPT -----------------------------------------------------------------------------------------------------------

#!/bin/bash -l
#SBATCH -J admixture
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/logs/admixture.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/logs/admixture.err
#SBATCH -t 00:30:00
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks 1
#SBATCH --mem=8G
#SBATCH --mail-user=ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# set the script to stop immediately if any command fails or any part of a pipeline fails, preventing silent errors and corrupted results
set -eo pipefail

# Load environment
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/admixure.env

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture

OUTPATH=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture
FILE=DP8.AN10.noScaff0003.mac2.thin20kb

for i in {3..14}
do
admixture --cv $OUTPATH/$FILE.bed $i > log${i}.out

done

grep -h CV log*.out > DP8.AN10.noScaff0003.mac2.thin20kb.cv.error

------------ END OF THE SCRIPT -----------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 3590269 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State

17_admixture_updated.sh
sacct -M biohpc_gen -j 3590510 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State


#Download the .dist file to local
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/DP8.AN10.noScaff0003.mac2.thin20kb.cv.error .
# .Q for ADMFILE in R script; admixture proportion for K populations
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/DP8.AN10.noScaff0003.mac2.thin20kb.6.Q .
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/DP8.AN10.noScaff0003.mac2.thin20kb.9.Q .
# .fam is from plink and it is important for list of samples in admixture
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/DP8.AN10.noScaff0003.mac2.thin20kb.fam .

------------ R SCRIPT -----------------------------------------------------------------------------------------------------------
library(ggplot2)

# Upload error file

error <- read.table("DP8.AN10.noScaff0003.mac2.thin20kb.cv.error")

error$V3 <- as.numeric(gsub(".*K=([0-9]+).*", "\\1", error$V3))

# plot all Ks
ggplot(error, aes(x = V3, y = V4)) +
  geom_point() +
  geom_line() +
  labs(x = "K", y = "CV error")

# plot k5-10
subset_data <- error[error$V3 >= 5 & error$V3 <= 10, ]

ggplot(subset_data, aes(x = V3, y = V4)) +
  geom_point() +
  geom_line() +
  labs(x = "K", y = "CV error")

### The best fit is K=6

# upload *.fam 
fam <- read.csv("DP8.AN10.noScaff0003.mac2.thin20kb.fam", sep=" ", header=F)

# upload *.Q 
tbl <- read.table("DP8.AN10.noScaff0003.mac2.thin20kb.6.Q")


# sample names
samples <- fam$V2
samples <- sub("\\.merged$", "", samples)

# margins for labels
par(mar = c(8, 4, 4, 2))  

# plot with sample names

barplot(t(as.matrix(tbl)), 
        col = rainbow(6),        # adjust to number of K
        xlab = NA,, 
        ylab = "Ancestry", 
        border = NA,  
        names.arg = samples,     
        las = 2,                 
        cex.names = 0.5)        

# ploting only LMUF samples

lmuf_idx <- grep("LMUF", samples)  
samples_lmuf <- samples[lmuf_idx]
tbl_lmuf <- tbl[lmuf_idx, ] 

barplot(t(as.matrix(tbl_lmuf)), 
        col = rainbow(6),        
        xlab = NA, 
        ylab = "Ancestry", 
        border = NA,  
        names.arg = samples_lmuf,     
        las = 2,                
        cex.names = 0.6,         
        space = 0.05)

------------ END OF R SCRIPT -----------------------------------------------------------------------------------------------------------
### OPTIONAL CONTINUATION

# remove F. exsecta from the dataset as well as unidentified references RN415-RN426

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/filt
# ind.list contains all samples (created in previous steps); 193 samples
# remove_samples.txt cantains samples which were removed for too high missingness; 7 samples

cp ind.list ind.list.tmp
sort ind.list.tmp > ind.list.sorted.tmp
# sort remove_samples.txt > remove_samples_sorted.txt             ## already done in 11.NeigborNet
grep -v -F -f remove_samples_sorted.txt ind.list.sorted.tmp > ind.list.filered.sorted.tmp  #186 samples; remove RNxxx manually - 176 samples

touch remove_exsecta_RNxxx.txt
vi remove_exsecta_RNxxx.txt    # paste names of the samples
grep -v -F -f remove_exsecta_RNxxx.txt ind.list.filered.sorted.tmp > ind.list.filered.sorted.noExsecta.txt # 172 samples

rm ind.list.tmp ind.list.sorted.tmp ind.list.filered.sorted.tmp 
mv ind.list.filered.sorted.noExsecta.txt ind.list.filtered.sorted.noExsecta.txt  #correcting a typo 



salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 01:00:00 --mem=2G --cpus-per-task=4
--------- INTERACTIVE JOB -----------------------------------------------------------------------------------

mamba activate admixure.env

VCF=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/filt/DP8.AN10.noScaff0003.mac2.thin20kb.vcf.gz

# Convert the VCF into PLINK format for ADMIXTURE analysis and keep only the desired individuals (in "LIST")
cp ind.list.filtered.sorted.noExsecta.txt /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admitxure/ind.list.filtered.sorted.noExsecta.txt
awk '{print $1, $1}' /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixure/ind.list.filtered.sorted.noExsecta.txt > ind.list.filtered.sorted.noExsecta.plink.txt

LIST=ind.list.filtered.sorted.noExsecta.plink.txt
OUTPATH=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture
OUTFILE=DP8.AN10.noScaff0003.mac2.thin20kb.noExsecta

plink --threads 4 --vcf $VCF --make-bed \
--double-id --allow-extra-chr --set-missing-var-ids @:# \
--keep $OUTPATH/$LIST \
--out $OUTPATH/$OUTFILE

# ADMIXTURE doesn't accept non-human chromosome names. Replace the first column by 0.

cd $OUTPATH
awk '{$1=0;print $0}' $OUTFILE.bim > $OUTFILE.bim.tmp
mv $OUTFILE.bim.tmp $OUTFILE.bim

------------- END OF INTERACTIVE SESSION -------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 3590242 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State

# 19_admixture_noExsecta.sh
------------- START OF THE SCRIPT -----------------------------------------------------------------------------------------------------------

#!/bin/bash -l
#SBATCH -J admixture_noExsecta
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/logs/admixture_noExsecta.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/logs/admixture_noExsecta.err
#SBATCH -t 00:30:00
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks 1
#SBATCH --mem=8G
#SBATCH --mail-user=ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# set the script to stop immediately if any command fails or any part of a pipeline fails, preventing silent errors and corrupted results
set -eo pipefail

# Load environment
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/admixure.env

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture

OUTPATH=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture
FILE=DP8.AN10.noScaff0003.mac2.thin20kb.noExsecta

for i in {3..14}
do
admixture --cv $OUTPATH/$FILE.bed $i > log${i}.out

done

grep -h CV log*.out > DP8.AN10.noScaff0003.mac2.thin20kb.noExsecta.cv.error

------------ END OF THE SCRIPT -----------------------------------------------------------------------------------------------------------

#Download the .dist file to local
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/DP8.AN10.noScaff0003.mac2.thin20kb.noExsecta.cv.error .
# .Q for ADMFILE in R script; admixture proportion for K populations
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/DP8.AN10.noScaff0003.mac2.thin20kb.noExsecta.5.Q .
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/DP8.AN10.noScaff0003.mac2.thin20kb.noExsecta.6.Q .
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/DP8.AN10.noScaff0003.mac2.thin20kb.noExsecta.7.Q .
# .fam is from plink and it is important for list of samples in admixture
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/DP8.AN10.noScaff0003.mac2.thin20kb.noExsecta.fam .

# download all Q files
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/DP8.AN10.noScaff0003.mac2.thin20kb.noExsecta.*.Q .


## Re-running K 2-14 (same script, just changed the loop "for i in {3..14} to for i in {2..14}")
sacct -M biohpc_gen -j 3617146 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State

#Download the .dist file to local
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/DP8.AN10.noScaff0003.mac2.thin20kb.noExsecta.cv.error .
# download all Q files
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/DP8.AN10.noScaff0003.mac2.thin20kb.noExsecta.*.Q .



# Dataset 3 - no references, only LMUF
# Dataset 4 - only LMUF F. aquilonia, F. polyctena and F. rufa
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Dataset 3
# keep only LMUF samples

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture
# list with no F. exsecta and no RNxxx "ind.list.filtered.sorted.noExsecta.txt"
grep '^LMUF' ind.list.filtered.sorted.noExsecta.txt > ind.list.onlyLMUF.txt
mkdir dataset3

salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 00:30:00 --mem=2G --cpus-per-task=4
--------- INTERACTIVE JOB -----------------------------------------------------------------------------------

mamba activate admixure.env

VCF=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/filt/DP8.AN10.noScaff0003.mac2.thin20kb.vcf.gz

# Convert the VCF into PLINK format for ADMIXTURE analysis and keep only the desired individuals (in "LIST")
awk '{print $1, $1}' /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/ind.list.onlyLMUF.txt > ind.list.onlyLMUF.plink.txt
mv ind.list.onlyLMUF.plink.txt dataset3/

LIST=ind.list.onlyLMUF.plink.txt
OUTPATH=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset3
OUTFILE=DP8.AN10.noScaff0003.mac2.thin20kb.onlyLMUF


plink --threads 4 --vcf $VCF --make-bed \
--double-id --allow-extra-chr --set-missing-var-ids @:# \
--keep $OUTPATH/$LIST \
--out $OUTPATH/$OUTFILE

# ADMIXTURE doesn't accept non-human chromosome names. Replace the first column by 0.

cd $OUTPATH
awk '{$1=0;print $0}' $OUTFILE.bim > $OUTFILE.bim.tmp
mv $OUTFILE.bim.tmp $OUTFILE.bim

------------- END OF INTERACTIVE SESSION -------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 3617744 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State

# 19_admixture_onlyLMUF.sh
------------- START OF THE SCRIPT -----------------------------------------------------------------------------------------------------------

#!/bin/bash -l
#SBATCH -J admixture_onlyLMUF
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/logs/admixture_onlyLMUF.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/logs/admixture_onlyLMUF.err
#SBATCH -t 00:30:00
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks 1
#SBATCH --mem=8G
#SBATCH --mail-user=ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# set the script to stop immediately if any command fails or any part of a pipeline fails, preventing silent errors and corrupted results
set -eo pipefail

# Load environment
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/admixure.env

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset3

OUTPATH=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset3
FILE=DP8.AN10.noScaff0003.mac2.thin20kb.onlyLMUF

for i in {2..10}
do
admixture --cv $OUTPATH/$FILE.bed $i > log${i}.out

done

grep -h CV log*.out > DP8.AN10.noScaff0003.mac2.thin20kb.onlyLMUF.cv.error

------------ END OF THE SCRIPT -----------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 3617846 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State

#Download the cv.error file to local
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset3/DP8.AN10.noScaff0003.mac2.thin20kb.onlyLMUF.cv.error .
# .Q for ADMFILEs in R script; admixture proportion for K populations
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset3/DP8.AN10.noScaff0003.mac2.thin20kb.onlyLMUF.*.Q .
# .fam is from plink and it is important for list of samples in admixture
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset3/DP8.AN10.noScaff0003.mac2.thin20kb.onlyLMUF.fam .


------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Dataset 4
# keep only LMUF F. aquilonia, F. polyctena and F. rufa

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture


touch ind.list.noPrat.noLug.noEx.txt
vi ind.list.noPrat.noLug.noEx.tmp    # paste names of the samples
awk '{print $0".merged"}' ind.list.noPrat.noLug.noEx.tmp > ind.list.noPrat.noLug.noEx.txt
rm ind.list.noPrat.noLug.noEx.tmp

mkdir dataset4

salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 00:30:00 --mem=2G --cpus-per-task=4
--------- INTERACTIVE JOB -----------------------------------------------------------------------------------

mamba activate admixure.env

VCF=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/filt/DP8.AN10.noScaff0003.mac2.thin20kb.vcf.gz

# Convert the VCF into PLINK format for ADMIXTURE analysis and keep only the desired individuals (in "LIST")
awk '{print $1, $1}' /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/ind.list.noPrat.noLug.noEx.txt > ind.list.noPrat.noLug.noEx.plink.txt
mv ind.list.noPrat.noLug.noEx.plink.txt dataset4/

LIST=ind.list.noPrat.noLug.noEx.plink.txt
OUTPATH=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset4
OUTFILE=DP8.AN10.noScaff0003.mac2.thin20kb.onlyLMUF.noPrat.noLug.noEx


plink --threads 4 --vcf $VCF --make-bed \
--double-id --allow-extra-chr --set-missing-var-ids @:# \
--keep $OUTPATH/$LIST \
--out $OUTPATH/$OUTFILE

# ADMIXTURE doesn't accept non-human chromosome names. Replace the first column by 0.

cd $OUTPATH
awk '{$1=0;print $0}' $OUTFILE.bim > $OUTFILE.bim.tmp
mv $OUTFILE.bim.tmp $OUTFILE.bim

------------- END OF INTERACTIVE SESSION -------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 3619140 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State

# 19_admixture_onlyLMUF.noPrat.noLug.noEx.sh
------------- START OF THE SCRIPT -----------------------------------------------------------------------------------------------------------

#!/bin/bash -l
#SBATCH -J admixture_onlyLMUF.noPrat.noLug.noEx
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/logs/admixture_onlyLMUF.noPrat.noLug.noEx.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/logs/admixture_onlyLMUF.noPrat.noLug.noEx.err
#SBATCH -t 00:30:00
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks 1
#SBATCH --mem=8G
#SBATCH --mail-user=ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# set the script to stop immediately if any command fails or any part of a pipeline fails, preventing silent errors and corrupted results
set -eo pipefail

# Load environment
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/admixure.env

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset4

OUTPATH=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset4
FILE=DP8.AN10.noScaff0003.mac2.thin20kb.onlyLMUF.noPrat.noLug.noEx

for i in {2..10}
do
admixture --cv $OUTPATH/$FILE.bed $i > log${i}.out

done

grep -h CV log*.out > DP8.AN10.noScaff0003.mac2.thin20kb.onlyLMUF.noPrat.noLug.noEx.cv.error

------------ END OF THE SCRIPT -----------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 3619142 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State

#Download the cv.error file to local
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset4/DP8.AN10.noScaff0003.mac2.thin20kb.onlyLMUF.noPrat.noLug.noEx.cv.error .
# .Q for ADMFILEs in R script; admixture proportion for K populations
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset4/DP8.AN10.noScaff0003.mac2.thin20kb.onlyLMUF.noPrat.noLug.noEx.*.Q .
# .fam is from plink and it is important for list of samples in admixture
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset4/DP8.AN10.noScaff0003.mac2.thin20kb.onlyLMUF.noPrat.noLug.noEx.fam .


------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Dataset 5
# All F. aquilonia, F. polyctena and F. rufa

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture

touch ind.list.Faq.Fru.Fpol.tmp
vi ind.list.Faq.Fru.Fpol.tmp    # paste names of the LMUF samples
awk '{print $0".merged"}' ind.list.Faq.Fru.Fpol.tmp > ind.list.Faq.Fru.Fpol.txt
vi ind.list.Faq.Fru.Fpol.txt # pastr names of the reference samples
sort ind.list.Faq.Fru.Fpol.txt > ind.list.Faq.Fru.Fpol.sorted.txt # 125 samples
rm ind.list.Faq.Fru.Fpol.tmp

mkdir dataset5

salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 00:30:00 --mem=2G --cpus-per-task=4
--------- INTERACTIVE JOB -----------------------------------------------------------------------------------

mamba activate admixure.env

VCF=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/filt/DP8.AN10.noScaff0003.mac2.thin20kb.vcf.gz

# Convert the VCF into PLINK format for ADMIXTURE analysis and keep only the desired individuals (in "LIST")
awk '{print $1, $1}' /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/ind.list.Faq.Fru.Fpol.sorted.txt > ind.list.Faq.Fru.Fpol.sorted.plink.txt
mv ind.list.Faq.Fru.Fpol.sorted.plink.txt dataset5/

LIST=ind.list.Faq.Fru.Fpol.sorted.plink.txt
OUTPATH=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset5
OUTFILE=DP8.AN10.noScaff0003.mac2.thin20kb.Faq.Fru.Fpol


plink --threads 4 --vcf $VCF --make-bed \
--double-id --allow-extra-chr --set-missing-var-ids @:# \
--keep $OUTPATH/$LIST \
--out $OUTPATH/$OUTFILE

# ADMIXTURE doesn't accept non-human chromosome names. Replace the first column by 0.

cd $OUTPATH
awk '{$1=0;print $0}' $OUTFILE.bim > $OUTFILE.bim.tmp
mv $OUTFILE.bim.tmp $OUTFILE.bim

------------- END OF INTERACTIVE SESSION -------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 3619140 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State

# 19_admixture.Faq.Fru.Fpol.sh
------------- START OF THE SCRIPT -----------------------------------------------------------------------------------------------------------

#!/bin/bash -l
#SBATCH -J admixture.Faq.Fru.Fpol
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/logs/admixture.Faq.Fru.Fpol.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/logs/admixture.Faq.Fru.Fpol.err
#SBATCH -t 00:30:00
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks 1
#SBATCH --mem=8G
#SBATCH --mail-user=ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# set the script to stop immediately if any command fails or any part of a pipeline fails, preventing silent errors and corrupted results
set -eo pipefail

# Load environment
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/admixure.env

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset5

OUTPATH=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset5
FILE=DP8.AN10.noScaff0003.mac2.thin20kb.Faq.Fru.Fpol

for i in {2..10}
do
admixture --cv $OUTPATH/$FILE.bed $i > log${i}.out

done

grep -h CV log*.out > DP8.AN10.noScaff0003.mac2.thin20kb.Faq.Fru.Fpol.cv.error

------------ END OF THE SCRIPT -----------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 3619941 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State

#### EDIT THE FILE NAMES
#Download the cv.error file to local
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset5/DP8.AN10.noScaff0003.mac2.thin20kb.Faq.Fru.Fpol.cv.error .
# .Q for ADMFILEs in R script; admixture proportion for K populations
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset5/DP8.AN10.noScaff0003.mac2.thin20kb.Faq.Fru.Fpol.*.Q .
# .fam is from plink and it is important for list of samples in admixture
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset5/DP8.AN10.noScaff0003.mac2.thin20kb.Faq.Fru.Fpol.fam .


------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Dataset 6
# removing samples from Russia and Scotland + F. exsecta samples + undetermined ref. samples
# 115-Flug (Russia), Lai_1w, Lai_2w, Loa_1w (Scotland), RN415, RN416, RN417, RN418, RN419, RN420, RN421, RN422, RN423, RN424, RN425, RN426, s353, s354 (undetermined samples)
# LMUF_00011a, LMUF_00030a, LMUF_00034a, LMUF_00427b

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/filt
cp remove_samples_sorted.txt remove_samples_D6.txt # this file contains samples which were removed because of high missingness
vi remove_samples_D6.txt # add sample names which should be filtered out
sort remove_samples_D6.txt > remove_samples_D6_sorted.txt
sort ind.list > ind.list.sorted
grep -v -F -f remove_samples_D6_sorted.txt ind.list.sorted > ind_list_D6 # 168 lines
# SANITY CHECK: ind.list.sorted (193) - remove_samples_D6_sorted.txt (25) = ind_list_D6 (168) 



salloc --clusters=biohpc_gen --partition=biohpc_gen_inter -t 01:00:00 --mem=2G --cpus-per-task=4
--------- INTERACTIVE JOB -----------------------------------------------------------------------------------

mamba activate admixure.env

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/filt

# Convert the VCF into PLINK format for ADMIXTURE analysis and keep only the desired individuals (in "LIST")
cp ind_list_D6 /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset6/ind_list_D6.txt

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset6
awk '{print $1, $1}' ind_list_D6.txt > ind_list_D6.plink.txt

VCF=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/filt/DP8.AN10.noScaff0003.mac2.thin20kb.vcf.gz
LIST=ind_list_D6.plink.txt
OUTPATH=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset6
OUTFILE=DP8.AN10.noScaff0003.mac2.thin20kb.D6

plink --threads 4 --vcf $VCF --make-bed \
--double-id --allow-extra-chr --set-missing-var-ids @:# \
--keep $OUTPATH/$LIST \
--out $OUTPATH/$OUTFILE

# ADMIXTURE doesn't accept non-human chromosome names. Replace the first column by 0.

awk '{$1=0;print $0}' $OUTFILE.bim > $OUTFILE.bim.tmp
mv $OUTFILE.bim.tmp $OUTFILE.bim

------------- END OF INTERACTIVE SESSION -------------------------------------------------------------------------------------------------


# 25_admixture_D6.sh
------------- START OF THE SCRIPT -----------------------------------------------------------------------------------------------------------
#!/bin/bash -l
#SBATCH -J admixture_D6
#SBATCH -o /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/logs/admixture_D6.out
#SBATCH -e /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/logs/admixture_D6.err
#SBATCH -t 00:30:00
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks 1
#SBATCH --mem=8G
#SBATCH --mail-user=ada.crhonkova@seznam.cz
#SBATCH --mail-type=END,FAIL

# set the script to stop immediately if any command fails or any part of a pipeline fails, preventing silent errors and corrupted results
set -eo pipefail

# Load environment
eval "$(mamba shell hook --shell bash)"
mamba activate /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/mamba_envs/admixure.env

cd /dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset6

echo "Start"
OUTPATH=/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset6
FILE=DP8.AN10.noScaff0003.mac2.thin20kb.D6

for i in {2..10}
do
admixture --cv $OUTPATH/$FILE.bed $i > log${i}.out

done

grep -h CV log*.out > DP8.AN10.noScaff0003.mac2.thin20kb.D6.cv.error

echo "End"
------------ END OF THE SCRIPT -----------------------------------------------------------------------------------------------------------
sacct -M biohpc_gen -j 4013937 --format=JobID,Elapsed,MaxRSS,AllocCPUS,State

#### EDIT THE FILE NAMES
#Download the cv.error file to local
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset6/DP8.AN10.noScaff0003.mac2.thin20kb.D6.cv.error .
# .Q for ADMFILEs in R script; admixture proportion for K populations
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset6/DP8.AN10.noScaff0003.mac2.thin20kb.D6.*.Q .
# .fam is from plink and it is important for list of samples in admixture
scp -r re98maw@cool.hpc.lrz.de:/dss/dsslegfs01/pn73qe/pn73qe-dss-0002/Formica_WGS/WGS_2024_2025/03.VCF/admixture/dataset6/DP8.AN10.noScaff0003.mac2.thin20kb.D6.fam .


