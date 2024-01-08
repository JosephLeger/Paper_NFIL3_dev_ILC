# NFIL3_dev_ILC
## Description
This repository contains codes and resources used to perform post-processing DNase-seq and Bulk RNA-seq data analysis.  
Before applying these scripts, raw data were processed using **Bulk_RNA** and **CUT_DNA** custom pipelines.  

## Commands used for pre-processing
```bash
### BULK RNA-SEQ
# Quality Check of raw FASTQ files
sh QC.sh Raw
sh MultiQC.sh QC/Raw

# Trimming to remove adapters
sh Trim.sh -S 4:15 -L 3 -T 3 -M 36 -I ./Ref/NexteraPE-PE_Clontech-TTT.fa:2:30:10 PE Raw

# Quality Check after trimming
sh QC.sh Trimmed/Trimmomatic
sh MultiQC.sh QC/Trimmed/Trimmomatic

# RSEM alignment
sh RSEM.sh PE Trimmed/Trimmomatic/Paired ./Ref/refdata-RSEM-mm39.108/mm39_108


### DNASE-SEQ
# Quality Check of raw FASTQ files
sh QC.sh Raw
sh MultiQC.sh QC/Raw

# Trimming to remove adapters and duplicated reads
sh Trim.sh -U 'Both' -S 4:15 -L 5 -T 5 -M 36 -I ./Ref/TruSeq3-SE_NexteraPE-PE.fa:2:30:7 -D True SE Raw

# Quality Check after trimming
sh QC.sh Trimmed/Trimmomatic
sh MultiQC.sh QC/Trimmed/Trimmomatic

# Mapping
sh Bowtie2.sh SE Trimmed/Trimmomatic ./Ref/refdata-Bowtie2-mm39/mm39
# Filtering
sh BowtieCheck.sh -N _sorted -T 10 -R false Mapped/mm39/BAM 

# Merge BAM files and BW generation
sh mergeBAM.sh -M DNAse-seq -N _Clum_Trimmed_sorted_filtered -R true Mapped/mm39/BAM SRA_list_DNAse-seq.csv
sh Bam2BW.sh -N _merged -F bigwig -M RPKM -R true Mapped/mm39/BAM

# Peak calling with MACS2
sh PeakyFinders.sh -U 'MACS2' -N _merged Mapped/mm39/BAM
sh PeakyFinders.sh -U 'MACS2' -N _Clum_Trimmed_sorted_filtered Mapped/mm39/BAM

```


