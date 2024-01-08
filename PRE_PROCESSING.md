# DATA PRE-PROCESSING

## DESCRIPTION
#===============================================================================
Before applying provided R scripts, Bulk RNA-seq and DNase-seq data were pre-processed using custom pipelines **Bulk_RNA_alignment** and **CUT_DNA**.
To reproduce the same pre-processing steps, we provided here the exact command lines used for these respective custom pipelines, or the equivalent command lines while using directly tools.

## BULK RNA-SEQ
### Requirements
```
Name                        Version
fastqc                      0.11.9
multiqc                     1.13
trimmomatic                 0.39
rsem                        1.3.2
star                        2.7.5a
```

### Using custom pipeline
```bash
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
```

### Using directly tools
```bash
# Trimming to remove adapters
trimmomatic PE -threads 4 $R1 $R2 ${outdir}/Paired/${P1} ${outdir}/Unpaired/${U1} ${outdir}/Paired/${P2} ${outdir}/Unpaired/${U2} \
SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36 ./Ref/NexteraPE-PE_Clontech-TTT.fa:2:30:10

# RSEM alignment 
rsem-calculate-expression -p 8 --paired-end --star --star-gzipped-read-file $R1 $R2 ./Ref/refdata-RSEM-mm39.108/mm39_108 RSEM/$output
```

## DNASE-SEQ
### Requirements
Name                        Version
fastqc                      0.11.9
multiqc                     1.13
trimmomatic                 0.39
bbmap                       39.00
bowtie2                     2.5.1
samtools                    1.15.1
picard                      2.23.5
bedtools                    2.30.0
ucsc-bedgraphtobigwig       377
gcc                         11.2.0
macs2                       2.2.7.1
homer                       4.11




```bash
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
