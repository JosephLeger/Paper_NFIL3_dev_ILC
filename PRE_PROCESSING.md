# DATA PRE-PROCESSING

## DESCRIPTION
Before applying provided R scripts, Bulk RNA-seq and DNase-seq data were pre-processed using custom pipelines **Bulk_RNA_alignment** and **CUT_DNA**.
To reproduce the same pre-processing steps, we provided here the exact command lines used for these respective custom pipelines, or the equivalent command line syntax while using directly tools.  
Mapping and annotation were realized using Ensembl [GRCm39](https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz) and [GTF release 108](https://ftp.ensembl.org/pub/release-108/gtf/mus_musculus/Mus_musculus.GRCm39.108.gtf.gz).

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
rsem-calculate-expression -p 8 --paired-end --star --star-gzipped-read-file $R1 $R2 ./Ref/refdata-RSEM-mm39.108/mm39_108 ${output}
```

## DNASE-SEQ
### Requirements
```
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
```

### Using custom pipeline
```bash
# Quality Check of raw FASTQ files
sh QC.sh Raw
sh MultiQC.sh QC/Raw

# Trimming to remove adapters and duplicated reads
sh Trim.sh -U 'Both' -S 4:15 -L 5 -T 5 -M 36 -I ./Ref/TruSeq3-SE_NexteraPE-PE.fa:2:30:7 -D True SE Raw

# Quality Check after trimming
sh QC.sh Trimmed/Trimmomatic
sh MultiQC.sh QC/Trimmed/Trimmomatic

# Mapping using Bowtie2
sh Bowtie2.sh SE Trimmed/Trimmomatic ./Ref/refdata-Bowtie2-mm39/mm39

# Filtering
sh BowtieCheck.sh -N _sorted -T 10 -R false Mapped/mm39/BAM 

# Merge BAM files and BW generation
sh mergeBAM.sh -M DNAse-seq -N _Clum_Trimmed_sorted_filtered -R true Mapped/mm39/BAM SRA_list_DNAse-seq.csv
sh Bam2BW.sh -N _merged -F bigwig -M RPKM -R true Mapped/mm39/BAM

# Peak calling with MACS2
sh PeakyFinders.sh -U 'MACS2' -N _merged Mapped/mm39/BAM
sh PeakyFinders.sh -U 'MACS2' -N _Clum_Trimmed_sorted_filtered Mapped/mm39/BAM

# Motif enrichment analysis and peak annotation with HOMER from patterns identified by R DiffBind
sh 8_Annotate.sh -N '[' -R 200 -L '8,10,12' -A true -M true NFIL3_dev_ILC/Peaks ./Ref/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa ./Ref/Mus_musculus.GRCm39.108.gtf
sh 8_Annotate.sh -N 'ALP_open' -R 200 -L '8,10,12' -A true -M true NFIL3_dev_ILC/Peaks ./Ref/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa ./Ref/Mus_musculus.GRCm39.108.gtf
```

### Using directly tools
```bash
# Trimming to remove duplicated reads and adapters
clumpify.sh in=${file} out=${output} dedupe=t subs=0
trimmomatic SE -threads 4 ${file} ${output} SLIDINGWINDOW:4:15 LEADING:5 TRAILING:5 MINLEN:36 ./Ref/TruSeq3-SE_NexteraPE-PE.fa:2:30:7

# Mapping using Bowtie2
bowtie2 -p 2 -N 0 -x ./Ref/refdata-Bowtie2-mm39/mm39 -U ${file} -S Mapped/mm39/SAM/${file}.sam
picard SortSam INPUT=Mapped/mm39/SAM/${file}.sam OUTPUT=Mapped/mm39/BAM/${file}_sorted.bam \
VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp SORT_ORDER=coordinate

# Filtering
samtools view -h ${file} | samtools view -b -Sq 10  > ${bam_filtered}
samtools index ${bam_filtered} ${bai_filtered}

# Peak calling with MACS2
macs2 callpeak -t ${file} -f BAM -g 1.87e9 --nomodel --shift 50 --extsize 100 -n ${output} --outdir ${outdir}

# Motif enrichment analysis and peak annotation with HOMER from patterns identified by R DiffBind
annotatePeaks.pl ${file} ./Ref/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa -gtf ./Ref/Mus_musculus.GRCm39.108.gtf > ${outdir}/${file}_annotated.txt
findMotifsGenome.pl ${file} ./Ref/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa ${current_tag} -size 200 -len '8,10,12' -S 40
```



