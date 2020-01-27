# Freebayes-GATK SNPs detection

This pipeline uses short sequencing reads to detect genomic variants (insertions, deletions and SNPs) in genome. Here, we combine the output of two known variant callers, freebayes and GATK, to generate a conservative list of variants. This specific example evaluates variants in two strains (A & B) of Neuropora crassa. 

## Create conda environment to run the programs needed

```
conda create neurospora
conda activate neurospora
conda install -c bioconda samtools freebayes vcftools bowtie2 bcftools bedops bedtools
```

## Reads quality filtering and mapping

1. Quality filtering using trimmomatic (100bp paired-end reads). 

```
java -jar trimmomatic-0.39.jar PE 712A_1.fastq.gz 712A_2.fastq.gz 712A_1_paired.fastq.gz 712A_1_unpaired.fastq.gz 712A_2_paired.fastq.gz 712A_2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:10:30 MINLEN:50 AVGQUAL:25
```

2. Index [genome of Neurospora crassa](https://www.ncbi.nlm.nih.gov/genome/?term=txid5141[orgn]).

```bowtie2-build GCF_000182925.2_NC12_genomic.fna Ncrassa```

3. Map filtered reads to N. crassa.

```
bowtie2  -x Ncrassa -1 712A_1_paired.fastq.gz -2 712A_2_paired.fastq.gz -S 712A_N.crassa.sam

bowtie2  -x Ncrassa -1 712B_1_paired.fastq.gz -2 712B_2_paired.fastq.gz -S 712B_N.crassa.sam
```

5. Convert sam file to bam, fix mate and sort.

```
for file in *.sam ; do base=${file##*/}; samtools view -S -b  $file > ${base%.*}.bam; samtools sort -n -o ${base%.*}.sorted.bam ${base%.*}.bam; samtools fixmate -m ${base%.*}.sorted.bam ${base%.*}.fixmate.bam;
samtools sort ${base%.*}.fixmate.bam -o ${base%.*}.fixmate.sorted.bam; samtools markdup -r ${base%.*}.fixmate.sorted.bam ${base%.*}.fixmate.sorted.dedup.bam; samtools sort -o ${base%.*}.sorted.bam ${base%.*}.fixmate.sorted.dedup.bam ;done 
```

## Freebayes

1. Run [freebayes](https://github.com/ekg/freebayes) for genomic variants detection.

```
freebayes -f GCF_000182925.2_NC12_genomic.fna -p 1 712A_N.crassa.fixmate.sorted.dedup.bam > 712A_N.crassa.vcf
freebayes -f GCF_000182925.2_NC12_genomic.fna -p 1 712B_N.crassa.fixmate.sorted.dedup.bam > 712B_N.crassa.vcf
```

2. Select variants unique to each strain.

```
bgzip -c 712B_N.crassa.vcf > 712B_N.crassa.vcf.gz
tabix -f -p vcf 712B_N.crassa.vcf.gz

bgzip -c 712A_N.crassa.vcf > 712A_N.crassa.vcf.gz
tabix -f -p vcf 712A_N.crassa.vcf.gz

perl /bin/vcftools/src/perl/vcf-isec -c 712A_N.crassa.vcf.gz 712B_N.crassa.vcf.gz > 712A_N.unique.vcf

perl /bin/vcftools/src/perl/vcf-isec -c 712B_N.crassa.vcf.gz 712A_N.crassa.vcf.gz > 712B_N.unique.vcf
```

3. Filter variants by coverage (> 10 mapping reads) and quality (> 30).

```
bcftools view -i 'FORMAT/DP>10' 712B_N.unique.vcf > 712B_N.unique_filtDP.vcf 
bcftools view -i 'FORMAT/DP>10' 712A_N.unique.vcf > 712A_N.unique_filtDP.vcf 

vcftools --vcf 712A_N.unique_filtDP.vcf  --recode --recode-INFO-all --minQ 30 --out 712A_N.unique_filtered
vcftools --vcf 712B_N.unique_filtDP.vcf  --recode --recode-INFO-all --minQ 30 --out 712B_N.unique_filtered
```

## GATK

1. Index genome file to run [GATK](https://gatk.broadinstitute.org/hc/en-us) for genomic variants detection.

```
gatk CreateSequenceDictionary -R GCF_000182925.2_NC12_genomic.fna

samtools faidx GCF_000182925.2_NC12_genomic.fna
```

2. Modify bam files. 

```
gatk AddOrReplaceReadGroups -I 712A_N.crassa.fixmate.sorted.dedup.bam -O 712A_N.crassa.fixmate.sorted.dedup.RG.bam --RGID 4  --RGLB lib1  --RGPL illumina  --RGPU unit1 --RGSM 20

gatk AddOrReplaceReadGroups -I 712B_N.crassa.fixmate.sorted.dedup.bam -O 712B_N.crassa.fixmate.sorted.dedup.RG.bam --RGID 4  --RGLB lib1  --RGPL illumina  --RGPU unit1 --RGSM 20

samtools index 712A_N.crassa.fixmate.sorted.dedup.RG.bam

samtools index 712B_N.crassa.fixmate.sorted.dedup.RG.bam
```

3. Run GATK. 

```
gatk HaplotypeCaller -R GCF_000182925.2_NC12_genomic.fna -I 712A_N.crassa.fixmate.sorted.dedup.RG.bam --sample-ploidy 1 -O 712A_N.crassa.vcf

gatk HaplotypeCaller -R GCF_000182925.2_NC12_genomic.fna -I 712B_N.crassa.fixmate.sorted.dedup.RG.bam --sample-ploidy 1 -O 712B_N.crassa.vcf
```

4. Select variants unique to each strain.

```
bgzip -c 712B_N.crassa.vcf > 712B_N.crassa.vcf.gz
tabix -f -p vcf 712B_N.crassa.vcf.gz

bgzip -c 712A_N.crassa.vcf > 712A_N.crassa.vcf.gz
tabix -f -p vcf 712A_N.crassa.vcf.gz

perl /bin/vcftools/src/perl/vcf-isec -c 712A_N.crassa.vcf.gz 712B_N.crassa.vcf.gz > 712A_N.unique.vcf

perl /bin/vcftools/src/perl/vcf-isec -c 712B_N.crassa.vcf.gz  712A_N.crassa.vcf.gz > 712B_N.unique.vcf
```

5. Filter variants by coverage (> 10 mapping reads) and quality (> 30).

```
bcftools view -i 'FORMAT/DP>10' 712B_N.unique.vcf > 712B_N.unique_filtDP.vcf 
bcftools view -i 'FORMAT/DP>10' 712A_N.unique.vcf > 712A_N.unique_filtDP.vcf 

vcftools --vcf 712A_N.unique_filtDP.vcf  --recode --recode-INFO-all --minQ 30 --out 712A_N.unique_filtered
vcftools --vcf 712B_N.unique_filtDP.vcf  --recode --recode-INFO-all --minQ 30 --out 712B_N.unique_filtered
```

## Merging freebayes and GATK results

1. Prepare files for merging.

```
bgzip -c 712A_N.unique_filtered.recode.vcf > 712A_N.unique_filtered.recode.vcf.gz
tabix -f -p vcf 712A_N.unique_filtered.recode.vcf.gz

bgzip -c 712B_N.unique_filtered.recode.vcf > 712B_N.unique_filtered.recode.vcf.gz
tabix -f -p vcf 712B_N.unique_filtered.recode.vcf.gz
```

2. Merge results.

```
perl /bin/vcftools/src/perl/vcf-isec /freebayes/712A_N.unique_filtered.recode.vcf.gz /gatk/712A_N.unique_filtered.recode.vcf.gz --force > 712A_N.unique_filt_merged.vcf

perl /bin/vcftools/src/perl/vcf-isec /freebayes/712B_N.unique_filtered.recode.vcf.gz /gatk/712B_N.unique_filtered.recode.vcf.gz --force > 712B_N.unique_filt_merged.vcf
```

## Intersecting results with gene information

1. Create bed file of gene position in genome.

```
convert2bed -i gff < GCF_000182925.2_NC12_genomic.gff > GCF_000182925.2_NC12_genomic.bed

awk '$8=="gene"{print}' GCF_000182925.2_NC12_genomic.bed > GCF_000182925.2_NC12_genes.bed
```

2. Intersect variants results with gene information.

```
bedtools intersect -a 712A_N.unique_filt_merged.vcf -b GCF_000182925.2_NC12_genes.bed -wb > 712A_N.unique_filt_merged_genes.bed

bedtools intersect -a 712B_N.unique_filt_merged.vcf -b GCF_000182925.2_NC12_genes.bed -wb > 712B_N.unique_filt_merged_genes.bed
```
