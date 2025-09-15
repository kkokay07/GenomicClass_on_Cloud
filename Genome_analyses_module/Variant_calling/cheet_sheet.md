# Simple Joint Calling Pipeline for Tetraploid Brassica - Server Execution

## Overview
This simplified pipeline (with 4 samples: it can be run for any number of samples) uses direct joint calling with CombineGVCFs for multiple Brassica samples (ploidy=4). All commands are formatted for server execution using `nohup`.

---

## Required Files Setup

```bash
# Create directory structure
mkdir -p brassica_analysis/{raw_data,trimmed,aligned,variants,logs}
cd brassica_analysis

# Sample list (modify with your actual sample names)
echo -e "sample105\nsample106\nsample107\nsample108" > sample_list.txt
```

---

## Step 1: Quality Control and Trimming

**Run for each sample with nohup:**

```bash
# Sample 105
nohup fastp \
  -i raw_data/sample105.R1.fastq.gz \
  -I raw_data/sample105.R2.fastq.gz \
  -o trimmed/sample105.R1.trim.fastq.gz \
  -O trimmed/sample105.R2.trim.fastq.gz \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --thread 8 \
  --html trimmed/sample105.fastp.html \
  --json trimmed/sample105.fastp.json \
  > logs/sample105.fastp.log 2>&1 &

# Sample 106
nohup fastp \
  -i raw_data/sample106.R1.fastq.gz \
  -I raw_data/sample106.R2.fastq.gz \
  -o trimmed/sample106.R1.trim.fastq.gz \
  -O trimmed/sample106.R2.trim.fastq.gz \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --thread 8 \
  --html trimmed/sample106.fastp.html \
  --json trimmed/sample106.fastp.json \
  > logs/sample106.fastp.log 2>&1 &

# Sample 107
nohup fastp \
  -i raw_data/sample107.R1.fastq.gz \
  -I raw_data/sample107.R2.fastq.gz \
  -o trimmed/sample107.R1.trim.fastq.gz \
  -O trimmed/sample107.R2.trim.fastq.gz \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --thread 8 \
  --html trimmed/sample107.fastp.html \
  --json trimmed/sample107.fastp.json \
  > logs/sample107.fastp.log 2>&1 &

# Sample 108
nohup fastp \
  -i raw_data/sample108.R1.fastq.gz \
  -I raw_data/sample108.R2.fastq.gz \
  -o trimmed/sample108.R1.trim.fastq.gz \
  -O trimmed/sample108.R2.trim.fastq.gz \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --thread 8 \
  --html trimmed/sample108.fastp.html \
  --json trimmed/sample108.fastp.json \
  > logs/sample108.fastp.log 2>&1 &

# Wait for all jobs to complete
wait
echo "Quality control completed for all samples"
```

---

## Step 2: Reference Genome Preparation

```bash
# Download and prepare reference (run once)
nohup bash -c "
datasets download genome accession GCA_015484525.1 --include genome
unzip ncbi_dataset.zip
mv ncbi_dataset/data/GCA_015484525.1/GCA_015484525.1_UDSC_Var_1.1_genomic.fna reference.fasta

# Index reference
bwa index reference.fasta
samtools faidx reference.fasta
gatk CreateSequenceDictionary -R reference.fasta -O reference.dict
" > logs/reference_prep.log 2>&1 &

wait
echo "Reference genome preparation completed"
```

---

## Step 3: Alignment and BAM Processing

**Sample 105:**
```bash
nohup bash -c "
# Alignment
bwa mem -t 16 -M reference.fasta \
  trimmed/sample105.R1.trim.fastq.gz \
  trimmed/sample105.R2.trim.fastq.gz | \
samtools view -@ 8 -bS | \
samtools sort -@ 8 -o aligned/sample105.sorted.bam

# Add read groups
gatk AddOrReplaceReadGroups \
  -I aligned/sample105.sorted.bam \
  -O aligned/sample105.RG.bam \
  -RGID sample105 \
  -RGLB lib1 \
  -RGPL ILLUMINA \
  -RGPU unit1 \
  -RGSM sample105

# Index BAM
samtools index aligned/sample105.RG.bam

# Clean up intermediate files
rm aligned/sample105.sorted.bam
" > logs/sample105.align.log 2>&1 &
```

**Sample 106:**
```bash
nohup bash -c "
# Alignment
bwa mem -t 16 -M reference.fasta \
  trimmed/sample106.R1.trim.fastq.gz \
  trimmed/sample106.R2.trim.fastq.gz | \
samtools view -@ 8 -bS | \
samtools sort -@ 8 -o aligned/sample106.sorted.bam

# Add read groups
gatk AddOrReplaceReadGroups \
  -I aligned/sample106.sorted.bam \
  -O aligned/sample106.RG.bam \
  -RGID sample106 \
  -RGLB lib1 \
  -RGPL ILLUMINA \
  -RGPU unit1 \
  -RGSM sample106

# Index BAM
samtools index aligned/sample106.RG.bam

# Clean up intermediate files
rm aligned/sample106.sorted.bam
" > logs/sample106.align.log 2>&1 &
```

**Sample 107:**
```bash
nohup bash -c "
# Alignment
bwa mem -t 16 -M reference.fasta \
  trimmed/sample107.R1.trim.fastq.gz \
  trimmed/sample107.R2.trim.fastq.gz | \
samtools view -@ 8 -bS | \
samtools sort -@ 8 -o aligned/sample107.sorted.bam

# Add read groups
gatk AddOrReplaceReadGroups \
  -I aligned/sample107.sorted.bam \
  -O aligned/sample107.RG.bam \
  -RGID sample107 \
  -RGLB lib1 \
  -RGPL ILLUMINA \
  -RGPU unit1 \
  -RGSM sample107

# Index BAM
samtools index aligned/sample107.RG.bam

# Clean up intermediate files
rm aligned/sample107.sorted.bam
" > logs/sample107.align.log 2>&1 &
```

**Sample 108:**
```bash
nohup bash -c "
# Alignment
bwa mem -t 16 -M reference.fasta \
  trimmed/sample108.R1.trim.fastq.gz \
  trimmed/sample108.R2.trim.fastq.gz | \
samtools view -@ 8 -bS | \
samtools sort -@ 8 -o aligned/sample108.sorted.bam

# Add read groups
gatk AddOrReplaceReadGroups \
  -I aligned/sample108.sorted.bam \
  -O aligned/sample108.RG.bam \
  -RGID sample108 \
  -RGLB lib1 \
  -RGPL ILLUMINA \
  -RGPU unit1 \
  -RGSM sample108

# Index BAM
samtools index aligned/sample108.RG.bam

# Clean up intermediate files
rm aligned/sample108.sorted.bam
" > logs/sample108.align.log 2>&1 &
```

```bash
# Wait for all alignment jobs to complete
wait
echo "Alignment completed for all samples"
```

---

## Step 4: Individual GVCF Generation

**Sample 105:**
```bash
nohup gatk HaplotypeCaller \
  -R reference.fasta \
  -I aligned/sample105.RG.bam \
  -O variants/sample105.g.vcf.gz \
  -ERC GVCF \
  --sample-ploidy 4 \
  --native-pair-hmm-threads 4 \
  --min-base-quality-score 20 \
  > logs/sample105.gvcf.log 2>&1 &
```

**Sample 106:**
```bash
nohup gatk HaplotypeCaller \
  -R reference.fasta \
  -I aligned/sample106.RG.bam \
  -O variants/sample106.g.vcf.gz \
  -ERC GVCF \
  --sample-ploidy 4 \
  --native-pair-hmm-threads 4 \
  --min-base-quality-score 20 \
  > logs/sample106.gvcf.log 2>&1 &
```

**Sample 107:**
```bash
nohup gatk HaplotypeCaller \
  -R reference.fasta \
  -I aligned/sample107.RG.bam \
  -O variants/sample107.g.vcf.gz \
  -ERC GVCF \
  --sample-ploidy 4 \
  --native-pair-hmm-threads 4 \
  --min-base-quality-score 20 \
  > logs/sample107.gvcf.log 2>&1 &
```

**Sample 108:**
```bash
nohup gatk HaplotypeCaller \
  -R reference.fasta \
  -I aligned/sample108.RG.bam \
  -O variants/sample108.g.vcf.gz \
  -ERC GVCF \
  --sample-ploidy 4 \
  --native-pair-hmm-threads 4 \
  --min-base-quality-score 20 \
  > logs/sample108.gvcf.log 2>&1 &
```

```bash
# Wait for all GVCF generation to complete
wait
echo "GVCF generation completed for all samples"

# Index all GVCF files
for sample in sample105 sample106 sample107 sample108; do
    gatk IndexFeatureFile -I variants/${sample}.g.vcf.gz
done
```

---

## Step 5: Simple Joint Calling

```bash
# Combine all GVCFs
nohup gatk CombineGVCFs \
  -R reference.fasta \
  --variant variants/sample105.g.vcf.gz \
  --variant variants/sample106.g.vcf.gz \
  --variant variants/sample107.g.vcf.gz \
  --variant variants/sample108.g.vcf.gz \
  -O variants/cohort.combined.g.vcf.gz \
  > logs/combine_gvcfs.log 2>&1 &

wait

# Joint genotyping
nohup gatk GenotypeGVCFs \
  -R reference.fasta \
  -V variants/cohort.combined.g.vcf.gz \
  -O variants/cohort.raw.vcf.gz \
  > logs/joint_genotyping.log 2>&1 &

wait
echo "Joint genotyping completed"
```

---

## Step 6: Hard Filtering

```bash
# Apply hard filters for tetraploid data
nohup gatk VariantFiltration \
  -R reference.fasta \
  -V variants/cohort.raw.vcf.gz \
  --filter-expression "QD < 2.0" \
  --filter-name "QD2" \
  --filter-expression "FS > 60.0" \
  --filter-name "FS60" \
  --filter-expression "MQ < 40.0" \
  --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5" \
  --filter-name "MQRankSum-12.5" \
  --filter-expression "ReadPosRankSum < -8.0" \
  --filter-name "ReadPosRankSum-8" \
  --filter-expression "SOR > 3.0" \
  --filter-name "SOR3" \
  -O variants/cohort.filtered.vcf.gz \
  > logs/filtering.log 2>&1 &

wait

# Select only PASS variants
nohup gatk SelectVariants \
  -R reference.fasta \
  -V variants/cohort.filtered.vcf.gz \
  --exclude-filtered \
  -O variants/cohort.final.vcf.gz \
  > logs/select_variants.log 2>&1 &

wait
echo "Filtering completed"
```

---

## Step 7: Final Statistics

```bash
# Generate statistics
nohup bash -c "
# Variant counts
bcftools stats variants/cohort.final.vcf.gz > variants/cohort.stats.txt

# Extract biallelic SNPs
gatk SelectVariants \
  -R reference.fasta \
  -V variants/cohort.final.vcf.gz \
  --select-type-to-include SNP \
  --restrict-alleles-to BIALLELIC \
  --min-allele-frequency 0.05 \
  -O variants/cohort.biallelic.snps.vcf.gz

echo 'Pipeline completed successfully!'
echo 'Final VCF: variants/cohort.final.vcf.gz'
echo 'Biallelic SNPs: variants/cohort.biallelic.snps.vcf.gz'
echo 'Statistics: variants/cohort.stats.txt'
" > logs/final_stats.log 2>&1 &
```

---

## Complete Automated Script

**Save as `run_pipeline.sh`:**

```bash
#!/bin/bash

# Set memory for GATK
export GATK_LOCAL_JAR="--java-options '-Xmx32g'"

echo "Starting Brassica variant calling pipeline..."
echo "Start time: $(date)"

# Create directories
mkdir -p brassica_analysis/{raw_data,trimmed,aligned,variants,logs}
cd brassica_analysis

# Check if sample files exist
samples=("sample105" "sample106" "sample107" "sample108")
for sample in "${samples[@]}"; do
    if [[ ! -f "raw_data/${sample}.R1.fastq.gz" ]]; then
        echo "Error: raw_data/${sample}.R1.fastq.gz not found"
        exit 1
    fi
done

echo "All input files found. Starting pipeline..."

# Step 1: Quality control (parallel)
echo "Step 1: Quality control - $(date)"
for sample in "${samples[@]}"; do
    nohup fastp \
      -i raw_data/${sample}.R1.fastq.gz \
      -I raw_data/${sample}.R2.fastq.gz \
      -o trimmed/${sample}.R1.trim.fastq.gz \
      -O trimmed/${sample}.R2.trim.fastq.gz \
      --detect_adapter_for_pe \
      --qualified_quality_phred 20 \
      --length_required 50 \
      --thread 8 \
      --html trimmed/${sample}.fastp.html \
      --json trimmed/${sample}.fastp.json \
      > logs/${sample}.fastp.log 2>&1 &
done
wait

# Step 2: Reference preparation
echo "Step 2: Reference preparation - $(date)"
if [[ ! -f "reference.fasta" ]]; then
    nohup bash -c "
    datasets download genome accession GCA_015484525.1 --include genome
    unzip ncbi_dataset.zip
    mv ncbi_dataset/data/GCA_015484525.1/GCA_015484525.1_UDSC_Var_1.1_genomic.fna reference.fasta
    bwa index reference.fasta
    samtools faidx reference.fasta
    gatk CreateSequenceDictionary -R reference.fasta -O reference.dict
    " > logs/reference_prep.log 2>&1 &
    wait
fi

# Step 3: Alignment (parallel)
echo "Step 3: Alignment - $(date)"
for sample in "${samples[@]}"; do
    nohup bash -c "
    bwa mem -t 4 -M reference.fasta \
      trimmed/${sample}.R1.trim.fastq.gz \
      trimmed/${sample}.R2.trim.fastq.gz | \
    samtools view -@ 2 -bS | \
    samtools sort -@ 2 -o aligned/${sample}.sorted.bam

    gatk AddOrReplaceReadGroups \
      -I aligned/${sample}.sorted.bam \
      -O aligned/${sample}.RG.bam \
      -RGID ${sample} \
      -RGLB lib1 \
      -RGPL ILLUMINA \
      -RGPU unit1 \
      -RGSM ${sample}

    samtools index aligned/${sample}.RG.bam
    rm aligned/${sample}.sorted.bam
    " > logs/${sample}.align.log 2>&1 &
done
wait

# Step 4: GVCF generation (parallel)
echo "Step 4: GVCF generation - $(date)"
for sample in "${samples[@]}"; do
    nohup gatk HaplotypeCaller \
      -R reference.fasta \
      -I aligned/${sample}.RG.bam \
      -O variants/${sample}.g.vcf.gz \
      -ERC GVCF \
      --sample-ploidy 4 \
      --native-pair-hmm-threads 2 \
      --min-base-quality-score 20 \
      > logs/${sample}.gvcf.log 2>&1 &
done
wait

# Index GVCFs
for sample in "${samples[@]}"; do
    gatk IndexFeatureFile -I variants/${sample}.g.vcf.gz
done

# Step 5: Joint calling
echo "Step 5: Joint calling - $(date)"
nohup gatk CombineGVCFs \
  -R reference.fasta \
  --variant variants/sample105.g.vcf.gz \
  --variant variants/sample106.g.vcf.gz \
  --variant variants/sample107.g.vcf.gz \
  --variant variants/sample108.g.vcf.gz \
  -O variants/cohort.combined.g.vcf.gz \
  > logs/combine_gvcfs.log 2>&1 &
wait

nohup gatk GenotypeGVCFs \
  -R reference.fasta \
  -V variants/cohort.combined.g.vcf.gz \
  -O variants/cohort.raw.vcf.gz \
  > logs/joint_genotyping.log 2>&1 &
wait

# Step 6: Filtering
echo "Step 6: Filtering - $(date)"
nohup gatk VariantFiltration \
  -R reference.fasta \
  -V variants/cohort.raw.vcf.gz \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 60.0" --filter-name "FS60" \
  --filter-expression "MQ < 40.0" --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  --filter-expression "SOR > 3.0" --filter-name "SOR3" \
  -O variants/cohort.filtered.vcf.gz \
  > logs/filtering.log 2>&1 &
wait

nohup gatk SelectVariants \
  -R reference.fasta \
  -V variants/cohort.filtered.vcf.gz \
  --exclude-filtered \
  -O variants/cohort.final.vcf.gz \
  > logs/select_variants.log 2>&1 &
wait

# Step 7: Statistics
echo "Step 7: Final statistics - $(date)"
bcftools stats variants/cohort.final.vcf.gz > variants/cohort.stats.txt

gatk SelectVariants \
  -R reference.fasta \
  -V variants/cohort.final.vcf.gz \
  --select-type-to-include SNP \
  --restrict-alleles-to BIALLELIC \
  --min-allele-frequency 0.05 \
  -O variants/cohort.biallelic.snps.vcf.gz

echo "Pipeline completed successfully!"
echo "End time: $(date)"
echo "Final outputs:"
echo "  - variants/cohort.final.vcf.gz (all filtered variants)"
echo "  - variants/cohort.biallelic.snps.vcf.gz (biallelic SNPs)"
echo "  - variants/cohort.stats.txt (statistics)"
```

**To run the complete pipeline:**
```bash
chmod +x run_pipeline.sh
nohup ./run_pipeline.sh > pipeline.log 2>&1 &

# Monitor progress
tail -f pipeline.log
```

---

## Key Features:

1. **All commands use nohup** for server execution
2. **Simple joint calling** with CombineGVCFs
3. **Four sample pipeline** ready to use
4. **Parallel processing** where possible
5. **Complete automation** with the master script
6. **Tetraploid optimized** (ploidy=4)
7. **Log files** for each step
8. **Memory management** for server environment

Just modify the sample names in the script to match your actual sample files!
