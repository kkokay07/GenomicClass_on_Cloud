Comprehensive GATK Variant Calling Pipeline for Tetraploid Brassica

## Overview

This document provides a complete methodology for variant calling in tetraploid Brassica species using the GATK (Genome Analysis Toolkit) pipeline. The workflow is optimized for ddRAD-seq data with ploidy level 4, ensuring accurate identification of genetic variants across multiple samples.

---

## Prerequisites

### Required Software
- **GATK** (≥4.2.0)
- **BWA** or **Bowtie2** for alignment
- **SAMtools** (≥1.10)
- **Picard Tools**
- **fastp** for quality control
- **Stacks** (for RAD-tag processing)

### Hardware Requirements
- Minimum 16 CPU cores recommended
- 32GB+ RAM for large datasets
- Sufficient storage for intermediate files

---

## Pipeline Workflow

### Step 1: Quality Control and Read Trimming

**Purpose**: Remove low-quality bases, adapters, and short reads to ensure high-quality input for alignment.

```bash
# Quality trimming with fastp
fastp \
  -i sample105.R1.fastq.gz \
  -I sample105.R2.fastq.gz \
  -o sample105.R1.trim.fastq.gz \
  -O sample105.R2.trim.fastq.gz \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --thread 8 \
  --html sample105.fastp.html \
  --json sample105.fastp.json
```

**Parameter Explanation**:
- `-i/-I`: Input forward/reverse reads
- `-o/-O`: Output trimmed forward/reverse reads
- `--detect_adapter_for_pe`: Automatically detects and removes adapters
- `--qualified_quality_phred 20`: Filters bases with quality score <20
- `--length_required 50`: Discards reads shorter than 50bp after trimming
- `--thread 8`: Uses 8 CPU threads for parallel processing

---

### Step 2: Reference Genome Preparation

**Purpose**: Prepare the reference genome for alignment and variant calling.

```bash
# Download Brassica reference genome
datasets download genome accession GCA_015484525.1 --include genome

# Extract and rename
unzip ncbi_dataset.zip
mv ncbi_dataset/data/GCA_015484525.1/GCA_015484525.1_UDSC_Var_1.1_genomic.fna reference.fasta

# Index reference genome for BWA
bwa index reference.fasta

# Alternative: Index for Bowtie2
bowtie2-build reference.fasta reference_index

# Create FASTA index for GATK
samtools faidx reference.fasta

# Create sequence dictionary for GATK
gatk CreateSequenceDictionary \
  -R reference.fasta \
  -O reference.dict
```

**Why these steps matter**:
- BWA/Bowtie2 indices enable rapid read alignment
- FASTA index allows random access to genome sequences
- Sequence dictionary provides chromosome information for GATK

---

### Step 3: RAD-tag Processing (Optional)

**Purpose**: Remove restriction enzyme sites and demultiplex samples if using RAD-seq data.

```bash
# Process RAD tags using Stacks
process_radtags \
  -1 sample105.R1.trim.fastq.gz \
  -2 sample105.R2.trim.fastq.gz \
  -o ./radtags_output/ \
  -b barcodes_file.txt \
  -e sphI \
  --renz_2 mluCI \
  -c -q -r \
  -t 140 \
  --paired \
  --inline_null
```

**Parameter Details**:
- `-e sphI --renz_2 mluCI`: Restriction enzymes used in ddRAD protocol
- `-c`: Remove reads with uncalled bases (Ns)
- `-q`: Discard low-quality reads
- `-r`: Rescue barcodes with sequencing errors
- `-t 140`: Truncate reads to 140bp (adjust based on quality)
- `--inline_null`: For inline barcodes without separate index

---

### Step 4: Read Alignment

**Purpose**: Map quality-filtered reads to the reference genome.

#### Option A: Using BWA-MEM (Recommended for longer reads)
```bash
# Align reads using BWA-MEM
bwa mem -t 16 -M \
  reference.fasta \
  sample105.R1.trim.fastq.gz \
  sample105.R2.trim.fastq.gz \
  > sample105.sam
```

#### Option B: Using Bowtie2 (Good for shorter reads)
```bash
# Align reads using Bowtie2
bowtie2 \
  -p 16 \
  --very-sensitive-local \
  -x reference_index \
  -1 sample105.R1.trim.fastq.gz \
  -2 sample105.R2.trim.fastq.gz \
  -S sample105.sam
```

**Parameter Explanation**:
- `-t 16/-p 16`: Use 16 CPU threads
- `-M`: Mark shorter split hits as secondary (BWA)
- `--very-sensitive-local`: Thorough local alignment (Bowtie2)

---

### Step 5: SAM/BAM Processing

**Purpose**: Convert, sort, and prepare alignment files for variant calling.

```bash
# Convert SAM to BAM and sort in one step
samtools view -@ 8 -bS sample105.sam | \
samtools sort -@ 8 -o sample105.sorted.bam

# Index the sorted BAM file
samtools index sample105.sorted.bam

# Add read groups (required for GATK)
gatk AddOrReplaceReadGroups \
  -I sample105.sorted.bam \
  -O sample105.RG.bam \
  -RGID sample105 \
  -RGLB lib1 \
  -RGPL ILLUMINA \
  -RGPU unit1 \
  -RGSM sample105

# Index the read group BAM
samtools index sample105.RG.bam

# Optional: Remove duplicates
gatk MarkDuplicates \
  -I sample105.RG.bam \
  -O sample105.dedup.bam \
  -M sample105.metrics.txt

samtools index sample105.dedup.bam
```

**Read Group Information**:
- `RGID`: Unique identifier for the read group
- `RGLB`: Library identifier
- `RGPL`: Sequencing platform
- `RGPU`: Platform unit (e.g., flowcell-barcode.lane)
- `RGSM`: Sample name

---

### Step 6: Quality Assessment

**Purpose**: Evaluate alignment quality before variant calling.

```bash
# Generate alignment statistics
samtools flagstat sample105.RG.bam > sample105.flagstat.txt

# Collect alignment summary metrics
gatk CollectAlignmentSummaryMetrics \
  -R reference.fasta \
  -I sample105.RG.bam \
  -O sample105.alignment_metrics.txt

# Validate SAM file
gatk ValidateSamFile \
  -I sample105.RG.bam \
  -MODE SUMMARY
```

---

### Step 7: Variant Calling (GVCF Generation)

**Purpose**: Call variants for each sample individually, accounting for tetraploid nature.

```bash
# Generate GVCF for tetraploid sample
gatk HaplotypeCaller \
  -R reference.fasta \
  -I sample105.RG.bam \
  -O sample105.g.vcf.gz \
  -ERC GVCF \
  --sample-ploidy 4 \
  --native-pair-hmm-threads 4 \
  --min-base-quality-score 20 \
  --dont-use-soft-clipped-bases true

# Index the GVCF file
gatk IndexFeatureFile -I sample105.g.vcf.gz
```

**Critical Parameters for Tetraploids**:
- `--sample-ploidy 4`: Essential for tetraploid Brassica
- `--native-pair-hmm-threads 4`: Parallel processing for HMM calculations
- `--min-base-quality-score 20`: Filter low-quality bases
- `--dont-use-soft-clipped-bases`: Improves accuracy for repetitive regions

---

### Step 8: Joint Genotyping

**Purpose**: Combine all samples for accurate joint variant calling.

```bash
# Create sample map file
echo -e "sample105\tsample105.g.vcf.gz" > cohort.sample_map
echo -e "sample106\tsample106.g.vcf.gz" >> cohort.sample_map
# ... add all samples

# Combine GVCFs using GenomicsDB
gatk GenomicsDBImport \
  --genomicsdb-workspace-path cohort_genomicsdb \
  --batch-size 50 \
  -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 \
  --sample-name-map cohort.sample_map \
  --tmp-dir tmp/

# Joint genotyping
gatk GenotypeGVCFs \
  -R reference.fasta \
  -V gendb://cohort_genomicsdb \
  -O cohort.raw.vcf.gz \
  --tmp-dir tmp/
```

**Alternative method for smaller datasets**:
```bash
# Direct GVCF combination (for <100 samples)
gatk CombineGVCFs \
  -R reference.fasta \
  --variant sample105.g.vcf.gz \
  --variant sample106.g.vcf.gz \
  -O cohort.combined.g.vcf.gz

gatk GenotypeGVCFs \
  -R reference.fasta \
  -V cohort.combined.g.vcf.gz \
  -O cohort.raw.vcf.gz
```

---

### Step 9: Variant Quality Score Recalibration (VQSR)

**Purpose**: Apply machine learning-based quality scoring for high-confidence variants.

**Note**: VQSR requires large datasets (>30 samples, >1M variants). For smaller datasets, use hard filtering (Step 10).

```bash
# Build recalibration model for SNPs
gatk VariantRecalibrator \
  -R reference.fasta \
  -V cohort.raw.vcf.gz \
  --resource:training,known=false,truth=true,prior=15.0 high_confidence_snps.vcf \
  -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
  -mode SNP \
  -O cohort.snps.recal \
  --tranches-file cohort.snps.tranches

# Apply recalibration to SNPs
gatk ApplyVQSR \
  -R reference.fasta \
  -V cohort.raw.vcf.gz \
  --recal-file cohort.snps.recal \
  --tranches-file cohort.snps.tranches \
  -mode SNP \
  --truth-sensitivity-filter-level 99.0 \
  -O cohort.snps.vcf.gz
```

---

### Step 10: Hard Filtering (Alternative to VQSR)

**Purpose**: Apply stringent filters for high-quality variants when VQSR is not applicable.

```bash
# Hard filtering for tetraploid data
gatk VariantFiltration \
  -R reference.fasta \
  -V cohort.raw.vcf.gz \
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
  -O cohort.filtered.vcf.gz

# Select only PASS variants
gatk SelectVariants \
  -R reference.fasta \
  -V cohort.filtered.vcf.gz \
  --exclude-filtered \
  -O cohort.pass.vcf.gz
```

**Filter Explanations**:
- **QD** (Quality by Depth): Variant confidence normalized by depth
- **FS** (Fisher Strand): Strand bias using Fisher's exact test
- **MQ** (Mapping Quality): Root mean square mapping quality
- **SOR** (Strand Odds Ratio): Strand bias using symmetric odds ratio

---

### Step 11: Final Quality Control and Statistics

```bash
# Generate variant statistics
gatk VariantEval \
  -R reference.fasta \
  -eval cohort.pass.vcf.gz \
  -O cohort.eval.txt

# Count variants
bcftools stats cohort.pass.vcf.gz > cohort.stats.txt

# Extract high-quality biallelic SNPs
gatk SelectVariants \
  -R reference.fasta \
  -V cohort.pass.vcf.gz \
  --select-type-to-include SNP \
  --restrict-alleles-to BIALLELIC \
  --min-allele-frequency 0.05 \
  -O cohort.biallelic.snps.vcf.gz
```

---

## Batch Processing Script

For processing multiple samples efficiently:

```bash
#!/bin/bash

# Sample list file (one sample name per line)
SAMPLES="sample_list.txt"
REF="reference.fasta"
THREADS=16

while read SAMPLE; do
    echo "Processing $SAMPLE..."
    
    # Quality control
    fastp -i ${SAMPLE}.R1.fastq.gz -I ${SAMPLE}.R2.fastq.gz \
          -o ${SAMPLE}.R1.trim.fastq.gz -O ${SAMPLE}.R2.trim.fastq.gz \
          --detect_adapter_for_pe --qualified_quality_phred 20 \
          --length_required 50 --thread 8
    
    # Alignment
    bwa mem -t $THREADS -M $REF \
        ${SAMPLE}.R1.trim.fastq.gz ${SAMPLE}.R2.trim.fastq.gz | \
    samtools view -@ 8 -bS | \
    samtools sort -@ 8 -o ${SAMPLE}.sorted.bam
    
    # Add read groups
    gatk AddOrReplaceReadGroups \
        -I ${SAMPLE}.sorted.bam -O ${SAMPLE}.RG.bam \
        -RGID $SAMPLE -RGLB lib1 -RGPL ILLUMINA \
        -RGPU unit1 -RGSM $SAMPLE
    
    samtools index ${SAMPLE}.RG.bam
    
    # Variant calling
    gatk HaplotypeCaller \
        -R $REF -I ${SAMPLE}.RG.bam \
        -O ${SAMPLE}.g.vcf.gz -ERC GVCF \
        --sample-ploidy 4 --native-pair-hmm-threads 4
    
    gatk IndexFeatureFile -I ${SAMPLE}.g.vcf.gz
    
    echo "$SAMPLE processing complete!"
    
done < $SAMPLES
```

---

## Troubleshooting Tips

### Common Issues and Solutions

1. **Memory Errors**: Increase JVM heap size
   ```bash
   export GATK_LOCAL_JAR="--java-options '-Xmx32g'"
   ```

2. **Slow Performance**: Use parallel processing and temporary directories
   ```bash
   --tmp-dir /fast/tmp/directory
   ```

3. **Ploidy-specific Issues**: Ensure `--sample-ploidy 4` is consistently used

4. **Low Variant Counts**: Check alignment quality and filtering stringency

### Quality Metrics to Monitor

- **Alignment Rate**: Should be >85% for good quality data
- **Duplicate Rate**: Should be <20% for diverse libraries
- **Ti/Tv Ratio**: Should be ~2.0-2.1 for high-quality SNPs
- **Het/Hom Ratio**: Expected to be higher in tetraploids

---

## Expected Outputs

- **Raw VCF**: All called variants before filtering
- **Filtered VCF**: High-confidence variants passing quality filters
- **Statistics Files**: Alignment and variant calling metrics
- **Quality Reports**: HTML/JSON reports from fastp and GATK

This comprehensive pipeline ensures robust variant calling for tetraploid Brassica species while maintaining computational efficiency and accuracy.
