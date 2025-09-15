# Detailed Description of Variant Calling Pipeline Steps

## Overview

This document provides comprehensive descriptions of each step in the GATK variant calling pipeline for tetraploid Brassica species. Understanding these steps is crucial for troubleshooting, optimization, and interpreting results.

---

## Step 1: Quality Control and Read Trimming

### Purpose
Raw sequencing data often contains low-quality bases, adapter sequences, and reads that are too short for reliable alignment. Quality control ensures that only high-quality sequence data proceeds to downstream analysis.

### What Happens
- **Adapter Detection**: Automatically identifies and removes adapter sequences from read ends
- **Quality Filtering**: Removes bases with quality scores below threshold (Phred score < 20)
- **Length Filtering**: Discards reads shorter than 50 bp after trimming
- **Quality Reports**: Generates HTML and JSON reports showing before/after statistics

### Key Parameters Explained
- `--qualified_quality_phred 20`: Sets minimum quality score (99% base calling accuracy)
- `--length_required 50`: Ensures sufficient sequence for reliable alignment
- `--detect_adapter_for_pe`: Automatically finds adapter sequences in paired-end data
- `--thread 8`: Uses 8 CPU cores for faster processing

### Expected Outcomes
- **Good Quality Data**: >90% of reads should pass filtering
- **Adapter Removal**: Should see significant reduction in adapter content
- **Length Distribution**: Most reads should retain >80% of original length
- **Base Quality**: Mean quality scores should improve significantly

### Troubleshooting
- **High Read Loss**: May indicate poor sequencing quality or incorrect adapter sequences
- **Low Quality Improvement**: Could suggest fundamental sequencing issues
- **Short Reads**: May require adjusting length threshold for ddRAD data

---

## Step 2: Reference Genome Preparation

### Purpose
The reference genome serves as a template for read alignment and variant calling. Proper indexing is essential for efficient processing and GATK compatibility.

### What Happens
- **Download**: Retrieves the latest Brassica reference genome from NCBI
- **BWA Indexing**: Creates suffix array indices for rapid read alignment
- **FASTA Indexing**: Enables random access to genome sequences
- **Dictionary Creation**: Provides chromosome information required by GATK

### Index Files Created
- `reference.fasta.amb, .ann, .bwt, .pac, .sa`: BWA alignment indices
- `reference.fasta.fai`: FASTA index for sequence retrieval
- `reference.dict`: Sequence dictionary with chromosome metadata

### Why Multiple Indices
- **BWA Index**: Enables fast approximate string matching for alignment
- **FASTA Index**: Allows tools to quickly access specific genomic regions
- **Dictionary**: Provides chromosome names, lengths, and checksums for validation

### Quality Checks
- Verify all index files are created without errors
- Check that reference.dict contains expected chromosome information
- Ensure file sizes are reasonable (indices should be fraction of genome size)

---

## Step 3: Read Alignment

### Purpose
Maps quality-filtered reads to specific positions in the reference genome, creating the foundation for variant detection.

### What Happens During BWA-MEM Alignment
1. **Seed Finding**: Identifies exact matches between reads and reference
2. **Seed Extension**: Extends matches allowing for mismatches and gaps
3. **Alignment Scoring**: Calculates alignment quality based on matches/mismatches
4. **Secondary Alignment**: Handles reads with multiple mapping locations

### SAM to BAM Conversion Process
- **SAM Format**: Human-readable text format with alignment information
- **BAM Format**: Binary compressed version for efficient storage and processing
- **Sorting**: Arranges alignments by genomic coordinates for indexed access
- **Indexing**: Creates lookup table for rapid region-specific queries

### Read Group Addition
Read groups are essential metadata that GATK requires for proper variant calling:
- **RGID**: Unique identifier for this sequencing run
- **RGLB**: Library preparation identifier
- **RGPL**: Sequencing platform (ILLUMINA)
- **RGPU**: Platform unit (flowcell-barcode.lane)
- **RGSM**: Sample identifier

### Quality Metrics to Monitor
- **Alignment Rate**: Should be >85% for good quality data
- **Properly Paired**: >95% of paired reads should align in proper orientation
- **Mapping Quality**: Mean mapping quality should be >20
- **Insert Size**: Should match expected library fragment size

### Common Issues
- **Low Alignment Rate**: May indicate wrong reference genome or poor quality data
- **High Secondary Alignments**: Could suggest repetitive sequences or contamination
- **Abnormal Insert Sizes**: Might indicate library preparation problems

---

## Step 4: Individual GVCF Generation

### Purpose
Calls variants for each sample individually while preserving information about non-variant sites, enabling accurate joint genotyping across multiple samples.

### GVCF vs Regular VCF
- **Regular VCF**: Only contains variant sites
- **GVCF**: Contains both variant and non-variant sites with confidence scores
- **Benefit**: Allows accurate joint calling by considering evidence from all samples

### HaplotypeCaller Process
1. **Local Assembly**: Reconstructs haplotypes in regions with variation
2. **Realignment**: Re-aligns reads to candidate haplotypes
3. **Genotype Likelihood**: Calculates probability of each possible genotype
4. **Variant Calling**: Identifies sites that differ from reference

### Tetraploid-Specific Considerations
- `--sample-ploidy 4`: Critical for accurate genotype calling in tetraploids
- **Allele Frequencies**: Can have 0/4, 1/4, 2/4, 3/4, 4/4 allele ratios
- **Complexity**: More possible genotype combinations than diploids
- **Coverage Requirements**: May need higher coverage for confident calling

### Quality Parameters
- `--min-base-quality-score 20`: Filters low-quality bases from analysis
- `--native-pair-hmm-threads 4`: Parallel processing for computational efficiency
- **Output Compression**: GVCF files are compressed to save storage space

### Expected File Sizes
- GVCF files are typically 2-5x larger than regular VCF files
- Size depends on genome coverage and variant density
- Files should be successfully indexed for downstream processing

---

## Step 5: Simple Joint Calling

### Purpose
Combines individual sample data to make more accurate variant calls by leveraging information across all samples in the cohort.

### Two-Step Process

#### Step 5a: CombineGVCFs
- **Function**: Merges individual GVCF files into a single multi-sample GVCF
- **Advantage**: Simpler than GenomicsDB for small datasets (<100 samples)
- **Output**: Single combined GVCF containing all samples
- **Memory Usage**: Moderate - scales linearly with sample number

#### Step 5b: GenotypeGVCFs
- **Function**: Performs joint genotyping across all samples simultaneously
- **Process**: Re-evaluates genotype calls using population-level information
- **Benefit**: More accurate calls, especially for rare variants
- **Population Genetics**: Can detect variants missed in individual calling

### Why Joint Calling is Superior
1. **Increased Sensitivity**: Detects variants with low coverage in individual samples
2. **Improved Accuracy**: Uses population allele frequencies for better genotyping
3. **Consistent Calls**: Ensures all samples are called at the same sites
4. **Better Filtering**: Population-level metrics improve variant quality assessment

### Computational Considerations
- **Memory Requirements**: Increases with sample number and genome size
- **Processing Time**: Longer than individual calling but more accurate
- **Temporary Files**: May create large temporary files during processing
- **I/O Intensive**: Reads all GVCF files simultaneously

### Quality Indicators
- **Variant Count**: Should be reasonable for species and population
- **Sample Coverage**: All samples should contribute variants
- **Missing Data**: Excessive missing calls may indicate processing issues

---

## Step 6: Hard Filtering

### Purpose
Applies stringent quality filters to retain only high-confidence variants while removing potential false positives.

### Filter Criteria Explained

#### QD (Quality by Depth) < 2.0
- **Meaning**: Variant quality score normalized by allelic depth
- **Purpose**: Removes variants with low quality relative to coverage
- **Why Important**: High coverage with low quality suggests systematic errors

#### FS (Fisher Strand) > 60.0
- **Meaning**: Phred-scaled probability of strand bias using Fisher's exact test
- **Purpose**: Identifies variants called predominantly on one DNA strand
- **Why Important**: True variants should be supported by both strands

#### MQ (Mapping Quality) < 40.0
- **Meaning**: Root mean square of mapping quality scores
- **Purpose**: Removes variants in regions with poor alignment quality
- **Why Important**: Low mapping quality suggests repetitive or problematic regions

#### MQRankSum < -12.5
- **Meaning**: Compares mapping quality between reference and alternate alleles
- **Purpose**: Identifies systematic differences in mapping quality
- **Why Important**: True variants shouldn't show mapping quality bias

#### ReadPosRankSum < -8.0
- **Meaning**: Tests for bias in variant position within reads
- **Purpose**: Removes variants consistently at read ends
- **Why Important**: End-of-read variants are often sequencing artifacts

#### SOR (Strand Odds Ratio) > 3.0
- **Meaning**: Symmetric odds ratio test for strand bias
- **Purpose**: More sensitive strand bias detection than Fisher Strand
- **Why Important**: Complements FS for comprehensive strand bias filtering

### Filter Strategy
- **Conservative Approach**: Better to remove borderline variants than include false positives
- **Tetraploid Adjustment**: Thresholds may need modification for polyploid-specific patterns
- **Population-Specific**: Different populations may require different thresholds

### Post-Filtering Steps
- **PASS Selection**: Only variants passing all filters are retained
- **Quality Assessment**: Statistics on filtered vs. unfiltered variants
- **Biallelic Selection**: Focus on simple two-allele variants for most analyses

---

## Step 7: Final Quality Control and Statistics

### Purpose
Validates the pipeline results and generates comprehensive statistics for publication and downstream analysis.

### Statistics Generated

#### Variant Counts
- **Total Variants**: SNPs, indels, and complex variants
- **SNP/Indel Ratio**: Should be approximately 3:1 to 5:1
- **Ti/Tv Ratio**: Transition to transversion ratio (should be ~2.0-2.1)
- **Heterozygote/Homozygote Ratio**: Expected to be higher in tetraploids

#### Quality Metrics
- **Mean Quality Score**: Average variant quality across dataset
- **Depth Distribution**: Coverage statistics per variant and sample
- **Missing Data**: Proportion of no-calls per sample and variant
- **Allele Frequency Spectrum**: Distribution of minor allele frequencies

#### Sample-Level Statistics
- **Variants per Sample**: Individual sample contribution to dataset
- **Genotype Quality**: Distribution of genotype confidence scores
- **Coverage Statistics**: Mean, median, and distribution of read depth
- **Callable Genome**: Proportion of genome with adequate coverage

### Biallelic SNP Selection
- **Purpose**: Creates clean dataset for population genetics analyses
- **Criteria**: Only two alleles per site, minimum allele frequency
- **Applications**: GWAS, population structure, phylogenetics
- **Quality**: Highest confidence subset of variants

### Expected Results for Brassica
- **Variant Density**: Approximately 1-10 SNPs per kb depending on diversity
- **Callable Sites**: 70-90% of genome depending on coverage
- **Ti/Tv Ratio**: 2.0-2.1 for high-quality SNP datasets
- **Missing Data**: <10% for well-covered samples

### Validation Checks
- **Mendelian Consistency**: If families available, check inheritance patterns
- **Hardy-Weinberg Equilibrium**: Population-level genotype frequencies
- **Linkage Disequilibrium**: Expected patterns for the organism
- **Sample Relationships**: Verify expected kinship relationships

---

## Pipeline Success Indicators

### Overall Quality Metrics
1. **High Alignment Rate**: >85% of reads successfully aligned
2. **Reasonable Variant Count**: 100K-10M variants depending on genome size and diversity
3. **Good Ti/Tv Ratio**: 2.0-2.1 for genome-wide SNPs
4. **Low Missing Data**: <10% no-calls across samples
5. **Consistent Sample Performance**: All samples contribute similar variant numbers

### Red Flags to Watch For
- **Very Low Alignment Rates**: May indicate wrong reference or contamination
- **Extreme Variant Counts**: Too few suggests poor calling; too many suggests false positives
- **Abnormal Ti/Tv Ratios**: <1.5 or >3.0 suggests quality issues
- **High Missing Data**: >20% suggests coverage or quality problems
- **Sample Outliers**: One sample with very different metrics

### Downstream Applications
- **Population Genetics**: Use biallelic SNP set for structure analysis
- **GWAS**: Apply additional filters for association studies
- **Phylogenetics**: Focus on high-confidence, low-missing data variants
- **Breeding**: Include structural variants and indels for comprehensive analysis

---

## Troubleshooting Common Issues

### Memory Problems
- **Symptoms**: Out of memory errors, job termination
- **Solutions**: Increase JVM heap size, reduce sample batch size, use more efficient tools

### Low Variant Quality
- **Symptoms**: Many variants filtered out, low Ti/Tv ratios
- **Solutions**: Check alignment quality, adjust filtering thresholds, verify reference genome

### Missing Data
- **Symptoms**: High no-call rates, uneven sample representation
- **Solutions**: Check coverage distribution, verify sample quality, adjust calling thresholds

### Performance Issues
- **Symptoms**: Very slow processing, high I/O wait times
- **Solutions**: Use parallel processing, optimize disk access, consider hardware upgrades

This comprehensive description provides the foundation for understanding, troubleshooting, and optimizing the variant calling pipeline for tetraploid Brassica genomics research.
