# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GenomicClass_on_Cloud is a collection of educational Google Colab notebooks for learning genomic analysis and bioinformatics. The primary runtime environment is **Google Colab**, not local development.

## Repository Structure

- **Root notebooks** (`.ipynb`): DNA repair mechanisms, NGS transcriptome analysis, NGS variant calling
- **GWAS/**: Genome-wide association studies tutorial
- **ML in genomics/**: Machine learning for cancer classification (logistic regression, decision trees, random forest)
- **Genome_analyses_module/Variant_calling/**: GATK variant calling pipeline for tetraploid species (Brassica) with detailed server execution scripts
- **Metagenomics/**: Project setup guide for metagenomic analysis workflows
- **Web_scraping/**: Biopython and web scraping tutorials

## Development Notes

### Notebooks
- All notebooks are designed for Google Colab with free GPU/TPU access
- Notebooks include installation cells for bioinformatics tools
- Follow the cell execution order - many cells depend on prior cells

### Bioinformatics Tools Referenced
- **QC/Trimming**: fastp
- **Alignment**: bwa, samtools
- **Variant Calling**: GATK (HaplotypeCaller, CombineGVCFs, GenotypeGVCFs)
- **Statistics**: bcftools
- **Reference Download**: NCBI datasets CLI

### Variant Calling Pipeline (Server)
The `Genome_analyses_module/Variant_calling/cheet_sheet.md` contains a complete GATK pipeline for tetraploid species:
1. Quality control with fastp
2. Reference genome indexing (bwa, samtools faidx, GATK CreateSequenceDictionary)
3. Read alignment with bwa-mem
4. GVCF generation with HaplotypeCaller (`--sample-ploidy 4`)
5. Joint calling with CombineGVCFs/GenotypeGVCFs
6. Hard filtering with standard GATK thresholds

All server commands use `nohup` for background execution with logging.

### Metagenomics Project Structure
Standard directory convention in `Metagenomics/readme.md`:
- `00_raw_data/` through `11_final_results/` numbered directories
- `logs/`, `scripts/`, `databases/`, `metadata/` supporting directories
