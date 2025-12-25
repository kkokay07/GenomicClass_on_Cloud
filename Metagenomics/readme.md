# Step 1: Create Project Directory Structure

## Purpose
Create an organized folder structure to keep all your metagenomic analysis files organized. This makes it easier to track data, find results, and reproduce your analysis.

## What You'll Learn
- How to organize a bioinformatics project
- Basic bash commands for directory creation
- Setting up environment variables

---

## Commands

### 1.1 Set Your Project Name and Location

```bash
# Define your project name
PROJECT_NAME="gut_metagenome_analysis"

# Define where you want to store your project (change this to your preference)
BASE_DIR="${HOME}/metagenomics_projects"

# Create the full project path
PROJECT_DIR="${BASE_DIR}/${PROJECT_NAME}"

# Display the project directory path
echo "Project will be created at: ${PROJECT_DIR}"
```

**What this does:**
- Sets a variable for your project name
- Defines where all projects will be stored
- Combines them to create the full path
- Shows you where the project will be created

---

### 1.2 Create Main Project Directory

```bash
# Create the main project directory
mkdir -p "${PROJECT_DIR}"

# Navigate into the project directory
cd "${PROJECT_DIR}"

# Confirm you're in the right place
pwd
```

**What this does:**
- `mkdir -p` creates the directory (and parent directories if needed)
- `cd` changes to the new directory
- `pwd` shows your current location

---

### 1.3 Create Subdirectories for Each Analysis Step

```bash
# Create all subdirectories at once
mkdir -p \
  00_raw_data \
  01_qc_reports \
  02_trimmed_reads \
  03_host_removed \
  04_assembly \
  05_mapping \
  06_binning \
  07_bin_quality \
  08_taxonomy \
  09_annotation \
  10_functional_analysis \
  11_final_results \
  logs \
  scripts \
  databases \
  metadata
```

**What this does:**
- Creates all necessary folders in one command
- `\` allows the command to continue on the next line
- Numbers help keep folders in order

---

### 1.4 Create Subdirectories Within Each Main Folder

```bash
# QC subdirectories
mkdir -p 01_qc_reports/{raw,trimmed,multiqc}

# Assembly subdirectories
mkdir -p 04_assembly/{megahit,metaspades,stats}

# Binning subdirectories
mkdir -p 06_binning/{metabat2,maxbin2,concoct,dastool,all_bins}

# Taxonomy subdirectories
mkdir -p 08_taxonomy/{kraken2,gtdbtk,visualization}

# Annotation subdirectories
mkdir -p 09_annotation/{prokka,prodigal,dram}

# Logs subdirectories (organized by date)
mkdir -p logs/{$(date +%Y-%m-%d)}
```

**What this does:**
- `{}` creates multiple subdirectories within a parent folder
- Organizes outputs from different tools
- Creates a date-stamped log folder

---

### 1.5 View Your Directory Structure

```bash
# Install tree command if not available (optional)
# sudo apt-get install tree  # Ubuntu/Debian
# brew install tree          # macOS

# View the directory tree (if tree is installed)
tree -L 2

# Or use ls to view directories
ls -lh
```

**What this does:**
- `tree -L 2` shows directory structure up to 2 levels deep
- `ls -lh` lists directories with human-readable sizes

---

### 1.6 Create a README File

```bash
# Create a README file for your project
cat > README.txt << 'EOF'
PROJECT: Gut Metagenome Analysis
DATE STARTED: $(date +%Y-%m-%d)
AUTHOR: [Your Name]

DESCRIPTION:
This project contains metagenomic analysis of gut microbiome samples.

DIRECTORY STRUCTURE:
00_raw_data/          - Raw FASTQ files from sequencing
01_qc_reports/        - Quality control reports (FastQC, MultiQC)
02_trimmed_reads/     - Quality-filtered and trimmed reads
03_host_removed/      - Reads after host contamination removal
04_assembly/          - Assembled contigs
05_mapping/           - Read mapping results
06_binning/           - Genome bins (MAGs)
07_bin_quality/       - Bin quality assessment results
08_taxonomy/          - Taxonomic classification results
09_annotation/        - Gene prediction and annotation
10_functional_analysis/ - Functional profiling results
11_final_results/     - Final outputs and figures
logs/                 - All log files
scripts/              - Analysis scripts
databases/            - Reference databases
metadata/             - Sample metadata and documentation

SAMPLE INFORMATION:
Sample ID: [To be added]
SRA Accession: [To be added]
Source: [To be added]
EOF

# View the README
cat README.txt
```

**What this does:**
- Creates a README file documenting your project
- `cat > file << 'EOF'` creates a multi-line file
- Documents the purpose of each directory

---

### 1.7 Create a Sample Tracking File

```bash
# Create a CSV file to track samples
cat > metadata/samples.csv << 'EOF'
sample_id,sra_accession,description,source,date_downloaded
sample_01,SRR1234567,Healthy gut microbiome,NCBI_SRA,
sample_02,SRR1234568,Disease gut microbiome,NCBI_SRA,
EOF

# View the file
cat metadata/samples.csv
```

**What this does:**
- Creates a CSV file to track sample information
- Helps organize multiple samples
- Can be updated as you download data

---

### 1.8 Create a Configuration File

```bash
# Create a configuration file for your analysis parameters
cat > config.txt << 'EOF'
# Metagenomic Analysis Configuration
# Edit these parameters for your analysis

# Computational Resources
THREADS=8
MEMORY=32G

# Quality Control Parameters
MIN_QUALITY=20
MIN_LENGTH=50

# Assembly Parameters
ASSEMBLER=megahit
MIN_CONTIG_LENGTH=1000

# Binning Parameters
MIN_BIN_SIZE=200000
MIN_COMPLETENESS=50
MAX_CONTAMINATION=10

# Database Paths (update these)
KRAKEN2_DB=/path/to/kraken2/database
GTDBTK_DB=/path/to/gtdbtk/database
CHECKM2_DB=/path/to/checkm2/database
EOF

# View the configuration
cat config.txt
```

**What this does:**
- Creates a central configuration file
- Stores all parameters in one place
- Makes it easy to reproduce analysis

---

### 1.9 Create an Analysis Log

```bash
# Create a log file to track your analysis steps
cat > logs/analysis_log.txt << EOF
METAGENOMIC ANALYSIS LOG
========================
Project: ${PROJECT_NAME}
Started: $(date)
User: $(whoami)
Working Directory: $(pwd)

ANALYSIS STEPS:
--------------
EOF

# View the log
cat logs/analysis_log.txt
```

**What this does:**
- Creates a log file to document your work
- Records when analysis started
- Helps track what you've done

---

### 1.10 Verify Directory Structure

```bash
# Count directories created
echo "Number of directories created:"
find . -type d | wc -l

# List main directories
echo -e "\nMain directories:"
ls -d */

# Check disk space available
echo -e "\nDisk space:"
df -h .

# Show current directory structure
echo -e "\nDirectory structure:"
ls -R | head -50
```

**What this does:**
- Counts all directories created
- Lists main directories
- Checks available disk space
- Shows the structure you've created

---

## Expected Output Structure

After completing Step 1, your directory should look like this:

```
gut_metagenome_analysis/
├── 00_raw_data/
├── 01_qc_reports/
│   ├── raw/
│   ├── trimmed/
│   └── multiqc/
├── 02_trimmed_reads/
├── 03_host_removed/
├── 04_assembly/
│   ├── megahit/
│   ├── metaspades/
│   └── stats/
├── 05_mapping/
├── 06_binning/
│   ├── metabat2/
│   ├── maxbin2/
│   ├── concoct/
│   ├── dastool/
│   └── all_bins/
├── 07_bin_quality/
├── 08_taxonomy/
│   ├── kraken2/
│   ├── gtdbtk/
│   └── visualization/
├── 09_annotation/
│   ├── prokka/
│   ├── prodigal/
│   └── dram/
├── 10_functional_analysis/
├── 11_final_results/
├── logs/
│   ├── 2024-12-25/
│   └── analysis_log.txt
├── scripts/
├── databases/
├── metadata/
│   └── samples.csv
├── README.txt
└── config.txt
```

---

## Check Your Work

Run these commands to verify Step 1 is complete:

```bash
# 1. Verify you're in the project directory
echo "Current directory: $(pwd)"

# 2. Check that all main directories exist
for dir in 00_raw_data 01_qc_reports 02_trimmed_reads 03_host_removed \
           04_assembly 05_mapping 06_binning 07_bin_quality \
           08_taxonomy 09_annotation 10_functional_analysis \
           11_final_results logs scripts databases metadata; do
    if [ -d "$dir" ]; then
        echo "✓ $dir exists"
    else
        echo "✗ $dir missing"
    fi
done

# 3. Check that README exists
if [ -f "README.txt" ]; then
    echo "✓ README.txt exists"
else
    echo "✗ README.txt missing"
fi

# 4. Check that config exists
if [ -f "config.txt" ]; then
    echo "✓ config.txt exists"
else
    echo "✗ config.txt missing"
fi
```

---

## Common Issues and Solutions

**Issue 1: Permission denied**
```bash
# Solution: Check directory permissions
ls -la ${BASE_DIR}

# Or create in a different location where you have permissions
BASE_DIR="${HOME}/projects"
```

**Issue 2: Directory already exists**
```bash
# Solution: Use a different project name or remove old directory
PROJECT_NAME="gut_metagenome_analysis_v2"
```

**Issue 3: Not enough disk space**
```bash
# Check available space
df -h

# Solution: Use a different location with more space
BASE_DIR="/mnt/large_disk/projects"
```

---

## Tips for Students

1. **Use meaningful names** - Choose project names that describe your data
2. **Stay organized** - Always put files in the correct directories
3. **Document everything** - Update README and logs as you work
4. **Check your location** - Use `pwd` to confirm you're in the right directory
5. **Back up important data** - Copy metadata and config files regularly

---

## Save This Setup for Future Use

```bash
# Create a script to recreate this structure
cat > scripts/create_project_structure.sh << 'EOF'
#!/bin/bash
# Script to create metagenomic project structure

PROJECT_NAME=$1
BASE_DIR="${HOME}/metagenomics_projects"
PROJECT_DIR="${BASE_DIR}/${PROJECT_NAME}"

mkdir -p "${PROJECT_DIR}"/{00_raw_data,01_qc_reports/{raw,trimmed,multiqc},02_trimmed_reads,03_host_removed,04_assembly/{megahit,metaspades,stats},05_mapping,06_binning/{metabat2,maxbin2,concoct,dastool,all_bins},07_bin_quality,08_taxonomy/{kraken2,gtdbtk,visualization},09_annotation/{prokka,prodigal,dram},10_functional_analysis,11_final_results,logs,scripts,databases,metadata}

echo "Project structure created at: ${PROJECT_DIR}"
EOF

# Make it executable
chmod +x scripts/create_project_structure.sh

# Test it (optional - for future projects)
# ./scripts/create_project_structure.sh new_project_name
```

---

## What's Next?

✅ **Step 1 Complete!** You now have an organized project structure.

**Next Step:** [Step 2 - Install Required Software →]

In Step 2, you'll learn how to install all the bioinformatics tools needed for the analysis.

---

## Summary of Commands Used

| Command | Purpose |
|---------|---------|
| `mkdir -p` | Create directories (including parents) |
| `cd` | Change directory |
| `pwd` | Print working directory |
| `cat > file << EOF` | Create a file with content |
| `tree` | Display directory structure |
| `ls -lh` | List files with details |
| `find` | Find files/directories |
| `df -h` | Check disk space |
| `chmod +x` | Make file executable |

---

**Time Required:** 5-10 minutes  
**Difficulty:** ⭐ Beginner  
**Prerequisites:** None
