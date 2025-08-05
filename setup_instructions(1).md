# Cancer Regulatory Elements Pipeline - Setup & Usage

## Installation

### 1. Install Snakemake and Conda/Mamba

```bash
# Install Miniconda/Mamba
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh

# Install Snakemake
mamba create -n snakemake snakemake pandas pyyaml
conda activate snakemake
```

### 2. Clone and Setup Pipeline

```bash
# Create project directory
mkdir cancer_pipeline && cd cancer_pipeline

# Create directory structure
mkdir -p {data,scripts,envs,results,logs}
mkdir -p data/{genome,chip_seq,expression,conservation,mutations}
mkdir -p data/depmap

# Download the pipeline files (Snakefile, config.yaml, scripts, envs)
# Place them in appropriate directories
```

### 3. Download Reference Data

```bash
# Download genome reference (choose hg19 or hg38)
cd data/genome

# For hg19
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
samtools faidx hg19.fa

# Chromosome sizes
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

# For hg38 
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
```

### 4. Download DepMap Data (Optional)

```bash
cd data/depmap
# Download latest CERES gene effect data
wget https://ndownloader.figshare.com/files/34008503 -O CERES_gene_effect.csv
```

## Configuration

### 1. Edit config.yaml for Your Analysis

```yaml
# Example: Melanoma SOX10 analysis
target_gene: "SOX10"
cancer_type: "melanoma"
genome_build: "hg19"
cell_lines:
  - "A375"
  - "SKMEL5" 
  - "SKMEL28"
  - "MALME3M"
  - "WM115"

# Example: Breast cancer ESR1 analysis  
# target_gene: "ESR1"
# cancer_type: "breast"
# cell_lines:
#   - "MCF7"
#   - "T47D"
#   - "MDA-MB-231"
```

### 2. Adjust Resource Requirements

```yaml
# In config.yaml
max_threads: 8        # Adjust based on your system
max_memory_gb: 32     # Adjust based on available RAM
```

## Usage

### 1. Dry Run (Test Configuration)

```bash
# Activate snakemake environment
conda activate snakemake

# Test the workflow
snakemake --configfile config.yaml --dry-run --cores 1
```

### 2. Run Full Pipeline

```bash
# Run with automatic environment creation
snakemake --configfile config.yaml --cores 8 --use-conda

# Run with specific conda frontend
snakemake --configfile config.yaml --cores 8 --use-conda --conda-frontend mamba
```

### 3. Run Specific Stages

```bash
# Run only ChIP-seq analysis
snakemake --configfile config.yaml --cores 4 \
  results/melanoma_SOX10/chip_seq/consensus_peaks.bed

# Run up to conservation analysis
snakemake --configfile config.yaml --cores 4 \
  results/melanoma_SOX10/conservation/high_conservation_regions.bed
```

### 4. Generate Only the Report

```bash
# If you have intermediate results and just want the report
snakemake --configfile config.yaml --cores 2 \
  results/melanoma_SOX10/final_report.html
```

## Advanced Usage

### 1. Run on Cluster (SLURM Example)

Create `cluster.yaml`:
```yaml
__default__:
  partition: "general"
  time: "04:00:00"
  mem: "8G"
  cores: 1

call_chip_peaks:
  time: "02:00:00" 
  mem: "16G"
  cores: 4

motif_scanning:
  time: "01:00:00"
  mem: "32G"
  cores: 8
```

Run on cluster:
```bash
snakemake --configfile config.yaml \
  --cluster-config cluster.yaml \
  --cluster "sbatch --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.cores}" \
  --jobs 10
```

### 2. Resume Failed Run

```bash
# Snakemake automatically resumes from where it left off
snakemake --configfile config.yaml --cores 8 --use-conda
```

### 3. Clean Up

```bash
# Clean intermediate files
snakemake clean --configfile config.yaml

# Clean everything
snakemake clean_all --configfile config.yaml
```

## Output Structure

```
results/
└── melanoma_SOX10/
    ├── genomic_region.bed                 # Analysis region coordinates
    ├── cell_lines_metadata.tsv           # Cell line information
    ├── chip_seq/
    │   ├── A375_H3K27ac_peaks.bed
    │   ├── consensus_peaks.bed            # Peaks across cell lines
    │   └── ...
    ├── conservation/
    │   ├── high_conservation_regions.bed  # Evolutionarily conserved regions
    │   └── conservation_scores.bigwig
    ├── motifs/
    │   ├── tf_binding_sites.bed           # Predicted TF binding sites
    │   └── motif_enrichment.tsv
    ├── expression/
    │   └── correlation_analysis.tsv       # Expression correlations
    ├── clinical/
    │   ├── mutation_analysis.tsv          # Clinical mutations
    │   └── therapeutic_targets.tsv        # Ranked drug targets
    └── final_report.html                  # Comprehensive report
```

## Troubleshooting

### 1. Common Issues

**Environment Creation Fails:**
```bash
# Update conda/mamba
mamba update -n base mamba
# Try with conda instead of mamba
snakemake --use-conda --conda-frontend conda
```

**Memory Issues:**
```bash
# Reduce memory requirements in config.yaml
max_memory_gb: 16
# Run with fewer cores
snakemake --cores 4
```

**Download Failures:**
```bash
# Check internet connection and retry
# Some files may need manual download
```

### 2. Debugging

```bash
# Verbose output
snakemake --configfile config.yaml --cores 4 --verbose

# Keep temporary files for debugging
snakemake --configfile config.yaml --cores 4 --notemp

# Print shell commands
snakemake --configfile config.yaml --cores 4 --printshellcmds
```

## Customization

### 1. Add New Cancer Types

Edit `scripts/04_identify_cancer_tfs.py`:
```python
cancer_tf_db = {
    'your_cancer': ['TF1', 'TF2', 'TF3'],
    # ... existing entries
}
```

### 2. Add New Analysis Steps

Add new rules to the Snakefile:
```python
rule your_analysis:
    input: "previous_output.txt"
    output: "your_output.txt" 
    script: "scripts/your_analysis.py"

# Update rule all to include your output
```

### 3. Modify Thresholds

Edit `config.yaml`:
```yaml
conservation_threshold: 0.8    # Stricter conservation
peak_qvalue: 0.01             # Stricter peak calling
motif_score_threshold: 0.9    # Stricter motif matching
```

## Example Analyses

### 1. Melanoma SOX10 (Original Study Replication)
```bash
# Use default config.yaml settings
snakemake --configfile config.yaml --cores 8 --use-conda
```

### 2. Breast Cancer ESR1
```bash
# Modify config.yaml for breast cancer
snakemake --configfile config_breast.yaml --cores 8 --use-conda
```

### 3. Multi-Cancer Comparison
```bash
# Run pipeline for multiple cancer types
for cancer in melanoma breast lung; do
    snakemake --configfile config_${cancer}.yaml --cores 4 --use-conda
done
```

## Citation

If you use this pipeline, please cite the original SOX10 study and relevant software:

- **Original Study**: [Your original SOX10 paper]
- **Snakemake**: Köster, Johannes and Rahmann, Sven. "Snakemake - A scalable bioinformatics workflow engine". Bioinformatics 2012.
- **MACS2**: Zhang et al. "Model-based analysis of ChIP-Seq (MACS)". Genome Biology 2008.
- **Bioconductor packages**: [Relevant citations for each package used]