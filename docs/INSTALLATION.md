# Installation Guide

This guide shows you how to install and use the Bioinformatics R 4.5.1 container on HPC systems.

---

## 📋 Prerequisites

Before installing, ensure you have:
- Access to an HPC system with Singularity/Apptainer installed
- Sufficient storage space (~2.5 GB for the SIF file)
- Internet connection (for pulling from DockerHub) OR ability to transfer files

Check your Singularity version:
```bash
singularity --version
# or
apptainer --version
```

**Minimum version:** 1.4+ (recommended: 1.4+)

---

## 🚀 Installation Methods

### **Method 1: Pull from DockerHub (Recommended)**

If your HPC has internet access:
```bash
# Pull and convert Docker image to SIF
singularity pull docker://gynecoloji/bioinformatic_r_4_5_1:v1

# This creates: bioinformatic_r_4_5_1_v1.sif
```

**Rename for cleaner naming (optional):**
```bash
mv bioinformatic_r_4_5_1_v1.sif bioinformatic_r_4.5.1.sif
```

---

### **Method 2: Transfer Pre-built SIF**

If your HPC has no internet access:

**On your local machine:**
```bash
# Pull the image locally first
singularity pull docker://gynecoloji/bioinformatic_r_4_5_1:v1

# Transfer to HPC via scp
scp bioinformatic_r_4_5_1_v1.sif username@hpc-server:/path/to/destination/
```

**On HPC:**
```bash
# Verify the transfer
ls -lh bioinformatic_r_4_5_1_v1.sif
md5sum bioinformatic_r_4_5_1_v1.sif  # Optional: verify integrity
```

---

## ✅ Verify Installation

Test that the container works:
```bash
# Test R installation
singularity exec bioinformatic_r_4_5_1_v1.sif R --version

# Check installed packages
singularity exec bioinformatic_r_4_5_1_v1.sif Rscript -e "cat(length(.packages(all.available=TRUE)), 'packages installed\n')"

# Interactive R session
singularity shell bioinformatic_r_4_5_1_v1.sif
# Then in the container:
R
```

**Expected output:** R version 4.5.1 with ~560 packages installed

---

## 📦 What's Included

This container includes:
- **Base:** R 4.5.1 (rocker/r-ver:4.5.1)
- **Packages:** ~560 bioinformatics packages including:
  - Epigenetics analysis tools (ChIP-seq, ATAC-seq, Cut&Run-seq)
  - RNA-seq analysis (bulk and single-cell)
  - Common dependencies and utilities

See `data/installed_packages.csv` for the complete package list with versions.

---

## 🎯 Quick Start

Once installed, run your analysis:
```bash
# Run an R script
singularity exec bioinformatic_r_4_5_1_v1.sif Rscript my_analysis.R

# Interactive session
singularity shell bioinformatic_r_4_5_1_v1.sif

# Bind mount your data directory
singularity exec --bind /path/to/data:/data bioinformatic_r_4_5_1_v1.sif Rscript /data/analysis.R
```

See `examples/` directory for more usage examples.

---

## 🐛 Troubleshooting

### **"command not found: singularity"**
- Your HPC may use `apptainer` instead - replace `singularity` with `apptainer`
- Load the module: `module load singularity` or `module load apptainer`

### **Pull fails with network error**
- Check internet connectivity: `ping google.com`
- Some HPCs restrict outbound connections - use Method 2 (transfer pre-built SIF)

### **"manifest unknown" error**
- Verify the tag exists: https://hub.docker.com/r/gynecoloji/bioinformatic_r_4_5_1/tags
- Current available tag: `v1`

### **Package not found when running R**
- Check if the package is included in the image

---

## 📞 Support

- **Issues:** [GitHub Issues](https://github.com/gynecoloji/apptainer/issues)
- **DockerHub:** https://hub.docker.com/r/gynecoloji/bioinformatic_r_4_5_1

---

## 🔄 Available Versions

Check DockerHub for available versions: https://hub.docker.com/r/gynecoloji/bioinformatic_r_4_5_1/tags

**Current version:** `v1`

To pull a specific version:
```bash
singularity pull docker://gynecoloji/bioinformatic_r_4_5_1:v1
```