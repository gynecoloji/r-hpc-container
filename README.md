# Bioinformatics R 4.5.1 Container

A production-ready Singularity/Apptainer container with R 4.5.1 and ~560 bioinformatics packages for High Performance Computing (HPC) systems.

[![Docker](https://img.shields.io/badge/docker-gynecoloji%2Fbioinformatic__r__4__5__1-blue)](https://hub.docker.com/r/gynecoloji/bioinformatic_r_4_5_1)
[![R Version](https://img.shields.io/badge/R-4.5.1-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

---

## 🎯 Overview

This container provides a complete, reproducible R environment for bioinformatics analysis on HPC systems. It includes packages for:

- **Epigenetics**: ChIP-seq, ATAC-seq, Cut&Run analysis (DiffBind, ChIPseeker, chromVAR)
- **Single-cell RNA-seq**: Seurat, monocle3, SingleR, scDblFinder
- **Bulk RNA-seq**: DESeq2, edgeR, limma
- **Genomics**: GenomicRanges, BSgenome, rtracklayer
- **Machine Learning**: mlr3 ecosystem, xgboost, caret
- **Visualization**: ggplot2, ComplexHeatmap, ggtree, plotly

**Total packages**: ~560 with all dependencies locked to specific versions for reproducibility.

---

## 📦 Quick Start

### Pull from DockerHub
```bash
singularity pull docker://gynecoloji/bioinformatic_r_4_5_1:v1
```

### Run an R script
```bash
singularity exec bioinformatic_r_4_5_1_v1.sif Rscript my_analysis.R
```

### Interactive R session
```bash
singularity shell bioinformatic_r_4_5_1_v1.sif
R
```

---

## 📚 Documentation

- **[Installation Guide](docs/INSTALLATION.md)** - Detailed installation instructions
- **[Usage Examples](examples/)** - Sample scripts and HPC job submission examples
- **[Package List](data/installed_packages.csv)** - Complete list of installed packages with versions

---

## 🚀 Features

✅ **Reproducible**: All packages locked to specific versions  
✅ **Portable**: Works across different HPC systems  
✅ **Complete**: Includes ~560 bioinformatics packages  
✅ **Ready-to-use**: No compilation or installation needed  
✅ **Documented**: Full package inventory with versions  

---

## 📋 Requirements

- Singularity/Apptainer 3.0+ (recommended: 3.8+)
- ~2.5 GB storage space
- HPC system or Linux environment

---

## 💡 Usage Examples

### Basic Analysis
```bash
# Run a script with your data mounted
singularity exec --bind /data:/data \
  bioinformatic_r_4_5_1_v1.sif \
  Rscript /data/analysis.R
```

### HPC Job Submission
```bash
# See examples/run_on_hpc.sh for SLURM example
sbatch examples/run_on_hpc.sh
```

### Check Installed Packages
```bash
singularity exec bioinformatic_r_4_5_1_v1.sif \
  Rscript -e "installed.packages()[,c('Package','Version')]"
```

---

## 🏗️ Building from Source

If you want to customize the container:
```bash
# Clone the repository
git clone https://github.com/gynecoloji/apptainer.git
cd apptainer

# Build using the provided script
bash scripts/build.sh
```

See `scripts/build.sh` for the complete build process.

---

## 📊 Included Package Categories

| Category | Key Packages | Count |
|----------|-------------|-------|
| Epigenetics | DiffBind, ChIPseeker, chromVAR, MACS2 | 20+ |
| Single-cell | Seurat, monocle3, scran, SingleR | 30+ |
| RNA-seq | DESeq2, edgeR, limma, tximport | 15+ |
| Genomics | GenomicRanges, BSgenome, Biostrings | 40+ |
| ML/Stats | mlr3*, xgboost, caret, VGAM | 50+ |
| Visualization | ggplot2, ComplexHeatmap, plotly | 30+ |

*Full list in [installed_packages.csv](data/installed_packages.csv)*

---

## 🔧 Troubleshooting

### Common Issues

**Container not found:**
```bash
# Make sure you're in the right directory
ls -lh *.sif
```

**Permission denied:**
```bash
# No sudo needed with Singularity - run as regular user
singularity exec bioinformatic_r_4_5_1_v1.sif R
```

**Package not available:**
```bash
# Check if package is included
singularity exec bioinformatic_r_4_5_1_v1.sif \
  Rscript -e "packageVersion('YourPackage')"
```

See [INSTALLATION.md](docs/INSTALLATION.md) for more troubleshooting.

---

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

For major changes, please open an issue first.

---

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 📞 Contact & Support

- **Issues**: [GitHub Issues](https://github.com/gynecoloji/apptainer/issues)
- **DockerHub**: [gynecoloji/bioinformatic_r_4_5_1](https://hub.docker.com/r/gynecoloji/bioinformatic_r_4_5_1)
- **Email**: gynecoloji@gmail.com

---

## 🙏 Acknowledgments

- Based on [rocker/r-ver:4.5.1](https://hub.docker.com/r/rocker/r-ver)
- Built for reproducible bioinformatics research
- Package selections curated for epigenetics and single-cell analysis

---

## 📅 Version History

- **v1** (2025-12) - Initial release with R 4.5.1 and ~560 packages

---

**Last Updated**: December 2025