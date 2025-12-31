#!/bin/bash
#SBATCH --job-name=r_analysis          # Job name
#SBATCH --output=r_analysis_%j.log     # Standard output log (%j = job ID)
#SBATCH --error=r_analysis_%j.err      # Standard error log
#SBATCH --ntasks=1                     # Number of tasks
#SBATCH --cpus-per-task=4              # Number of CPU cores per task
#SBATCH --mem=16G                      # Memory per node
#SBATCH --time=02:00:00                # Time limit (HH:MM:SS)
#SBATCH --partition=general            # Partition name (adjust for your HPC)

# ==============================================================================
# HPC Job Script for Bioinformatics R 4.5.1 Container
# ==============================================================================
#
# This script demonstrates how to run R analyses using the container on HPC
# systems with SLURM scheduler. Modify the SBATCH parameters above according
# to your HPC system's configuration and resource requirements.
#
# Usage:
#   sbatch run_on_hpc.sh
#
# ==============================================================================

# Exit on error
set -e

# Print job information
echo "=========================================="
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"
echo "=========================================="
echo ""

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------

# Container location (modify this path)
CONTAINER="/path/to/bioinformatic_r_4_5_1_v1.sif"

# Data directories (modify these paths)
DATA_DIR="/path/to/your/data"
OUTPUT_DIR="/path/to/your/output"
SCRIPT_DIR="/path/to/your/scripts"

# R script to run
R_SCRIPT="example_analysis.R"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# ------------------------------------------------------------------------------
# LOAD MODULES (if needed)
# ------------------------------------------------------------------------------

# Uncomment and modify if your HPC requires module loading
# module load singularity/3.8.0
# module load gcc/9.3.0  # or other dependencies

echo "Modules loaded:"
module list
echo ""

# ------------------------------------------------------------------------------
# VERIFY CONTAINER
# ------------------------------------------------------------------------------

echo "Verifying container..."
if [ ! -f "$CONTAINER" ]; then
    echo "ERROR: Container not found at $CONTAINER"
    exit 1
fi

echo "Container found: $CONTAINER"
echo "Container size: $(du -h $CONTAINER | cut -f1)"
echo ""

# Test container
echo "Testing container..."
singularity exec "$CONTAINER" R --version
echo ""

# ------------------------------------------------------------------------------
# RUN ANALYSIS
# ------------------------------------------------------------------------------

echo "=========================================="
echo "Starting R analysis..."
echo "=========================================="
echo ""

# Run R script with container
# Bind mount directories to make them accessible inside container
singularity exec \
    --bind "$DATA_DIR:/data" \
    --bind "$OUTPUT_DIR:/output" \
    --bind "$SCRIPT_DIR:/scripts" \
    "$CONTAINER" \
    Rscript "/scripts/$R_SCRIPT"

# Alternative: Run the example script included in the repository
# singularity exec "$CONTAINER" Rscript /workspace/examples/example_analysis.R

# ------------------------------------------------------------------------------
# VERIFY OUTPUT
# ------------------------------------------------------------------------------

echo ""
echo "=========================================="
echo "Checking outputs..."
echo "=========================================="
echo ""

# List output files
echo "Output files created:"
ls -lh "$OUTPUT_DIR/"
echo ""

# ------------------------------------------------------------------------------
# CLEANUP (optional)
# ------------------------------------------------------------------------------

# Uncomment if you want to clean up intermediate files
# echo "Cleaning up temporary files..."
# rm -f /tmp/R_session_*.RData

# ------------------------------------------------------------------------------
# JOB SUMMARY
# ------------------------------------------------------------------------------

echo ""
echo "=========================================="
echo "Job completed at: $(date)"
echo "Total runtime: $SECONDS seconds"
echo "=========================================="

exit 0