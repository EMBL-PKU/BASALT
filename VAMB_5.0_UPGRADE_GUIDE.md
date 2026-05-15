# VAMB 5.0 Upgrade Guide

## Overview

VAMB 5.0 brings Python 3.12 compatibility and improved performance to BASALT. This guide helps you upgrade from older VAMB versions or install VAMB 5.0 for the first time.

## What's New in VAMB 5.0

- ✅ **Python 3.12 Support**: Full compatibility with the latest Python version
- ✅ **Improved Performance**: Faster binning and reduced memory usage
- ✅ **Better Error Messages**: More informative error reporting
- ✅ **Updated Command Syntax**: New `vamb bin default` command structure
- ✅ **Enhanced Multi-threading**: Better thread control with `-p` parameter

## Installation Options

### Option 1: Fresh BASALT Installation (Recommended)

If you're setting up BASALT for the first time or want a clean installation:

```bash
# Clone or download BASALT
cd BASALT

# Create environment with VAMB 5.0 included
conda env create -f basalt_new_environment.yml

# Activate environment
conda activate basalt

# Verify VAMB installation
vamb --version
# Should show: Vamb 5.x.x
```

### Option 2: Upgrade Existing BASALT Environment

If you already have BASALT installed and want to add/upgrade VAMB:

```bash
# Activate your existing BASALT environment
conda activate basalt_env  # or your environment name

# Remove old VAMB if present
pip uninstall vamb -y

# Install VAMB 5.0 (GPU version)
pip install "vamb>=5.0.0"

# Verify installation
vamb --version
```

### Option 3: CPU-Only Installation

For systems without GPU or with CUDA compatibility issues:

```bash
# Activate BASALT environment
conda activate basalt_env

# Install VAMB 5.0 without dependencies
pip install "vamb>=5.0.0" --no-deps

# Install PyTorch CPU-only version
conda install pytorch cpuonly -c pytorch

# Verify installation
vamb --version
```

## Verifying VAMB 5.0 Installation

Run these commands to verify your installation:

```bash
# Check VAMB version
vamb --version

# Check if VAMB can access help
vamb bin default --help

# Test basic functionality (optional)
# This will show available subcommands
vamb --help
```

Expected output for VAMB 5.0:
```
Vamb 5.x.x
```

## Using VAMB 5.0 with BASALT

No changes needed! BASALT automatically detects VAMB 5.0 and uses the correct command syntax:

```bash
# Standard BASALT command with VAMB
BASALT -a assembly.fa -s reads_R1.fq,reads_R2.fq -e v -t 60 -m 200 --mode new -o output

# Combine with other binners
BASALT -a assembly.fa -s reads_R1.fq,reads_R2.fq -e m,v,l -t 60 -m 200 --mode new -o output
```

## Troubleshooting

### Issue: "vamb: command not found"

**Solution:**
```bash
# Ensure VAMB is installed
pip install "vamb>=5.0.0"

# Check if pip installed to correct environment
which vamb
```

### Issue: "ModuleNotFoundError: No module named 'torch'"

**Solution:**
```bash
# Install PyTorch
conda install pytorch -c pytorch

# Or for CPU-only
conda install pytorch cpuonly -c pytorch
```

### Issue: CUDA errors with GPU version

**Solution:**
```bash
# Switch to CPU-only version
pip uninstall vamb torch torchvision -y
pip install "vamb>=5.0.0" --no-deps
conda install pytorch cpuonly -c pytorch
```

### Issue: "vamb bin default: command not found"

This means you have VAMB 4.x installed. Either:

**Option A: Upgrade to VAMB 5.0**
```bash
pip install --upgrade "vamb>=5.0.0"
```

**Option B: Continue with VAMB 4.x**
BASALT's error handling will catch the command failure and log it. You can continue using other binners.

### Issue: VAMB fails during BASALT run

Check the log file:
```bash
# View error details
cat Basalt_log.txt | grep -A 5 "vamb"

# Common issues:
# - Insufficient memory: Increase -m parameter
# - Missing BAM files: Check read mapping step
# - CUDA errors: Switch to CPU version
```

## Performance Comparison

| Version | Python Support | Speed | Memory Usage | GPU Support |
|---------|---------------|-------|--------------|-------------|
| VAMB 4.x | 3.7-3.10 | Baseline | Baseline | Yes |
| VAMB 5.0 | 3.8-3.12 | ~20% faster | ~15% less | Yes |

## Migration Notes

### Command Changes

**VAMB 4.x:**
```bash
vamb --outdir output --fasta assembly.fa --bamfiles *.bam --minfasta 500000
```

**VAMB 5.0:**
```bash
vamb bin default --outdir output --fasta assembly.fa --bamfiles *.bam --minfasta 500000 -p 60
```

BASALT handles this automatically - no manual changes needed!

### Output Format

VAMB 5.0 maintains backward compatibility with the output format:
- `clusters.tsv` - Same format as VAMB 4.x
- Bin FASTA files - Same naming convention
- Quality metrics - Compatible with CheckM/CheckM2

## Recommended Setup for Different Use Cases

### High-Performance Computing (HPC) with GPU
```bash
# Install GPU version
pip install "vamb>=5.0.0"
conda install pytorch pytorch-cuda=11.8 -c pytorch -c nvidia

# Use with BASALT
BASALT -a assembly.fa -s reads.fq -e v -t 80 -m 500 --mode new -o output
```

### Local Workstation (CPU-only)
```bash
# Install CPU version
pip install "vamb>=5.0.0" --no-deps
conda install pytorch cpuonly -c pytorch

# Use with BASALT
BASALT -a assembly.fa -s reads.fq -e v -t 16 -m 64 --mode new -o output
```

### Cloud Computing (AWS/GCP/Azure)
```bash
# Check GPU availability first
nvidia-smi

# If GPU available, use GPU version
pip install "vamb>=5.0.0"

# If no GPU, use CPU version
pip install "vamb>=5.0.0" --no-deps
conda install pytorch cpuonly -c pytorch
```

## Getting Help

If you encounter issues with VAMB 5.0:

1. **Check BASALT logs**: `cat Basalt_log.txt`
2. **Check VAMB version**: `vamb --version`
3. **Test VAMB directly**: `vamb bin default --help`
4. **Report issues**: https://github.com/PKU-EMBL/BASALT/issues

Include in your report:
- VAMB version (`vamb --version`)
- Python version (`python --version`)
- BASALT command used
- Error messages from `Basalt_log.txt`
- System info (OS, RAM, GPU if applicable)

## Additional Resources

- **VAMB Documentation**: https://github.com/RasmussenLab/vamb
- **BASALT Documentation**: https://github.com/PKU-EMBL/BASALT
- **PyTorch Installation**: https://pytorch.org/get-started/locally/
