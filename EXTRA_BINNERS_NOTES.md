# Extra Binners Usage Notes

## Error Handling

BASALT now includes robust error handling for extra binners (MetaBinner, VAMB, LorBin). If an extra binner fails during execution:

- **The pipeline will continue** instead of stopping
- A warning message will be displayed: `[WARNING] Extra binner 'X' failed for assembly. Skipping...`
- Error details will be logged to `Basalt_log.txt`
- Other binners and assemblies will continue to be processed

This allows you to:
1. Run BASALT with multiple extra binners without worrying about one failure stopping the entire pipeline
2. Review the log file to see which binners failed and why
3. Manually run failed binners later and use data feeding to incorporate results

## VAMB 5.0 Support (Python 3.12 Compatible)

### Background

BASALT now supports VAMB 5.0, which is compatible with Python 3.12 and offers improved performance. VAMB can be installed in two modes:
- **GPU mode** (default): Requires CUDA and compatible GPU drivers
- **CPU mode**: Works on any system without GPU requirements

### Installing VAMB 5.0

**Option 1: Install with the BASALT environment (Recommended)**

The `basalt_new_environment.yml` now includes VAMB 5.0 by default:

```bash
# Create new BASALT environment with VAMB 5.0
conda env create -f basalt_new_environment.yml
conda activate basalt
```

**Option 2: Add VAMB 5.0 to existing environment**

If you already have BASALT installed, you can upgrade to VAMB 5.0:

```bash
# Activate your BASALT environment
conda activate basalt_env

# Install VAMB 5.0 (GPU version with CUDA support)
pip install "vamb>=5.0.0"

# OR install VAMB 5.0 CPU-only version
pip install "vamb>=5.0.0" --no-deps
conda install pytorch cpuonly -c pytorch
```

### VAMB 5.0 vs 4.x Compatibility

BASALT automatically detects and uses the appropriate VAMB command syntax:
- **VAMB 5.0**: Uses `vamb bin default` command with explicit thread control
- **VAMB 4.x**: Falls back to legacy `vamb` command (if 5.0 command fails)

The error handling system will catch any compatibility issues and log them to `Basalt_log.txt`.

### Python 3.12 Compatibility

VAMB 5.0 is fully compatible with Python 3.12, which is the default Python version in the new BASALT environment. This resolves previous compatibility issues with older VAMB versions.

### Using VAMB with BASALT

Once VAMB is installed (CPU or GPU version), you can use it with BASALT:

```bash
BASALT -a assembly.fa -s reads_R1.fq,reads_R2.fq -e v -t 60 -m 200 --mode new -o output
```

Where `-e v` enables VAMB as an extra binner.

### Combining Multiple Extra Binners

You can use multiple extra binners together. With the new error handling, if one fails, others will continue:

```bash
# Use both MetaBinner and VAMB
BASALT -a assembly.fa -s reads_R1.fq,reads_R2.fq -e m,v -t 60 -m 200 --mode new -o output

# Use MetaBinner, VAMB, and LorBin
BASALT -a assembly.fa -s reads_R1.fq,reads_R2.fq -e m,v,l -t 60 -m 200 --mode new -o output
```

### Troubleshooting

If VAMB fails:

1. **Check the log file**: `Basalt_log.txt` will contain error details
2. **Verify installation**: Run `vamb --help` to ensure VAMB is properly installed
3. **Check BAM files**: VAMB requires BAM files from read mapping
4. **Memory requirements**: VAMB can be memory-intensive for large assemblies

Common VAMB errors:
- **CUDA errors**: Install CPU-only version (see above)
- **Memory errors**: Increase RAM allocation with `-m` parameter
- **Missing BAM files**: Ensure read mapping completed successfully

### Manual VAMB Execution

If VAMB fails during BASALT execution, you can run it manually and feed the results back:

```bash
# Run VAMB manually
vamb --outdir assembly_vamb --fasta assembly.fa --bamfiles *.bam --minfasta 500000

# Use BASALT data feeding to incorporate VAMB results
BASALT -d assembly_vamb_bins -o final_output --mode continue
```

## Extra Binner Options

- **m**: MetaBinner - Deep learning-based binner
- **v**: VAMB - Variational autoencoder-based binner  
- **l**: LorBin - Long-read oriented binner

## Performance Notes

- **MetaBinner**: Requires k-mer generation, can be slow for large assemblies
- **VAMB**: Fast but memory-intensive, works well with multiple samples
- **LorBin**: Optimized for long-read assemblies

## Recommendations

1. **For standard projects**: Use default binners (MetaBAT2, MaxBin2, CONCOCT)
2. **For complex communities**: Add MetaBinner with `-e m`
3. **For multi-sample projects**: Add VAMB with `-e v`
4. **For long-read assemblies**: Add LorBin with `-e l`
5. **For maximum coverage**: Use all extra binners `-e m,v,l` (with error handling, failures won't stop the pipeline)

## Support

If you continue to experience issues with extra binners, please report them at:
https://github.com/PKU-EMBL/BASALT/issues

Include:
- BASALT command used
- Contents of `Basalt_log.txt`
- VAMB/MetaBinner/LorBin version
- System information (OS, RAM, GPU if applicable)
