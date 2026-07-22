#!/bin/bash
# BASALT 1.2 Installation Script

set -euo pipefail
 
# 1. Ensure the Conda environment is activated
if [ -z "${CONDA_PREFIX:-}" ]; then
    echo "Error: Please activate your Conda environment first."
    exit 1
fi

if [ ! -d "$CONDA_PREFIX/bin" ] || [ ! -w "$CONDA_PREFIX/bin" ]; then
    echo "Error: $CONDA_PREFIX/bin is not a writable directory."
    echo "Use a user-owned Conda/Micromamba environment; do not use chmod 777."
    exit 1
fi

if [ ! -x "$CONDA_PREFIX/bin/jgi_summarize_bam_contig_depths" ]; then
    echo "Error: metabat2 did not provide jgi_summarize_bam_contig_depths."
    echo "Recreate or repair the environment from basalt_environment.yml."
    exit 1
fi

# 2. Copy source files directly into the bin directory. The legacy bundled
# jgi binary is deliberately skipped: the maintained metabat2 package supplies
# the platform-compatible executable. Python and Perl files receive mode 0755
# (owner rwx, group/others r-x); never use world-writable mode 0777 here.
echo "Copying files to $CONDA_PREFIX/bin/ ..."
for source in BASALT/*; do
    case "$(basename "$source")" in
        BASALT|jgi_summarize_bam_contig_depths|__pycache__|*.pyc|*.pyo)
            continue
            ;;
    esac
    target="$CONDA_PREFIX/bin/$(basename "$source")"
    cp -R "$source" "$target"
    case "$source" in
        *.py|*.pl)
            chmod 0755 "$target"
            ;;
    esac
done

# 3. Create the 'BASALT' launcher command
# This allows you to type 'BASALT' in the terminal to execute BASALT.py located in bin
printf '%s\n' \
    '#!/bin/bash' \
    'set -euo pipefail' \
    'bin_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"' \
    'exec "$bin_dir/python" "$bin_dir/BASALT.py" "$@"' \
    > "$CONDA_PREFIX/bin/BASALT"
chmod 0755 "$CONDA_PREFIX/bin/BASALT"

echo "Installation complete. Files have been placed directly in $CONDA_PREFIX/bin/"
