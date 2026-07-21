# BASALT implementation map

This directory contains the installed BASALT command-line workflow and its internal stage modules. These modules are not a stable Python API.

Primary entry points:

- `BASALT.py`: CLI parsing and workflow dispatch.
- `BASALT_main_d.py`: default CheckM2 orchestration.
- `BASALT_main_c.py` and `BASALT_main_c_*`: legacy CheckM orchestration.
- `Data_feeding.py`: preparation of external binsets.

Numbered `S1`–`S10` files implement candidate generation, bin comparison, contig screening and retrieval, reassembly, and OLC stages. Matching `_checkm.py` files support the legacy quality-control branch.

See the maintained [developer reference](../BASALT_Guide/docs/api/index.md) and [pipeline methods](../BASALT_Guide/docs/pipeline.md) for the complete map, stability boundaries, and documentation build instructions.
