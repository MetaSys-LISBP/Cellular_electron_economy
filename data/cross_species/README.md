# Cross-species model snapshots

This directory contains the exact genome-scale model snapshots and configuration files used for the cross-species ETN analysis.

The paired main comparison includes *S. cerevisiae*, *B. subtilis*, *S. enterica* and *K. phaffii*, together with the primary *E. coli* iML1515 model stored one directory above. iJO1366 is included as an independent *E. coli* reconstruction robustness test.


`model_metadata.tsv` and `checksums.sha256` record provenance and exact bundled-file hashes. `condition_config.tsv` documents model-specific physiological configurations. The shared ETN reconstruction algorithm applies the same redox rules to every organism.
