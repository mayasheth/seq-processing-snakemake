# seq-processing-snakemake
process sequencing data in accordance with ENCODE best practices


# run pipeline
```
    conda activate run_snakemake9
    snakemake --configfile config/config.yml --profile .snakemake_profile/slurm

```