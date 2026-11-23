# seq-processing-snakemake
process sequencing data in accordance with ENCODE best practices


# run pipeline
```
    conda activate run_snakemake9
    snakemake --configfile config/config.yml --profile .snakemake_profile/slurm

```
# reference
The bowtie2 index in `config/config.yml` (`hg38_resources/hg38_bowtie2_index`) was deleted on 2026-11-23 to save space. Rebuild it before running:
```
    bowtie2-build --threads 8 /oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa /oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38_bowtie2_index
```
