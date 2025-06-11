rule download_fastq:
    """
    Download with fasterq-dump and gzip
    """
    output:
        fq1 = temp("{scratch}/{sra}_1.fastq.gz"),
        fq2 = temp("{scratch}/{sra}_2.fastq.gz")
    params:
        scratch = SCRATCH
    log:
        "{scratch}/logs/download_{sra}.log"
    threads: THREADS
    resources:
        mem_mb=2000
    conda:
        "envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.scratch}/logs
        fasterq-dump {wildcards.sra} \
            -O {params.scratch} \
            --split-files --skip-technical \
            &> {log}
        gzip -f {params.scratch}/{wildcards.sra}_1.fastq \
                {params.scratch}/{wildcards.sra}_2.fastq
        mv {params.scratch}/{wildcards.sra}_1.fastq.gz {output.fq1}
        mv {params.scratch}/{wildcards.sra}_2.fastq.gz {output.fq2}
        """