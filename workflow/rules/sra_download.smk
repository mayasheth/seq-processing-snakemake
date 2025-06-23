rule download_fastq:
    """
    Download with fasterq-dump and gzip
    """
    output:
        fq1 = os.path.join(SCRATCH_DIR, "{sample_name}", "{assay}", "{accession}_1.fastq.gz"),
        fq2 = os.path.join(SCRATCH_DIR, "{sample_name}", "{assay}", "{accession}_2.fastq.gz")
    params:
        scratch = SCRATCH_DIR,
        results = RESULTS_DIR
    log:
        os.path.join(RESULTS_DIR, "logs", "download_{sample_name}_{accession}_{assay}.log")
    threads: THREADS
    resources:
        mem_mb=16000,
        runtime=720
    conda:
        "../envs/seq_tools.yml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.results}/logs
        fasterq-dump {wildcards.accession} \
            -O {params.scratch} \
            --split-files --skip-technical \
            &> {log}
        gzip -f {params.scratch}/{wildcards.accession}_1.fastq \
                {params.scratch}/{wildcards.accession}_2.fastq
        mv {params.scratch}/{wildcards.accession}_1.fastq.gz {output.fq1}
        mv {params.scratch}/{wildcards.accession}_2.fastq.gz {output.fq2}
        """