ule align:
    """
    Bowtie2 → sorted BAM
    """
    input:
        fq1 = rules.download_fastq.output.fq1,
        fq2 = rules.download_fastq.output.fq2
    output:
        bam = temp("{scratch}/{sra}.bam")
    params:
        index = INDEX
    log:
        "{scratch}/logs/align_{sra}.log"
    threads: THREADS
    resources:
        mem_mb=8000
    conda:
        "envs/pipeline.yaml"
    shell:
        r"""
        bowtie2 -X 2000 --mm \
          -x {params.index} \
          -1 {input.fq1} -2 {input.fq2} \
          -p {threads} 2> {log} | \
        samtools view -h -@ {threads} - | \
        samtools sort -@ {threads} -o {output.bam}
        """

rule post_filter:
    """
    Filter MAPQ, proper pairs, fixmate, resort
    """
    input:
        bam = rules.align.output.bam
    output:
        filtered = "{out}/{sra}.filtered.sorted.bam"
    params:
        scratch = SCRATCH,
        mapq    = MAPQ
    log:
        "{out}/logs/filter_{sra}.log"
    threads: THREADS
    resources:
        mem_mb=8000
    conda:
        "envs/pipeline.yaml"
    shell:
        r"""
        mkdir -p {SCRATCH}/logs {OUT}/logs
        samtools view -F 1804 -f 2 -q {params.mapq} -u {input.bam} | \
          samtools sort -n -@ {threads} -o {params.scratch}/tmp.{wildcards.sra}.nmsrt.bam
        samtools fixmate -r -m \
          {params.scratch}/tmp.{wildcards.sra}.nmsrt.bam \
          {params.scratch}/tmp.{wildcards.sra}.fixmate.bam
        samtools view -F 1804 -f 2 -u {params.scratch}/tmp.{wildcards.sra}.fixmate.bam | \
          samtools sort -@ {threads} -o {output.filtered}
        rm {params.scratch}/tmp.{wildcards.sra}.nmsrt.bam \
           {params.scratch}/tmp.{wildcards.sra}.fixmate.bam
        """

rule markdup:
    """
    Mark duplicates, remove them, index, and cleanup
    """
    input:
        filtered = rules.post_filter.output.filtered
    output:
        dedup = "{out}/{sra}.filtered.sorted.dedup.bam",
        bai   = "{out}/{sra}.filtered.sorted.dedup.bam.bai"
    params:
        scratch = SCRATCH
    log:
        "{out}/logs/markdup_{sra}.log"
    threads: THREADS
    resources:
        mem_mb=8000
    conda:
        "envs/pipeline.yaml"
    shell:
        r"""
        samtools markdup -f {out}/{wildcards.sra}.markdup.qc \
          {input.filtered} \
          {params.scratch}/tmp.{wildcards.sra}.filtered.sorted.bam \
          &> {log}
        mv {params.scratch}/tmp.{wildcards.sra}.filtered.sorted.bam {input.filtered}
        samtools view -F 1804 -f 2 -b {input.filtered} > {output.dedup}
        samtools index {output.dedup}
        rm {input.filtered}
        """