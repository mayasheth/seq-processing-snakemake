rule align_se_fastq:
	"""
	Bowtie2 → sorted BAM
	"""
	input:
		fq1 = lambda wildcards: get_fastq_files(wildcards.accession, SAMPLE_CONFIG, config, SCRATCH_DIR)[0]
		scratch_index = os.path.join(SCRATCH_DIR, "bowtie2_index", os.path.basename(config["bowtie2_index"]) + ".1.bt2")
	params:
		results = RESULTS_DIR,
		index = os.path.join(SCRATCH_DIR, "bowtie2_index", os.path.basename(config["bowtie2_index"]))
	output:
		bam = temp(os.path.join(SCRATCH_DIR, "{sample_name}", "{assay}", "{accession}.bam"))
	log:
		os.path.join(RESULTS_DIR, "logs", "align_{sample_name}_{assay}_{accession}.log")
	threads: THREADS
	resources:
		mem_mb = determine_mem_mb
	conda:
		"../envs/seq_tools.yml"
	shell:
		"""
		set -euo pipefail

		echo "Start: $(date)" > {log}
		echo "Threads: {threads}" >> {log}

		# Align
		bowtie2 --mm \
			-x {params.index} \
			-U {input.fq} \
			-p {threads} 2>> {log} \
		samtools sort -@ {threads} -o {output.bam} - >> {log} 2>&1

		echo "Finished: $(date)" >> {log}
		"""

rule post_filter_se_bam:
	"""
	Filter MAPQ 
	"""
	input:
		bam = rules.align_se_fastq.output.bam
	output:
		filtered = temp(os.path.join(SCRATCH_DIR, "{sample_name}", "{assay}", "{accession}.filtered.sorted.bam"))
	params:
		scratch = SCRATCH_DIR,
		mapq    = MAPQ
	log:
		os.path.join(RESULTS_DIR, "logs", "filter_{sample_name}_{assay}_{accession}.log")
	threads: THREADS
	resources:
		mem_mb = determine_mem_mb
	conda:
		"../envs/seq_tools.yml"
	shell:
		"""
		samtools view -F 1804 -q {params.mapq} -b ${input.bam} -o ${output.filtered}
		"""

rule markdup_se_bam:
	"""
	Mark duplicates, remove them, index, and cleanup
	"""
	input:
		filtered = rules.post_filter_se_bam.output.filtered
	output:
		dedup = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.filtered.sorted.dedup.bam"),
		bai = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.filtered.sorted.dedup.bam.bai")
	params:
		scratch = SCRATCH_DIR
	log:
		os.path.join(RESULTS_DIR, "logs", "markdup_{sample_name}_{assay}_{accession}.log")
	threads: THREADS
	resources:
		mem_mb = determine_mem_mb
	conda:
		"../envs/seq_tools.yml"
	shell:
		"""
		samtools markdup -@ {threads} -f {params.scratch}/{wildcards.accession}.markdup.qc \
			{input.filtered} \
			{params.scratch}/tmp.{wildcards.accession}.filtered.sorted.bam \
			&> {log}
		
		mv {params.scratch}/tmp.{wildcards.accession}.filtered.sorted.bam {input.filtered}
		
		samtools view -F 1804 -b {input.filtered} > {output.dedup}
		samtools index {output.dedup}
		"""