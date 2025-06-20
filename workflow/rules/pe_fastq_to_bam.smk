rule align_pe_fastq:
	"""
	Bowtie2 → sorted BAM
	"""
	input:
		fq1 = lambda wildcards: get_fastq_files(wildcards.accession, SAMPLE_CONFIG, config, SCRATCH_DIR)[0],
		fq2 = lambda wildcards: get_fastq_files(wildcards.accession, SAMPLE_CONFIG, config, SCRATCH_DIR)[1],
	params:
		results = RESULTS_DIR,
		index = config["bowtie2_index"]
	output:
		bam = temp(os.path.join(SCRATCH_DIR, "{sample_name}", "{assay}", "{accession}.bam"))
	log:
		 os.path.join(RESULTS_DIR, "logs", "align_{sample_name}_{assay}_{accession}.log")
	threads: THREADS
	resources:
		mem_mb = determine_mem_mb,
		runtime_min = 12 * 60
	conda:
		"../envs/seq_tools.yml"
	shell:
		"""
		mkdir -p {params.results}/logs

		bowtie2 -X 2000 --mm \
		  -x {params.index} \
		  -1 {input.fq1} -2 {input.fq2} \
		  -p {threads} 2> {log} | \
		samtools view -h -@ {threads} - | \
		samtools sort -@ {threads} -o {output.bam}
		"""

rule post_filter_pe_bam:
	"""
	Filter MAPQ, proper pairs, fixmate, resort
	"""
	input:
		bam = rules.align_pe_fastq.output.bam
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
		samtools view -F 1804 -f 2 -q {params.mapq} -u {input.bam} | \
		  samtools sort -n -@ {threads} -o {params.scratch}/tmp.{wildcards.accession}.nmsrt.bam
		samtools fixmate -r -m \
		  {params.scratch}/tmp.{wildcards.accession}.nmsrt.bam \
		  {params.scratch}/tmp.{wildcards.accession}.fixmate.bam
		samtools view -F 1804 -f 2 -u {params.scratch}/tmp.{wildcards.accession}.fixmate.bam | \
		  samtools sort -@ {threads} -o {output.filtered}
		rm {params.scratch}/tmp.{wildcards.accession}.nmsrt.bam \
		   {params.scratch}/tmp.{wildcards.accession}.fixmate.bam
		"""

rule markdup_pe:
	"""
	Mark duplicates, remove them, index, and cleanup
	"""
	input:
		filtered = rules.post_filter_pe_bam.output.filtered
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
		samtools markdup -f {params.scratch}/{wildcards.accession}.markdup.qc \
		  {input.filtered} \
		  {params.scratch}/tmp.{wildcards.accession}.filtered.sorted.bam \
		  &> {log}
		mv {params.scratch}/tmp.{wildcards.accession}.filtered.sorted.bam {input.filtered}
		samtools view -F 1804 -f 2 -b {input.filtered} > {output.dedup}
		samtools index {output.dedup}
		"""

rule bam_to_tagalign:
	"""
	Convert filtered BAM to shifted, sorted, and indexed tagAlign file.
	"""
	input:
		bam = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.filtered.sorted.dedup.bam")
	output:
		tagalign_init = temp(os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "initial_{accession}.tn5.sorted.tagAlign.gz")),
		tagalign = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.tn5.sorted.tagAlign.gz"),
		index = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.tn5.sorted.tagAlign.gz.tbi")
	params:
		run_type = lambda wildcards: SAMPLE_CONFIG.loc[wildcards.accession, "run_type"],
		scratch = SCRATCH_DIR
	log:
		os.path.join(RESULTS_DIR, "logs", "tagalign_{sample_name}_{assay}_{accession}.log")
	threads: THREADS
	resources:
		mem_mb = determine_mem_mb
	conda:
		"../envs/seq_tools.yml"
	shell:
		"""
		set -euo pipefail

		echo "Working on {wildcards.accession}..." > {log}

		if [[ "{params.run_type}" == "SINGLE" ]]; then
			samtools sort -n {input.bam} | \
				bedtools bamtobed -i stdin | \
				awk 'BEGIN{{OFS="\\t"}}{{$4="N";$5="1000";print $0}}' | \
				gzip -nc > {output.tagalign_init}

		elif [[ "{params.run_type}" == "PAIRED" ]]; then
			samtools sort -n {input.bam} | \
				bedtools bamtobed -bedpe -mate1 -i stdin | \
				awk 'BEGIN{{OFS="\\t"}}{{printf "%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n", $1, $2, $3, $9, $4, $5, $6, $10}}' | \
				gzip -nc > {output.tagalign_init}
		else
			echo "Unsupported run_type: {params.run_type}" >&2
			exit 1
		fi

		zcat {output.tagalign_init} | \
			awk 'BEGIN{{FS="\\t"; OFS="\\t"}}{{ if ($6 == "+") {{$2 = $2 + 4}} else if ($6 == "-") {{$3 = $3 - 5}} print $0 }}' | \
			sort -k1,1V -k2,2n -k3,3n --parallel {threads} | \
			bgzip -c > {output.tagalign}

		tabix -p bed {output.tagalign}
		"""
