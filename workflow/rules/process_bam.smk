
rule bam_to_tagalign:
	"""
	Convert filtered BAM (pe or se) to sorted, and indexed tagAlign file 
	"""
	input:
		bam = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.filtered.sorted.dedup.bam")
	output:
		tagalign = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.sorted.tagAlign.gz"),
		index = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.sorted.tagAlign.gz.tbi")
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
				sort -k1,1V -k2,2n -k3,3n --parallel {threads} | \
				bgzip -c > {output.tagalign}

		elif [[ "{params.run_type}" == "PAIRED" ]]; then
			samtools sort -n {input.bam} | \
				bedtools bamtobed -bedpe -mate1 -i stdin | \
				awk 'BEGIN{{OFS="\\t"}}{{printf "%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n", $1, $2, $3, $9, $4, $5, $6, $10}}' | \
				sort -k1,1V -k2,2n -k3,3n --parallel {threads} | \
				bgzip -c > {output.tagalign}
		else
			echo "Unsupported run_type: {params.run_type}" >&2
			exit 1
		fi

		# index
		tabix -p bed {output.tagalign}

		"""

rule tn5_shift_tagalign:
	"""
	Apply Tn5 shift to tagAlign file for ATAC-seq
	"""
	input:
		tagalign = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.sorted.tagAlign.gz"),
	output:
		tagalign_shift = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.tn5.sorted.tagAlign.gz"),
		index_shift = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.tn5.sorted.tagAlign.gz.tbi")
	params:
		run_type = lambda wildcards: SAMPLE_CONFIG.loc[wildcards.accession, "run_type"],
		scratch = SCRATCH_DIR
	log:
		os.path.join(RESULTS_DIR, "logs", "shift_tagalign_{sample_name}_{assay}_{accession}.log")
	threads: THREADS
	resources:
		mem_mb = determine_mem_mb
	conda:
		"../envs/seq_tools.yml"
	shell:
		"""
		set -euo pipefail

		echo "Working on {wildcards.accession}..." > {log}

		if [[ "{wildcards.assay}" != "ATAC" ]]; then
			echo "Only apply Tn5 shift to ATAC-seq tagAligns!"
			exit 1
		fi

		zcat {input.tagalign} | \
			awk 'BEGIN{{FS="\\t"; OFS="\\t"}}{{ if ($6 == "+") {{$2 = $2 + 4}} else if ($6 == "-") {{$3 = $3 - 5}} print $0 }}' | \
			sort -k1,1V -k2,2n -k3,3n --parallel {threads} | \
			bgzip -c > {output.tagalign_shift}

		# index
		tabix -p bed {output.tagalign_shift}
		"""

def get_tagalign_file(assay, sample_name, accession):
	if assay == "ATAC":
		return os.path.join(RESULTS_DIR, sample_name, assay, f"{accession}.tn5.sorted.tagAlign.gz")
	else:
		return os.path.join(RESULTS_DIR, sample_name, assay, f"{accession}.sorted.tagAlign.gz")


def get_macs2_params(assay, param):
	if assay in [config["MACS2"]["accessibility"]["assays"]]:
		extsize = config["MACS2"]["accessibility"]["extsize"]
		shift = config["MACS2"]["accessibility"]["shift"]
	else:
		extsize = "fraglen"
		shift = 0

	if param == "extsize":
		return extsize
	elif param == "shift":
		return shift
	else:
		raise Exception(f"{param} is not a valid parameter.")


rule macs2_call_peaks:
	"""
	Call peaks with MACS2 with settings for ChIP-seq or accessibility depending on assay; only restrict based on p-value
	"""
	input:
		tagalign = lambda wildcards: get_tagalign_file(wildcards.assay, wildcards.sample_name, wildcards.accession)
	output:
		narrowpeak = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.macs2_peaks.narrowPeak"),
		macs2_cmd = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.macs2_command.txt")
	params:
		results = RESULTS_DIR,
		scratch = SCRATCH_DIR,
		genome = config["MACS2"]["genome"],
		p_threshold = config["MACS2"]["p_value"],
		extsize = lambda wildcards: get_macs2_params(wildcards.assay, "extsize"),
		shift = lambda wildcards: get_macs2_params(wildcards.assay, "shift")
	log:
		os.path.join(RESULTS_DIR, "logs", "macs2_{sample_name}_{assay}_{accession}.log")
	threads: THREADS
	resources:
		mem_mb = determine_mem_mb
	conda:
		"../envs/seq_tools.yml"
	shell:
		"""
		out_dir="{params.results}/{wildcards.sample_name}/{wildcards.assay}"
		prefix="{wildcards.accession}.macs2"
		
		# estimate or set extension size
		if [[ "{params.extsize}" == "fraglen" ]]; then
			echo "Estimating fragment length..." >> {log} 2>&1
			predictd_output="$out_dir/${prefix}_predictd.txt"

			macs2 predictd -i {input.tagalign} \
				-f BED \
				-g {params.genome} \
				-n $prefix \
				--outdir $out_dir \
				> $predictd_output 2>&1

			extsize=$(grep -oP "predicted fragment length is \K[0-9]+" $predictd_output)
			echo "Estimated d: $extsize" >> {log} 2>&1

		else
			extsize="{params.extsize}"
			echo "Using preset fragment length: $extsize" >> {log} 2>&1
		fi

		# call peaks
		macs2_cmd="macs2 callpeak -t {input.tagalign} \
			-f BED \
			-g {params.genome} \
			-n $prefix \
			-p {params.p_threshold} \
			--nomodel \
			--extsize $extsize \
			--shift {params.shift} \
			--keep-dup all \
			--call-summits \
			--outdir $out_dir"

		echo "$macs2_cmd" > {output.macs2_cmd}
		echo "Running: $macs2_cmd" >> {log} 2>&1
		eval "$macs2_cmd" >> {log} 2>&1

		"""
