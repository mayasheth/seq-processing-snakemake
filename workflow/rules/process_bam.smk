## Rules that are common to single and paired-ended runs

rule prepare_bowtie2_index:
  """
  Copy Bowtie2 index files to scratch for faster access
  """
  input:
    expand("{index}.1.bt2", index=config["bowtie2_index"])  # or .1.bt2l for large indexes
  output:
    index_base = os.path.join(SCRATCH_DIR, "bowtie2_index", os.path.basename(config["bowtie2_index"]) + ".1.bt2")
  params:
    index = config["bowtie2_index"],
    scratch = os.path.join(SCRATCH_DIR, "bowtie2_index")
  shell:
    """
    mkdir -p {params.scratch}
    cp {params.index}.* {params.scratch}/
    """


rule samtools_stats:
  input:
    bam = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.{runtype}.filtered.sorted.dedup.bam")
  output:
    stats = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.{runtype}.samtools_stats.txt")
  resources:
    mem_mb=3000,
    runtime=30
  conda:
    "../envs/seq_tools.yml"
  shell:
    """
    samtools stats {input.bam} > {output.stats}
    """

rule bam_to_bigwig: 
  input:
    bam = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.{runtype}.filtered.sorted.dedup.bam"),
    chrom_sizes = config["chrom_sizes"]
  output:
    bam_filt = temp(os.path.join(SCRATCH_DIR, "{sample_name}", "{assay}", "{accession}.{runtype}.chrom_sizes_filt.bam")),
    bam_filt_index = temp(os.path.join(SCRATCH_DIR, "{sample_name}", "{assay}", "{accession}.{runtype}.chrom_sizes_filt.bam.bai")),
    bg = temp(os.path.join(SCRATCH_DIR, "{sample_name}", "{assay}", "{accession}.{runtype}.bg")),
    bw = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.{runtype}.bw")
  params:
    se_fragment_size = config["se_fragment_size"]
  resources:
    mem_mb=determine_mem_mb,
    runtime=360
  log:
    os.path.join(RESULTS_DIR, "logs", "bam_to_bw.{sample_name}_{assay}_{accession}_{runtype}.log")
  conda:
    "../envs/ucsc_tools.yml"
  shell:
    """
    set -euo pipefail

    echo "Starting BAM to BigWig conversion for {wildcards.accession}..." > {log}

    # Filter BAM file to chrom sizes
    echo "Filtering BAM to chromosomes in chrom_sizes..." >> {log}
    samtools view -b {input.bam} -o {output.bam_filt} 2>> {log}
    
    # Filter to only chromosomes in chrom_sizes file
    samtools view -h {output.bam_filt} | \
    awk 'BEGIN{{while((getline < "{input.chrom_sizes}") > 0) chroms[$1]=1}} 
         /^@/ || chroms[$3]' | \
    samtools view -b -o {output.bam_filt}.tmp - 2>> {log}
    mv {output.bam_filt}.tmp {output.bam_filt}
    
    samtools index {output.bam_filt} 2>> {log}

    # Generate bedGraph
    echo "Generating bedGraph..." >> {log}
    if [[ {wildcards.runtype} == "se" ]]; then
        echo "Using fragment size of {params.se_fragment_size} bp" >> {log}
        bedtools genomecov -ibam {output.bam_filt} -bg -fs {params.se_fragment_size} > {output.bg} 2>> {log}
    elif [[ {wildcards.runtype} == "pe" ]]; then
        echo "Using paired-end mode" >> {log}
        bedtools genomecov -ibam {output.bam_filt} -bg -pc > {output.bg} 2>> {log}
    fi

    # Sort bedGraph (required for bedGraphToBigWig)
    echo "Sorting bedGraph..." >> {log}
    sort -k1,1 -k2,2n {output.bg} > {output.bg}.sorted 2>> {log}
    mv {output.bg}.sorted {output.bg}

    # Convert to bigWig
    echo "Converting bedGraph to bigWig..." >> {log}
    bedGraphToBigWig {output.bg} {input.chrom_sizes} {output.bw} 2>> {log}
    
    echo "Completed BAM to BigWig conversion" >> {log}
    """

rule bam_to_tagalign:
  """
  Convert filtered BAM (pe or se) to sorted, and indexed tagAlign file 
  """
  input:
    bam = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.{runtype}.filtered.sorted.dedup.bam")
  output:
    tagalign = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.{runtype}.sorted.tagAlign.gz"),
    index = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.{runtype}.sorted.tagAlign.gz.tbi")
  params:
    scratch = SCRATCH_DIR
  log:
    os.path.join(RESULTS_DIR, "logs", "tagalign.{sample_name}_{assay}_{accession}_{runtype}.log")
  threads: THREADS
  resources:
    mem_mb = determine_mem_mb
  conda:
    "../envs/seq_tools.yml"
  shell:
    """
    set -euo pipefail

    echo "Working on {wildcards.accession}..." > {log}

    if [[ "{wildcards.runtype}" == "se" ]]; then
      samtools sort -n {input.bam} | \
        bedtools bamtobed -i stdin | \
        awk 'BEGIN{{OFS="\\t"}}{{$4="N";$5="1000";print $0}}' | \
        sort -k1,1V -k2,2n -k3,3n --parallel {threads} | \
        bgzip -c > {output.tagalign}

    elif [[ "{wildcards.runtype}" == "pe" ]]; then
      samtools sort -n {input.bam} | \
        bedtools bamtobed -bedpe -mate1 -i stdin | \
        awk 'BEGIN{{OFS="\\t"}}{{printf "%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n", $1, $2, $3, $9, $4, $5, $6, $10}}' | \
        sort -k1,1V -k2,2n -k3,3n --parallel {threads} | \
        bgzip -c > {output.tagalign}
    else
      echo "Unsupported runtype: {wildcards.runtype}" >&2
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
    tagalign = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.{runtype}.sorted.tagAlign.gz"),
  output:
    tagalign_shift = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.{runtype}.tn5.sorted.tagAlign.gz"),
    index_shift = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.{runtype}.tn5.sorted.tagAlign.gz.tbi")
  params:
    run_type = lambda wildcards: SAMPLE_CONFIG.loc[wildcards.accession, "run_type"],
    scratch = SCRATCH_DIR
  log:
    os.path.join(RESULTS_DIR, "logs", "shift_tagalign.{sample_name}_{assay}_{accession}_{runtype}.log")
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

def get_tagalign_file(assay, sample_name, accession, runtype):
  if assay == "ATAC":
    return os.path.join(RESULTS_DIR, sample_name, assay, f"{accession}.{runtype}.tn5.sorted.tagAlign.gz")
  else:
    return os.path.join(RESULTS_DIR, sample_name, assay, f"{accession}.{runtype}.sorted.tagAlign.gz")


def get_macs2_params(assay, param):
  if assay in [config["MACS2"]["accessibility_params"]["assays"]]:
    extsize = config["MACS2"]["accessibility_params"]["extsize"]
    shift = config["MACS2"]["accessibility_params"]["shift"]
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
    tagalign = lambda wildcards: get_tagalign_file(wildcards.assay, wildcards.sample_name, wildcards.accession, wildcards.runtype)
  output:
    narrowpeak = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.{runtype}.macs2_peaks.narrowPeak"),
    macs2_cmd = os.path.join(RESULTS_DIR, "{sample_name}", "{assay}", "{accession}.{runtype}.macs2_command.txt")
  params:
    results = RESULTS_DIR,
    scratch = SCRATCH_DIR,
    genome = config["MACS2"]["genome"],
    p_threshold = config["MACS2"]["p_value"],
    extsize = lambda wildcards: get_macs2_params(wildcards.assay, "extsize"),
    shift = lambda wildcards: get_macs2_params(wildcards.assay, "shift")
  log:
    os.path.join(RESULTS_DIR, "logs", "macs2.{sample_name}_{assay}_{accession}_{runtype}.log")
  threads: THREADS
  resources:
    mem_mb = determine_mem_mb
  conda:
    "../envs/seq_tools.yml"
  shell:
    """
    out_dir="{params.results}/{wildcards.sample_name}/{wildcards.assay}"
    prefix="{wildcards.accession}.{wildcards.runtype}.macs2"
    
    # estimate or set extension size
    if [[ "{params.extsize}" == "fraglen" ]]; then
      echo "Estimating fragment length..." >> {log} 2>&1
      predictd_output="$out_dir/${{prefix}}_predictd.txt"

      macs2 predictd -i {input.tagalign} \
        -f BED \
        -g {params.genome} \
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
