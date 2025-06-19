## Update paths in the config obj to absolute path
def make_paths_absolute(obj, base_path):
	"""
	Use absolute paths to be compatible with github submodules
	Recursively go through the dictionary and convert relative paths to absolute paths.
	"""
	if isinstance(obj, dict):
		for key, value in obj.items():
			obj[key] = make_paths_absolute(value, base_path)
	elif isinstance(obj, str):
		# We assume all strings are paths. If converting the string
		# to an absolute path results in a valid file, then the str was a path
		new_file = os.path.join(base_path, obj)
		if os.path.exists(new_file):
			return new_file
	return obj

def determine_mem_mb(wildcards, input, attempt, min_gb=8):
	# Memory resource calculator for snakemake rules
	input_size_mb = input.size_mb
	if ".gz" in str(input):
		input_size_mb *= 8  # assume gz compressesed the file <= 8x
	attempt_multiplier = 2 ** (attempt - 1)  # Double memory for each retry
	mem_to_use_mb = attempt_multiplier *  max(4 * input_size_mb, min_gb * 1000)
	return min(mem_to_use_mb, MAX_MEM_MB)

def validate_sample_config(sample_config):
	required_cols = {"accession", "sample_name", "assay"}
	missing = required_cols - set(SAMPLE_CONFIG.columns)
	if missing:
		raise ValueError(f"Missing required columns in sample config: {missing}")

def get_fastq_files(accession, sample_config, config, scratch_dir):
    if accession not in sample_config.index:
        raise KeyError(f"Accession '{accession}' not found in sample metadata")

    sample = sample_config.loc[accession]
    sample_name = sample["sample_name"]
    assay = sample["assay"]
    run_type = sample["run_type"]

    if config["sra_download"]:
        if run_type == "PAIRED":
            return [
                os.path.join(scratch_dir, sample_name, assay, f"{accession}_1.fastq.gz"),
                os.path.join(scratch_dir, sample_name, assay, f"{accession}_2.fastq.gz"),
            ]
        elif run_type == "SINGLE":
            return [
                os.path.join(scratch_dir, sample_name, assay, f"{accession}.fastq.gz")
            ]
        else:
            raise ValueError(f"Unknown run_type: {run_type}")
    else:
        # Use pre-computed fastq paths
        fastq_str = sample.get("fastq_files")
        if not fastq_str:
            raise ValueError(f"No fastq_files entry for {accession}")

        paths = [s.strip() for s in fastq_str.split(",") if s.strip()]

        if run_type == "PAIRED":
            if len(paths) != 2:
                raise ValueError(f"Expected 2 fastq files for PAIRED run_type for {accession}, got: {paths}")
            return paths
        elif run_type == "SINGLE":
            if len(paths) != 1:
                raise ValueError(f"Expected 1 fastq file for SINGLE run_type for {accession}, got: {paths}")
            return paths
        else:
            raise ValueError(f"Unknown run_type: {run_type}")
