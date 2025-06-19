rule save_reference_configs:
    """
    Save resolved config file
    Save input sample config with per-sample output files
    """
    params:
        main_config = config,
        sample_config = SAMPLE_CONFIG,
        out_dir = os.path.join(RESULTS_DIR, "config"), 
        results_dir = RESULTS_DIR,
        scratch_dir = SCRATCH_DIR
    output:
        main_config_out = os.path.join(RESULTS_DIR, "config", "reference_config.yml"),
        sample_config_out = os.path.join(RESULTS_DIR, "config", "sample_metadata.tsv")
    run:
        import os, yaml, pandas as pd

        # make sure the output dir exists
        os.makedirs(params.out_dir, exist_ok=True)

        # 1) save main config
        with open(output.main_config_out, 'w') as fh:
            yaml.safe_dump(params.main_config, fh)

        # 2) augment biosample table
        df = params.sample_config

        # helper to construct file paths for each sample
            # fastq if necessary
            # bam 
            # peaks if necessary
        def mkpaths(row):
            accession = row["accession"] 
            sample_name = row["sample_name"]
            assay = row["assay"]
            run_type = row["run_type"]
            out = {}


            # bam files
            out["bam"] = os.path.join(params.scratch_dir, sample_name, assay, f"{accession}.filtered.sorted.dedup.bam.bai")

            # fastq if indicated
            if params.main_config["sra_download"]:
                if run_type == "PAIRED":
                    f12 = [os.path.join(params.scratch_dir, sample_name, assay, f"{accession}_{n}.fastq") for n in [1, 2]]
                    out["fastq_download"] = ",".join(f12)
                elif run_type == "SINGLE":
                    out["fastq_download"] = os.path.join(params.scratch_dir, sample_name, assay, f"{accession}.fastq")

            # peaks if indicated
            if params.main_config["call_peaks"]:
                out["peaks"] = os.path.join(params.scratch_dir, sample_name, assay, f"{accession}.sorted.narrowPeak")

            return pd.Series(out)

        df = pd.concat([df, df.apply(mkpaths, axis=1)], axis=1)
        df.to_csv(output.sample_config_out, sep='\t', index=False)