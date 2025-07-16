def get_period_run_info_mgc(filename):
    """
    Extract period and run information from the filename.
    Looks for patterns like 'p##' for period and 'r###' for run.
    """
    import re
    parts = filename.split('/')
    
    period = None
    run = None
    
    # Look for period pattern (p followed by digits)
    for part in parts:
        if re.match(r'^p\d+$', part):
            period = part
            break
    
    # Look for run pattern (r followed by digits) 
    for part in parts:
        if re.match(r'^r\d+$', part):
            run = part
            break
    
    # If we didn't find the patterns, try to extract from filename itself
    if period is None or run is None:
        basename = parts[-1]  # Get just the filename
        # Look for period pattern in filename
        period_match = re.search(r'p\d+', basename)
        if period_match:
            period = period_match.group()
        
        # Look for run pattern in filename  
        run_match = re.search(r'r\d+', basename)
        if run_match:
            run = run_match.group()
    
    return period, run


rule make_sum_mgc_per_run:
    """
    Collect mgc.lh5 files for a specific period/run and concatenate them.
    """
    input:
        barrier="/tmp/mgc_all.done"
    output:
        "gen/mu_hpge_coinc/{period}/sum_mgc_{run}.lh5"
    run:
        import os
        import lgdo.lh5 as lh5
        import lgdo.types as types
        import awkward as ak
        from tqdm import tqdm
        from ge77m_search_workflow import utils as ut

        period = wildcards.period
        run = wildcards.run
        output_file = str(output[0])
        
        if os.path.exists(output_file):
            print(f"Output file {output_file} already exists. Skipping processing.")
            return

        # Find all mgc files for this specific period/run
        mgc_pattern = f"gen/mu_hpge_coinc/{period}/{run}/*.lh5"
        mgc_files = sorted(glob.glob(mgc_pattern))
        
        if not mgc_files:
            print(f"Warning: No mgc files found for pattern {mgc_pattern}")
            # Create empty file to satisfy Snakemake
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            with open(output_file, 'w') as f:
                pass
            return

        print(f"Processing {len(mgc_files)} files for {period}/{run}")
        
        concatenated_data = []
        for file in tqdm(mgc_files, desc=f"Processing {period}/{run}"):
            try:
                data = lh5.read_as("mgc", file, "ak")
                concatenated_data.append(data)
            except Exception as e:
                print(f"Warning: Could not read {file}: {e}")
                continue
        
        if concatenated_data:
            final_data = types.Table(ak.concatenate(concatenated_data))
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            lh5.write(final_data, name="mgc", lh5_file=output_file)
        else:
            print(f"Warning: No valid data found for {period}/{run}")
            # Create empty file to satisfy Snakemake
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            with open(output_file, 'w') as f:
                pass


rule make_sum_mgc_global:
    """
    Collect all per-run sum files and create global sum_mgc.lh5.
    """
    input:
        sum_files=lambda wildcards: [
            f"gen/mu_hpge_coinc/{period}/sum_mgc_{run}.lh5"
            for period, run in {(d["lvl1"], d["lvl2"]) for d in initial}
        ]
    output:
        "gen/mu_hpge_coinc/sum_mgc.lh5"
    run:
        import os
        import lgdo.lh5 as lh5
        import lgdo.types as types
        import awkward as ak
        from tqdm import tqdm

        output_file = str(output[0])
        
        if os.path.exists(output_file):
            print(f"Output file {output_file} already exists. Skipping processing.")
            return

        valid_sum_files = [f for f in input.sum_files if os.path.exists(f) and os.path.getsize(f) > 0]
        
        if not valid_sum_files:
            raise ValueError("No valid sum files found for global concatenation. This should not happen if dependencies are correct.")

        print(f"Processing {len(valid_sum_files)} sum files for global concatenation")
        
        concatenated_data = []
        for file in tqdm(valid_sum_files, desc="Processing sum files"):
            try:
                data = lh5.read_as("mgc", file, "ak")
                concatenated_data.append(data)
            except Exception as e:
                print(f"Warning: Could not read {file}: {e}")
                continue
        
        if concatenated_data:
            final_data = types.Table(ak.concatenate(concatenated_data))
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            lh5.write(final_data, name="mgc", lh5_file=output_file)
        else:
            raise ValueError("No valid data found for global concatenation")
