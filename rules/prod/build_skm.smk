
sum_mdc_files = sorted(glob.glob("gen/mu_delayed_coinc/*/sum_mdc_*.lh5"))

def get_period_run_info_skm(file_path):
    """
    Extract period and run information from the file path.
    Assumes the format is gen/mu_delayed_coinc/{period}/sum_mdc_{run}.lh5
    """
    parts = file_path.split('/')
    period = parts[-2]
    run = parts[-1].replace("sum_mdc_", "").replace(".lh5", "")
    return period, run



rule make_skm:
    """
    Collect all mgc.lh5 files from the mu_hpge_coinc directory and concatenate them into sum_mgc.lh5.
    """

    input:
        mdc_sum_file= config["out_root"] + "/mu_delayed_coinc/sum_mdc.lh5",
    output:
        os.path.join("gen/skm/skm_" + "{psd_ms_type}" + ".lh5")
    params:
        psd_ms_type="{psd_ms_type}",
    run:
        import lgdo.lh5 as lh5
        import lgdo.types as types
        import awkward as ak
        from tqdm import tqdm
        from ge77m_search_workflow import utils as ut
        import legendmeta as lm
        import numpy as np

        lmeta = lm.LegendMetadata(config["metadata"])

        if "-" in params.psd_ms_type:
            psd_ms_type = params.psd_ms_type.split("-")
        else:
            psd_ms_type = params.psd_ms_type

        runwise_skm_files = []
        for file in sum_mdc_files:
            period, run = get_period_run_info_skm(file)
            skm_file = f"gen/skm/{period}/skm_{run}_" + "{}.lh5"
            runwise_skm_files.append(skm_file)

        if not os.path.exists("gen/skm"):
            os.makedirs("gen/skm", exist_ok=True)

        runwise_skm_files = ([f.format(params.psd_ms_type) for f in runwise_skm_files])

        physics_dataset = lmeta.datasets.runlists[config["physics_dataset"]]["phy"]

        for sum_file, skm_file in tqdm(zip(sum_mdc_files, runwise_skm_files), total=len(sum_mdc_files)):

            if os.path.exists(skm_file):
                continue

            # Read the sum_mdc file
            sum_mdc = lh5.read_as("mdc/",sum_file,"ak")
            
            mask_univ = (
                (sum_mdc["prompt"]["geds"]["energy"] > config["selection"]["main"]["prompt_e_range"][0])
                & (sum_mdc["prompt"]["geds"]["energy"] < config["selection"]["main"]["prompt_e_range"][1])
                & (sum_mdc["prompt"]["mu"]["coinc_flags"]["muon_offline"])
                & (sum_mdc["delayed"]["geds"]["quality"]["quality_is_bb_like"])
            )

            period, run = get_period_run_info_skm(sum_file)

            runs_in_physics_period = ut.expand_range_token(physics_dataset[period])

            if run not in runs_in_physics_period:
                continue

            if not os.path.exists("gen/skm/"+period):
                os.makedirs("gen/skm/"+period, exist_ok=True)

            runinfo = lmeta.datasets.runinfo[period][run].phy
            chmap = ut.generate_channel_map(runinfo.start_key, metadata=lmeta)

            data_streams_hpge = ut.select_datastreams(chmap, "HPGE",psd_ms_type=psd_ms_type)

            mask_psd_ms_type = np.sum(np.array([
                rawid == sum_mdc["delayed"]["geds"]["id"]["hit_table"] for rawid in list(data_streams_hpge)
            ]).T,axis=-1) > 0

            mask = mask_univ & mask_psd_ms_type

            # Apply mask to get filtered data
            filtered_data = sum_mdc[mask]
            
            # Check if we have any data after filtering
            if len(filtered_data) == 0:
                # No data after filtering, create empty structure
                output_data = filtered_data  # Use the original empty filtered data structure
            else:
                # Group by unique delayed events (hit_table + timestamp)
                # First, sort by hit_table
                sorted_indices = ak.argsort(filtered_data["delayed"]["geds"]["id"]["hit_table"])
                data_sorted_by_hit_table = filtered_data[sorted_indices]
            
                # Group by hit_table
                hit_table_vals = data_sorted_by_hit_table["delayed"]["geds"]["id"]["hit_table"]
                hit_table_breaks = np.where(np.diff(hit_table_vals, prepend=hit_table_vals[0]-1, append=hit_table_vals[-1]+1) != 0)[0]
                hit_table_lengths = np.diff(hit_table_breaks)
                data_grouped_by_hit_table = ak.unflatten(data_sorted_by_hit_table, hit_table_lengths)
                
                # Within each hit_table group, sort by timestamp and group by unique timestamps
                data_grouped_by_hit_table_sorted_by_ts = data_grouped_by_hit_table[
                    ak.argsort(data_grouped_by_hit_table["delayed"]["geds"]["id"]["timestamp"], axis=1)
                ]
                
                # Group by unique timestamps within each hit_table group
                restructured_groups = []
                for hit_table_group in data_grouped_by_hit_table_sorted_by_ts:
                    if len(hit_table_group) == 0:
                        continue
                        
                    timestamps = hit_table_group["delayed"]["geds"]["id"]["timestamp"]
                    ts_breaks = np.where(np.diff(timestamps, prepend=timestamps[0]-1, append=timestamps[-1]+1) != 0)[0]
                    ts_lengths = np.diff(ts_breaks)
                    
                    if len(ts_lengths) > 0:
                        ts_groups = ak.unflatten(hit_table_group, ts_lengths)
                        restructured_groups.append(ts_groups)
                
                # Restructure the data into the desired format: one record per unique delayed event
                if restructured_groups:
                    # Concatenate all groups
                    all_groups = ak.concatenate(restructured_groups)
                    
                    # Create final records list
                    final_records = []
                    
                    # Process each group (each group contains events with the same delayed hit_table + timestamp)
                    for group in all_groups:
                        if len(group) == 0:
                            continue
                            
                        # Get the unique delayed entry (should be the same for all entries in this group)
                        delayed_entry = group["delayed"][0]  # Take the first one since they're all the same
                        
                        # Get all the dT values for this delayed event
                        dT_list = ak.to_list(group["dT"])
                        
                        # Get all the prompt entries for this delayed event
                        prompt_list = ak.to_list(group["prompt"])
                        
                        # Create the record
                        record = {
                            "delayed": delayed_entry,
                            "dT": dT_list,
                            "prompts": prompt_list
                        }
                        
                        final_records.append(record)
                    
                    # Convert back to awkward array format for LH5 writing
                    if final_records:
                        # Create arrays from the records
                        delayed_array = ak.Array([record["delayed"] for record in final_records])
                        dT_array = ak.Array([record["dT"] for record in final_records])
                        
                        # Extract prompt ID information as separate arrays
                        prompt_hit_table_array = ak.Array([[p["geds"]["id"]["hit_table"] for p in record["prompts"]] for record in final_records])
                        prompt_timestamp_array = ak.Array([[p["coinc"]["id"]["timestamp"] for p in record["prompts"]] for record in final_records])
                        
                        output_data = ak.Array({
                            "delayed": delayed_array,
                            "dT": dT_array,
                            "prompt_hit_table": prompt_hit_table_array,
                            "prompt_timestamp": prompt_timestamp_array
                        })
                    else:
                        output_data = filtered_data[:0]  # Empty array with same structure
                else:
                    output_data = filtered_data[:0]  # Empty array with same structure
            print(output_data)
            if len(output_data) == 0:
                print(f"No valid output data available for {skm_file}.")
                

            lh5.write(
                types.Table(output_data),
                name="skm",
                lh5_file=skm_file
            )

        output_data = []
        for skm_file in runwise_skm_files:
            if not os.path.exists(skm_file):
                print(f"File {skm_file} does not exist, skipping.")
                continue
            tmp = lh5.read_as("skm", skm_file, "ak")
            if len(tmp) == 0:
                print(f"No data in {skm_file}, skipping.")
                continue
            output_data.append(tmp)
        print(output_data)
        final_data = types.Table(ak.concatenate(output_data))
        lh5.write(final_data, name="skm", lh5_file=str(output))
        