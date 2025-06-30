
sum_mdc_files = sorted(glob.glob("gen/mu_delayed_coinc/*/sum_mdc_*.lh5"))

def get_period_run_info(file_path):
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
            period, run = get_period_run_info(file)
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

            period, run = get_period_run_info(sum_file)

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

            output_data = sum_mdc[mask]

            lh5.write(
                types.Table(output_data),
                name="mdc",
                lh5_file=skm_file
            )

        output_data = []
        for skm_file in runwise_skm_files:
            tmp = lh5.read_as("mdc", skm_file, "ak")
            output_data.append(tmp)
        final_data = types.Table(ak.concatenate(output_data))
        lh5.write(final_data, name="mdc", lh5_file=str(output))
        