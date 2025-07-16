# rules/delayed_coinc.smk

configfile: "config.yaml"

# Format input_root with the default_ref_version from config
input_root = config["input_root"].format(default_ref_version=config["default_ref_version"])


rule delayed_coinc:
    """
    Take the mu_hpge_coinc output and produce
    tier_mdc.lh5 in gen/mu_delayed_coinc/...
    """

    wildcard_constraints:
        base=".*tier_mdc.*"

    input:
        mgc_sum_file=(lambda wc: config["out_root"] + "/mu_hpge_coinc/{lvl1}/sum_mgc_{lvl2}.lh5".format(
            lvl1=wc["lvl1"],
            lvl2=wc["lvl2"]
        )),
        pht_files=(lambda wc: input_root + "/{lvl1}/{lvl2}/{base}.lh5".format(
            lvl1=wc["lvl1"],
            lvl2=wc["lvl2"],
            base=wc["base"].replace("tier_mdc", "tier_pht")
        )),
        barrier="/tmp/mgc_all.done",
        #waveform_block="/tmp/mgc_plots_all.done"
    output:
        os.path.join(
            config["out_root"],
            "mu_delayed_coinc",
            "{lvl1}",
            "{lvl2}",
            "{base}" + ".lh5"
        )
    run:
        import os
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        from ge77m_search_workflow.mu_dc import process_mu_delayed_coinc

        if os.path.exists(output[0]):
            print(f"Output file {output[0]} already exists. Skipping processing.")
            return
        process_mu_delayed_coinc(
            input, 
            output,
            default_ref_version=config["default_ref_version"],
            fallback_ref_version=config["fallback_ref_version"],
            metadata=config["metadata"],
            max_delta_sec=config["selection"]["pre"]["max_delta_sec"],
            min_cuspEmax=config["selection"]["pre"]["e_min"]
        )

