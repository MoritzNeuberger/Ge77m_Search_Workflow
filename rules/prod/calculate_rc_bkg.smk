# rules/calculate_exposure.smk
rule calculate_rc_bkg:
    input: 
        #os.path.join("gen/mu_delayed_coinc/sum_mdc.lh5"),
        os.path.join("gen/skm/skm_" + "{psd_ms_type}" + ".lh5"),
        os.path.join("gen/exposure/exposures_" + "{psd_ms_type}" + ".yaml")
    output:
        "gen/rc_bkg/rates_{psd_ms_type}.yaml"
    params:
        psd_ms_type="{psd_ms_type}"
    run:
        from ge77m_search_workflow.rc_bkg import calculate_rc_bkg

        output_distr = str(output[0].replace("rates_", "distr_").replace(".yaml", ".lh5"))

        print(input[0], input[1], output[0], output_distr)

        calculate_rc_bkg(
            skm_file=input[0],
            exposure_file=input[1],
            output_rates=str(output[0]),
            output_distr=output_distr,
            delayed_e_range=config["selection"]["main"]["delayed_e_range"],
            on_dataset=config["selection"]["main"]["on_dataset"],
            off_dataset=config["selection"]["main"]["off_dataset"],
            psd_ms_type=params.psd_ms_type,
        )

