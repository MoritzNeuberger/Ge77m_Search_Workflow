# rules/calculate_exposure.smk
rule calculate_exposure:
    output:
        "gen/exposure/exposures_{psd_ms_type}.yaml"
    params:
        psd_ms_type="{psd_ms_type}"
    run:
        from ge77m_search_workflow.exposure import calculate_exposure

        if "-" in params.psd_ms_type:
            psd_ms_type = params.psd_ms_type.split("-")
        else:
            psd_ms_type = params.psd_ms_type

        print("psd_ms_type:", psd_ms_type)

        calculate_exposure(
            output_path=output[0],
            metadata_path=config["metadata"],
            psd_ms_type=psd_ms_type,
            ignore_periods=config["ignore_periods"],
            physics_dataset=config["physics_dataset"]
        )