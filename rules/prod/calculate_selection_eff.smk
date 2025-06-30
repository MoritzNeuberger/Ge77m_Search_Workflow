# rules/calculate_exposure.smk
rule calculate_selection_eff:
    input: 
        "input/sim/sim_prompt.hf5",
        "input/sim/sim_delayed.hf5"
    output:
        "gen/sim/selection_eff.yaml",
    run:
        from ge77m_search_workflow.selection_eff import calculate_selection_eff

        calculate_selection_eff(
            sim_file_prompt=input[0],
            sim_file_delayed=input[1],
            output_file=output[0],
            selection=config["selection"]
        )

