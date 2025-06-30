# rules/calculate_exposure.smk
rule calculate_ideal_main_selection:
    input: 
        rc_rates="gen/rc_bkg/rates_bb_like.yaml",
        sel_eff="gen/sim/selection_eff.yaml"
    output:
        "gen/ideal_selection.yaml"
    run:
        from ge77m_search_workflow.ideal_selection import ideal_selection

        ideal_selection(
            rc_rate_file=input[0],
            selection_eff_file=input[1],
            output_file=output[0]
        )




