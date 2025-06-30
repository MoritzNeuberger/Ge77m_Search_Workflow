# rules/calculate_mle_and_ci.smk
rule calculate_mle_and_ci:
    input:
        rc_rates="gen/rc_bkg/rates_bb_like.yaml",
        sel_eff="gen/sim/selection_eff.yaml",
        selected_candidates="gen/selected_candidates.lh5",
        exposure="gen/exposure/exposures_bb_like.yaml"
    output:
        "gen/mle_and_ci.yaml",
        "gen/figs/pval_distr.pdf"
    run:
        from ge77m_search_workflow.freq_calc import calculate_mle_and_ci

        calculate_mle_and_ci(
            output_path=output[0],
            rc_rates_path=input.rc_rates,
            sel_eff_path=input.sel_eff,
            selected_candidates_path=input.selected_candidates,
            exposure_path=input.exposure
        )