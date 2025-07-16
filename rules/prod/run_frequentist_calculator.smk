# rules/run_frequentist_calculator.smk
rule run_frequentist_calculator:
    input: 
        "gen/rc_bkg/rates_bb_like.yaml",
        "gen/sim/selection_eff.yaml",
        "gen/selected_candidates.lh5",
        "gen/exposure/exposures_bb_like.yaml",
    output:
        "gen/freq_calc/{state}-{sys}.lh5"
    params:
        n_threads=20,
        sample_size=50000,
        range_max=2,
        range_points=200,
        state="{state}",
        sys="{sys}"
    run:
        from ge77m_search_workflow.freq_calc import calculate_mle_and_ci

        if str(params.state) not in ["Ge77m", "Ge77_and_Ge77m"]:
            raise ValueError(f"Invalid state: {params.state}. Must be 'Ge77m' or 'Ge77_and_Ge77m'.")
        
        parameter = "r_" + str(params.state)

        if str(params.sys) not in ["w_sys", "wo_sys"]:
            raise ValueError(f"Invalid systematics: {params.sys}. Must be 'w_sys' or 'wo_sys'.")

        if str(params.sys) == "w_sys":
            fix_systematics = False
        else:
            fix_systematics = True

        if parameter == "r_Ge77_and_Ge77m":
            range_max = params.range_max * 2 
        else:
            range_max = params.range_max

        calculate_mle_and_ci(
            output[0],
            input[0],
            input[1],
            input[2],
            input[3],
            n_threads=params.n_threads,
            sample_size=params.sample_size,
            range_max=range_max,
            range_points=params.range_points,
            fix_systematics=fix_systematics,
            param_type=parameter
        )
