# rules/delayed_coinc.smk

rule delayed_coinc:
    """
    Take the mu_hpge_coinc output and produce
    tier_dc.lh5 in gen/delayed_coinc/...
    """

    wildcard_constraints:
        base=".*dc.*"

    input:
        # depends on mu_hpge_coinc output
        os.path.join(
            config["out_root"],
            "mu_hpge_coinc",
            "{lvl1}",
            "{lvl2}",
            "{base}".replace("tier_dc", "tier_mgc") + ".lh5"
        )
    output:
        # e.g. gen/delayed_coinc/p03/r000/l200-p03-r000-phy-…-tier_dc.lh5
        os.path.join(
            config["out_root"],
            "delayed_coinc",
            "{lvl1}",
            "{lvl2}",
            "{base}" + ".lh5"
        )
    run:
        import os
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        from ge77m_search_workflow.dc import process_delayed_coinc
        process_delayed_coinc(input[0], output[0])
