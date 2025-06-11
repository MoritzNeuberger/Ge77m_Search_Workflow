# rules/mu_hpge_coinc.smk

configfile: "config.yaml"

# Format input_root with the default_ref_version from config
input_root = config["input_root"].format(default_ref_version=config["default_ref_version"])

rule mu_hpge_coinc:
    """
    Take .lh5 file from input_root and produce
    tier_mgc.lh5 in gen/mu_hpge_coinc/...
    """

    wildcard_constraints:
        base=".*mgc.*"

    # input is the previous-layer file
    input:
        lambda wc: input_root + "/{lvl1}/{lvl2}/{base}.lh5".format(
            lvl1=wc["lvl1"],
            lvl2=wc["lvl2"],
            base=wc["base"].replace("tier_mgc", "tier_pht")
        )
    output:
        # e.g. gen/mu_hpge_coinc/p03/r000/l200-p03-r000-phy-…-tier_mgc.lh5
        os.path.join(
            config["out_root"],
            "mu_hpge_coinc",
            "{lvl1}",
            "{lvl2}",
            "{base}" + ".lh5"
        )
    params:
        lvl1="{lvl1}",
        lvl2="{lvl2}",
        base="{base}"
    run:
        import os
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        from ge77m_search_workflow.mu_hpge import process_mu_hpge_coinc
        process_mu_hpge_coinc(input, output, default_ref_version=config["default_ref_version"], fallback_ref_version=config["fallback_ref_version"])
