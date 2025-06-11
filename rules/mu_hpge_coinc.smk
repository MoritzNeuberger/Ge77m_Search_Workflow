# rules/mu_hpge_coinc.smk

rule mu_hpge_coinc:
    """
    Take .lh5 file from input_root and produce
    tier_mgc.lh5 in gen/mu_hpge_coinc/...
    """

    # declare Conda environment for reproducibility
    conda:
        "../envs/mu_hpge_coinc.yaml"

    # input is the previous-layer file
    input:
        # e.g. /…/phy/p03/r000/l200-p03-r000-phy-…-tier_pht.lh5
        lambda wc:  config["input_root"] + "/{lvl1}/{lvl2}/{base}.lh5".format(**wc)
    output:
        # e.g. gen/mu_hpge_coinc/p03/r000/l200-p03-r000-phy-…-tier_mgc.lh5
        os.path.join(
            config["out_root"],
            "mu_hpge_coinc",
            "{lvl1}",
            "{lvl2}",
            "{base}".replace("tier_pht", "tier_mgc") + ".lh5"
        )
    params:
        lvl1="{lvl1}",
        lvl2="{lvl2}",
        base="{base}"
    threads: 1
    run:
        # you would import and call your function here
        from ge77m_search_workflow.mu_hpge import process_mu_hpge_coinc
        process_mu_hpge_coinc(input[0], output[0])
