# rules/make_skm.smk

rule make_skm:
    """
    Collect all delayed_coinc .lh5 files under a specific lvl2 folder and merge into a single output file.
    Output: gen/make_skm/{lvl1}/{lvl2}.lh5
    """

    wildcard_constraints:
        base="*dc.*"
    input:
        lambda wc: sorted(
            glob.glob(
                os.path.join(
                    config["out_root"],
                    "delayed_coinc",
                    wc.lvl1,
                    wc.lvl2,
                    "*.lh5"
                )
            )
        )
    output:
        os.path.join(
            config["out_root"],
            "make_skm",
            "{lvl1}",
            "{lvl1}-{lvl2}.lh5"
        )
    params:
        lvl1="{lvl1}",
        lvl2="{lvl2}"
    run:
        import os
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        # Import and call your merge function here
        from ge77m_search_workflow.make_skm import merge_lh5_files
        merge_lh5_files(list(input), output[0])
