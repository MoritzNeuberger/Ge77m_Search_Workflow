rule make_skm_mgc:
    """
    Collect all mgc.lh5 files from the mu_hpge_coinc directory and concatenate them into skm_mgc.lh5.
    """

    input:
        barrier="/tmp/mgc_all.done"
    output:
        "gen/mu_hpge_coinc/skm_mgc.lh5"
    run:
        import lgdo.lh5 as lh5
        import lgdo.types as types
        import awkward as ak
        from tqdm import tqdm
        from ge77m_search_workflow import utils as ut

        mgc_files = sorted(glob.glob("gen/mu_hpge_coinc/*/*/*.lh5"))
        concatenated_data = []

        for file in tqdm(mgc_files):
            data = lh5.read_as("mgc", file, "ak")
            concatenated_data.append(data)

        #final_data = types.Table(ut.dict_to_lgdo(ak.concatenate(concatenated_data).to_list()))
        final_data = types.Table(ak.concatenate(concatenated_data))
        lh5.write(final_data, name="mgc", lh5_file=str(output))
