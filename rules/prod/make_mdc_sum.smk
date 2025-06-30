def get_period_run_info(filename):
    """
    Extract period and run information from the filename.
    """
    parts = filename.split('/')
    period = parts[-3]  # e.g., p01
    run = parts[-2]     # e.g., r000
    return period, run


mdc_files = sorted(glob.glob("gen/mu_delayed_coinc/*/*/*.lh5"))

structured_files = {}

for file in mdc_files:
    p,r = get_period_run_info(file)
    if p not in structured_files:
        structured_files[p] = {}
    if r not in structured_files[p]:
        structured_files[p][r] = []
    structured_files[p][r].append(file)


sum_files = []
for p in structured_files.keys():
    for r in structured_files[p].keys():
        structured_files[p][r] = sorted(structured_files[p][r])
        concatenated_data = []
        output_file_name = f"gen/mu_delayed_coinc/{p}/sum_mdc_{r}.lh5"
        sum_files.append(output_file_name)


rule make_sum_mdc:
    """
    Collect all mdc.lh5 files from the mu_delayed_coinc directory and concatenate them into sum_mdc.lh5.
    """

    input:
        barrier="/tmp/mdc_all.done"
    output:
        "gen/mu_delayed_coinc/sum_mdc.lh5",
        *sum_files
    run:
        import lgdo.lh5 as lh5
        import lgdo.types as types
        import awkward as ak
        from tqdm import tqdm
        from ge77m_search_workflow import utils as ut

        mdc_files = sorted(glob.glob("gen/mu_delayed_coinc/*/*/*.lh5"))
        
        structured_files = {}

        for file in mdc_files:
            p,r = get_period_run_info(file)
            if p not in structured_files:
                structured_files[p] = {}
            if r not in structured_files[p]:
                structured_files[p][r] = []
            structured_files[p][r].append(file)


        sum_files = []
        for p in structured_files.keys():
            for r in structured_files[p].keys():
                structured_files[p][r] = sorted(structured_files[p][r])
                concatenated_data = []
                output_file_name = f"./gen/mu_delayed_coinc/{p}/sum_mdc_{r}.lh5"
                if os.path.exists(output_file_name):
                    print(f"Skipping {output_file_name} as it already exists.")
                    continue
                for file in tqdm(structured_files[p][r], desc=f"Processing {p}/{r}"):
                    data = lh5.read_as("mdc", file, "ak")
                    concatenated_data.append(data)
                final_data = types.Table(ak.concatenate(concatenated_data))
                lh5.write(final_data, name="mdc", lh5_file=str(output_file_name))
                sum_files.append(output_file_name)

        if os.path.exists("./gen/mu_delayed_coinc/sum_mdc.lh5"):
            print(f"Skipping {output} as it already exists.")
            return

        concatenated_data = []
        for file in tqdm(sum_files):
            data = lh5.read_as("mdc", file, "ak")
            concatenated_data.append(data)
        #final_data = types.Table(ut.dict_to_lgdo(ak.concatenate(concatenated_data).to_list()))
        final_data = types.Table(ak.concatenate(concatenated_data))
        lh5.write(final_data, name="mdc", lh5_file="gen/mu_delayed_coinc/sum_mdc.lh5")
