# rules/delayed_coinc.smk

configfile: "config.yaml"



def process_one_file(inp):
    import os
    import lgdo.lh5 as lh5
    import matplotlib.pyplot as plt
    import awkward as ak
    from ge77m_search_workflow import vis
    from glob import glob
    from ge77m_search_workflow import utils as ut
    from tqdm import tqdm
    import concurrent.futures
    import numpy as np

    low_energy_counter = 0

    lvl1 = inp.split("/")[-3]
    lvl2 = inp.split("/")[-2]

    os.makedirs("gen/figs/mu_hpge_coinc/{}/{}/".format(lvl1, lvl2), exist_ok=True)

    tmp = lh5.read_as("mgc", str(inp), "ak")
    path_pht = inp.replace("gen/mu_hpge_coinc",input_root).replace("mgc","pht")
    paths = ut.generate_paths_of_different_tiers_from_pht(path_pht,config["default_ref_version"], config["fallback_ref_version"], config["input_raw"])
    if len(tmp) == 0:
        return  # Skip empty files
    
    data_sorted_by_timestamp = tmp[ak.argsort(tmp["coinc"]["id"]["timestamp"])]
    lengths = np.diff(np.where(np.diff(data_sorted_by_timestamp["coinc"]["id"]["timestamp"], prepend=data_sorted_by_timestamp["coinc"]["id"]["timestamp"][0],append=data_sorted_by_timestamp["coinc"]["id"]["timestamp"][-1]+1) > 0)[0],prepend=0)
    data_grouped = ak.unflatten(data_sorted_by_timestamp, lengths)

    for i, event in enumerate(data_grouped):

        event_id = i
        base = os.path.basename(inp).replace(".lh5", "")

        output_file = "gen/figs/mu_hpge_coinc/{}/{}/{}_{:.9f}.pdf".format(lvl1, lvl2, base, event[0]["coinc"]["id"]["timestamp"])

        if os.path.exists(output_file):
            continue  # Skip if the file already exists
        
        n_triggered_ged = len(event)

        # smartly divide the figure in n_triggered_ged subplots
        n_cols = 3
        n_rows = (n_triggered_ged + n_cols - 1) // n_cols  # Ceiling division
        fig, axs = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
        axs = axs.flatten()  # Flatten the array
        for j, ax in enumerate(axs):
            if j < n_triggered_ged:
                vis.draw_event_single(event[j], paths, ax=ax)
            else:
                ax.axis("off")
        plt.tight_layout()
        plt.savefig(output_file, bbox_inches="tight")
        plt.close(fig)

# Format input_root with the default_ref_version from config
input_root = config["input_root"].format(default_ref_version=config["default_ref_version"])

waveform_block = "/tmp/mgc_plots_all.done"

rule draw_mgc_candidates:    
    """
    Take the mu_hpge_coinc output and produce plots for each entry in each file meeting certain conditions.
    """

    input:
        skm_mgc="gen/mu_hpge_coinc/skm_mgc.lh5"
    output:
        # e.g. gen/mu_delayed_coinc/p03/r000/l200-p03-r000-phy-…-tier_mdc.lh5
        waveform_block
    run:

        import os
        import lgdo.lh5 as lh5
        import matplotlib.pyplot as plt
        import awkward as ak
        from ge77m_search_workflow import vis
        from glob import glob
        from ge77m_search_workflow import utils as ut
        from tqdm import tqdm
        import concurrent.futures
        import numpy as np

        files = glob("../Ge77m_Search_Workflow/gen/mu_hpge_coinc/*/*/*.lh5")

        with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:
            results = list(executor.map(process_one_file, files))        
        #for inp in tqdm(files, desc="Processing files"):
        #    process_one_file(inp)

        # create pdf combining all plots


        output_files = glob("gen/figs/mu_hpge_coinc/*/*/*.pdf")
        output_files_low_energy = np.sort(glob("gen/figs/mu_hpge_coinc/*/*/*_low_energy.pdf"))

        output_files = np.sort([f for f in output_files if not f.endswith("_low_energy.pdf")])

        if not os.path.exists("gen/figs/mu_hpge_coinc/summary.pdf"):
            vis.generate_summary_pdf(
                output_files,
                "gen/figs/mu_hpge_coinc/summary.pdf"
            )
        if not os.path.exists("gen/figs/mu_hpge_coinc/summary_low_energy.pdf"):
            vis.generate_summary_pdf(
                output_files_low_energy,
                "gen/figs/mu_hpge_coinc/summary_low_energy.pdf"
            )

        os.system(f"touch {waveform_block}")  # Create the barrier file to indicate completion

