# rules/delayed_coinc.smk

configfile: "config.yaml"

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

        files = glob("../Ge77m_Search_Workflow/gen/mu_hpge_coinc/*/*/*.lh5")

        for inp in tqdm(files):

            lvl1 = inp.split("/")[-3]
            lvl2 = inp.split("/")[-2]

            os.makedirs("gen/figs/mu_hpge_coinc/{}/{}/".format(lvl1, lvl2), exist_ok=True)

            tmp = lh5.read_as("mgc", str(inp), "ak")
            path_pht = inp.replace("gen/mu_hpge_coinc",input_root).replace("mgc","pht")
            paths = ut.generate_paths_of_different_tiers_from_pht(path_pht,config["default_ref_version"], config["fallback_ref_version"], config["input_raw"])
            if len(tmp) == 0:
                continue  # Skip empty files
            
            for i, event in enumerate(tmp):
                if event["coinc"]["mu_diff"] < -2000 or event["coinc"]["mu_diff"] > 5000:
                    continue  # Skip events outside the mu_diff range
                
                if event["geds"]["energy"] < 500:
                    continue  # Skip events with energy below threshold
                
                # Draw the event
                event_id = i
                base = os.path.basename(inp).replace(".lh5", "")
                
                output_file = "gen/figs/mu_hpge_coinc/{}/{}/{}_{}_{:.9f}.pdf".format(lvl1, lvl2, base, event["geds"]["id"]["hit_table"], event["coinc"]["id"]["timestamp"])
                if os.path.exists(output_file):
                    continue  # Skip if the file already exists
                    
                fig, ax = plt.subplots(figsize=(10, 6))
                vis.draw_event(event, paths, ax=ax)
                plt.savefig(output_file, bbox_inches="tight")
                plt.close(fig)


        os.system(f"touch {waveform_block}")  # Create the barrier file to indicate completion
