# rules/select_candidates.smk

import numpy as np

def select_topology_candidates(skm):
    mask = (skm["delayed"]["coinc"]["multiplicity"] == 1 
        & ~(skm["delayed"]["coinc"]["spm_coinc"]) 
        & skm["delayed"]["geds"]["quality"]["psd_is_bb_like"]
    )
    return skm[mask]

def select_energy_candidates(skm, selection):

    mask = np.zeros(len(skm), dtype=bool)

    for ran in selection["delayed_e_range"]:
        mask = (mask + (
            (skm["delayed"]["geds"]["energy"] > ran[0]) &
            (skm["delayed"]["geds"]["energy"] < ran[1])
        ) ) > 0

    return skm[mask]

rule select_candidates:
    input:
        "gen/skm/skm_bb_like.lh5",
        "gen/rc_bkg/rates_bb_like.yaml"
    output:
        "gen/selected_candidates.lh5"
    run:
        import lgdo.lh5 as lh5
        import lgdo.types as types
        import matplotlib.pyplot as plt
        import numpy as np
        import awkward as ak 
        import yaml

        skm = lh5.read_as("mdc",input[0],"ak")
        selection = config["selection"]["main"]

        with open(input[1], 'r') as f:
            rc_bkg = yaml.safe_load(f)
        relative_fraction = rc_bkg["scaling"]

        # select time in on
        on_dataset = selection["on_dataset"]
        mask_on_data = (skm["dT"] > on_dataset[0]) & (skm["dT"] < on_dataset[1])
        skm_on = skm[mask_on_data]
        skm_on_topology = select_topology_candidates(skm_on)

        # select time in off 
        off_dataset = selection["off_dataset"]
        mask_off_data = (skm["dT"] > off_dataset[0]) & (skm["dT"] < off_dataset[1])
        skm_off = skm[mask_off_data]
        skm_off_topology = select_topology_candidates(skm_off)


        skm_on_selected = select_energy_candidates(skm_on_topology, selection)
        
        print(skm_on_selected)

        output = lh5.write(
            types.Table(skm_on_selected),name="mdc", lh5_file=str(output[0])
        )

        # Plotting the energy distribution of the selected candidates on and off (off scaled with relative_fraction)

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(
            skm_on_topology["delayed"]["geds"]["energy"],
            bins=np.linspace(0, 3000, 151),
            label="On candidates",
            histtype="step",
            color="black"
        )
        ax.hist(
            skm_off_topology["delayed"]["geds"]["energy"],
            bins=np.linspace(0, 3000, 151),
            label="Off candidates (scaled)",
            histtype="step",
            color="tab:orange",
            weights=np.full(len(skm_off_topology["delayed"]["geds"]["energy"]), relative_fraction)
        )   
        for ran in selection["delayed_e_range"]:
            ax.axvspan(ran[0], ran[1], color="grey", alpha=0.3, label=f"Selection range {ran[0]}-{ran[1]} keV")
        ax.set_xlabel("Energy (keV)")
        ax.set_ylabel("Counts")
        ax.set_title("Energy distribution of selected candidates")
        ax.set_yscale("log")
        ax.legend()
        plt.savefig("gen/figs/selected_candidates_energy_distribution.pdf", bbox_inches="tight")
        plt.close()
