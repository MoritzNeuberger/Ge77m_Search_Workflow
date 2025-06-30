# rules/delayed_coinc.smk

configfile: "config.yaml"

# Format input_root with the default_ref_version from config
input_root = config["input_root"].format(default_ref_version=config["default_ref_version"])


rule check_muon_tagging:    
    """
    Count the number of tagged muon per run via different methods and summarize to textfile.
    """

    output:
        # e.g. gen/mu_delayed_coinc/p03/r000/l200-p03-r000-phy-…-tier_mdc.lh5
        "gen/muon_tagging_stats.txt"
    run:

        import lgdo.lh5 as lh5
        import awkward as ak
        from glob import glob
        from tqdm import tqdm
        import json
        import numpy as np

        pet_files = glob(input_root.replace("pht","pet") + "/*.lh5")

        summary = []

        for pet_file in tqdm(pet_files):
            data = lh5.read_as("evt/",pet_file,"ak")

            n_mu = np.sum(data["coincident"]["muon"])
            n_mu_offline = np.sum(data["coincident"]["muon_offline"])
            n_mu_either = np.sum(data["coincident"]["muon"] | data["coincident"]["muon_offline"])

            summary.append(
                {
                    "runinfo": {"p": pet_file.rsplit("-")[-4], "r": pet_file.rsplit("-")[-3]},
                    "n_mu": n_mu,
                    "n_mu_offline": n_mu_offline,
                    "n_mu_either": n_mu_either,
                    "ratio_offline_either": n_mu_offline / n_mu_either if n_mu_either > 0 else 0,
                    "ratio_mu_either": n_mu / n_mu_either if n_mu_either > 0 else 0
                }
            )

        class NumpyEncoder(json.JSONEncoder):
            """ Special json encoder for numpy types """
            def default(self, obj):
                if isinstance(obj, np.integer):
                    return int(obj)
                elif isinstance(obj, np.floating):
                    return float(obj)
                elif isinstance(obj, np.ndarray):
                    return obj.tolist()
                return json.JSONEncoder.default(self, obj)

        with open("gen/muon_tagging_stats.txt", "w") as f:
            f.write(json.dumps(summary, indent=4, cls=NumpyEncoder))

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        x = [f"{s['runinfo']['p']}-{s['runinfo']['r']}" for s in summary]
        argsort = np.argsort(x)
        x = np.array(x)[argsort]
        y_mu = np.array([s["n_mu"] for s in summary])[argsort]
        y_mu_offline = np.array([s["n_mu_offline"] for s in summary])[argsort]
        y_mu_either = np.array([s["n_mu_either"] for s in summary])[argsort]
        ax.plot(x, y_mu, label="Muon", marker='o', color="black")
        ax.plot(x, y_mu_offline, label="Muon Offline", marker='o', color="blue")
        ax.plot(x, y_mu_either, label="Muon Either", marker='o', color="red")
        ax.set_xlabel("Run")
        ax.set_ylabel("Number of Muons")
        plt.savefig("gen/figs/muon_tagging_stats.pdf", bbox_inches="tight")
        plt.close()

        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        x = [f"{s['runinfo']['p']}-{s['runinfo']['r']}" for s in summary]
        argsort = np.argsort(x)
        x = np.array(x)[argsort]
        y_mu_frac = np.array([s["ratio_mu_either"] for s in summary])[argsort]
        y_mu_frac_offline = np.array([s["ratio_offline_either"] for s in summary])[argsort]
        ax.plot(x, y_mu_frac, label="Muon", marker='o', color="black")
        ax.plot(x, y_mu_frac_offline, label="Muon Offline", marker='o', color="blue")
        ax.set_xlabel("Run")
        ax.set_ylabel("Fraction of Muons")
        plt.savefig("gen/figs/muon_tagging_stats_frac.pdf", bbox_inches="tight")
        plt.close()


