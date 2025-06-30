import h5py
import numpy as np
import awkward as ak

def load_sim_file(sim_file):
    with h5py.File(sim_file, "r") as f:
        group = f["awkward"]
        reconstituted = ak.from_buffers(
            ak.forms.from_json(group.attrs["form"]),
            group.attrs["length"],
            {k: np.asarray(v) for k, v in group.items()},
        )
    return ak.Array(reconstituted)

rule draw_mgc_data_vs_sim:
    """
    Once sum_mgc.lh5 is created, this rule draws the mgc rate plot.
    """

    input:
        sum_mgc="gen/mu_hpge_coinc/sum_mgc.lh5",
        sim="input/sim/sim_all_muons_with_hpge.hf5",
        exposure="gen/exposure/exposures_any.yaml"
    output:
        "gen/figs/mgc_data_vs_sim_energy_log.pdf"
    run:
        
        import lgdo.lh5 as lh5
        import numpy as np
        import matplotlib.pyplot as plt
        import lgdo.types as types
        from ge77m_search_workflow import utils as ut
        import yaml
        # Read the data from the LH5 file
        data = lh5.read_as("mgc", str(input.sum_mgc), "ak")
        mask_is_good = data["geds"]["quality"]["is_good_channel"]
        data_energies = np.array(ak.ravel(data["geds"]["energy"][mask_is_good]))

        # Filter out energies below 25 keV
        mask_data_energy = data_energies > 25
        data_energies = data_energies[mask_data_energy]

        # project all energies > 6.5 MeV to 7 MeV
        data_energies = np.where(data_energies > 6.5e3, 7e3, data_energies)

        with open(input.exposure, "r") as f:
            exposure = yaml.safe_load(f) 

        sim = load_sim_file(input.sim)
        n_files = 998
        n_mu = n_files * 50_000
        time = n_mu/0.1248 / 60 / 60 / 24 / 365.25  # in years
        sim_exp = time * 140

        mask_time = sim["ged_w_t"] < 3e6
        sim_energies = np.array(ak.ravel(sim["ged_etot_a"][mask_time])) * 1e3
        mask_sim_energy = sim_energies > 25
        sim_energies = sim_energies[mask_sim_energy] 

        # project all energies > 6.5 MeV to 7 MeV
        sim_energies = np.where(sim_energies > 6.5e3, 7e3, sim_energies)

        fig, ax = plt.subplots(2,1,figsize=(10, 6), height_ratios=(1, 0.33), sharex=True)
        bins = 10**np.linspace(1.4,4,401)
        #bin_width = bins[1] - bins[0]
        d_a,d_b,_ = ax[0].hist(data_energies,
            bins=bins,
            label="Data",
            histtype="step",
            weights=np.full(len(data_energies), 1.0 / exposure["total_exposure"]))
        s_a,s_b,_ = ax[0].hist(sim_energies,
            bins=bins,
            label="Sim",
            histtype="step",
            weights=np.full(len(sim_energies), 1.0 / sim_exp))

        diff = d_a - s_a
        ax[0].step(bins[:-1], diff, where="post", label="Difference", color="black")
        ax[0].axvline(500, color="red", linestyle="--", label="500 keV threshold")

        ax[0].set_xscale("log")
        ax[0].set_yscale("log")
        ax[0].set_ylim(1e-3, 2e2)
        #ax[0].set_xlabel("Energy (keV)")
        ax[0].set_ylabel("Counts / (kg yr)")
        ax[0].legend()

        # Draw the ratio plot
        ratio = np.divide(d_a, s_a, out=np.zeros_like(d_a), where=s_a!=0)
        ax[1].step(bins[:-1], ratio, where="post", color="black")
        ax[1].axhline(1, color="red", linestyle="--")
        ax[1].set_xlabel("Energy (keV)")
        ax[1].set_ylabel("Data / Sim")
        ax[1].set_xscale("log")
        ax[1].set_yscale("log")

        plt.savefig("./gen/figs/mgc_data_vs_sim_energy_log.pdf", bbox_inches="tight")
