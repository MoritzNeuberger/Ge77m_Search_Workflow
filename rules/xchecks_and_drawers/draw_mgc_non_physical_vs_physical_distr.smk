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


rule draw_mgc_non_physical_vs_physical_distr:
    """
    Once sum_mgc.lh5 is created, this rule draws the mgc rate plot.
    """

    input:
        sum_mgc="gen/mu_hpge_coinc/sum_mgc.lh5",
        sim="input/sim/sim_all_muons_with_hpge.hf5"
    output:
        "gen/figs/mu_non_physical_vs_physical_distr_energy_lin.pdf"
    run:
        
        import lgdo.lh5 as lh5
        import numpy as np
        import matplotlib.pyplot as plt
        import lgdo.types as types
        from ge77m_search_workflow import utils as ut
        # Read the data from the LH5 file
        tmp = lh5.read_as("mgc", str(input.sum_mgc), "ak")

        physics_mu_diff = [-1250,-400]
        non_physics_mu_diff = [-400, 1000]

        mask = tmp["mu"]["coinc_flags"]["muon"]


        sim = load_sim_file(input.sim)
        n_files = 998
        n_mu = n_files * 50_000
        time = n_mu/0.1248 / 60 / 60 / 24 / 365.25  # in years
        sim_exp = time * 140

        mask_time = sim["ged_w_t"] < 3e6
        sim_energies = np.array(ak.ravel(sim["ged_etot_a"][mask_time])) * 1e3
        mask_sim_energy = sim_energies > 25
        sim_energies = sim_energies[mask_sim_energy] 


        fig, ax = plt.subplots(figsize=(10, 6))

        ax.hist2d(np.array(tmp["coinc"]["mu_diff"][mask]), np.array(tmp["geds"]["energy"][mask]), 
           bins=[np.linspace(-2000,2000, 101), 10**np.linspace(1.1,4,101)], 
           cmin=1)
        for mu_diff in physics_mu_diff:
            ax.axvline(mu_diff, color="red", linestyle="--")
        for mu_diff in non_physics_mu_diff:
            ax.axvline(mu_diff, color="blue", linestyle="--")
        ax.set_yscale("log")
        ax.set_xlabel("HPGe - Muon triggers (ns)")
        ax.set_ylabel("Energy (keV)")
        plt.savefig("./gen/figs/mu_non_physical_vs_physical_distr_mu_diff_vs_energy.pdf", bbox_inches="tight")
        plt.close()

        mask_physics = (
            (tmp["coinc"]["mu_diff"][mask] > physics_mu_diff[0]) & 
            (tmp["coinc"]["mu_diff"][mask] < physics_mu_diff[1])
        )
        mask_non_physics = (
            (tmp["coinc"]["mu_diff"][mask] > non_physics_mu_diff[0]) & 
            (tmp["coinc"]["mu_diff"][mask] < non_physics_mu_diff[1])
        )

        energies_physics = np.array(tmp["geds"]["energy"][mask][mask_physics])
        energies_non_physics = np.array(tmp["geds"]["energy"][mask][mask_non_physics])

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(energies_physics, 
           bins=10**np.linspace(1,4,1001),
           label="Physics mu_diff",
           histtype="step")
        ax.hist(energies_non_physics, 
           bins=10**np.linspace(1,4,1001),
           label="Non-physics mu_diff",
           histtype="step")
        ax.set_xscale("log")
        ax.set_xlabel("Energy (keV)")
        ax.set_ylabel("Counts")
        ax.legend()
        plt.savefig("./gen/figs/mu_non_physical_vs_physical_distr_energy_log.pdf", bbox_inches="tight")
        plt.close()


        fig, ax = plt.subplots(2,1,figsize=(10, 6),height_ratios=(1, 0.33), sharex=True)
        ax[0].hist(energies_physics, 
           bins=np.linspace(0,1000,201),
           label="Physics mu_diff",
           histtype="step")
        ax[0].hist(energies_non_physics, 
           bins=np.linspace(0,1000,201),
           label="Non-physics mu_diff",
           histtype="step")
        ax[0].set_yscale("log")
        ax[0].set_xlabel("Energy (keV)")
        ax[0].set_ylabel("Counts")
        ax[0].legend()

        ax[1].hist(sim_energies, 
           bins=np.linspace(0,1000,201),
           label="Simulation",
           histtype="step", density=True,color="black")
        ax[1].set_xlabel("Energy (keV)")
        ax[1].set_ylabel("au")
        ax[1].set_yscale("log")
        ax[1].legend()

        plt.savefig("./gen/figs/mu_non_physical_vs_physical_distr_energy_lin.pdf", bbox_inches="tight")
        plt.close()