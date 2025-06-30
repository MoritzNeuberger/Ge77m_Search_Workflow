# rules/delayed_coinc.smk

configfile: "config.yaml"

# Format input_root with the default_ref_version from config
input_root = config["input_root"].format(default_ref_version=config["default_ref_version"])


rule check_saturation_sim_vs_data:    
    input: 
        data_path="gen/mu_hpge_coinc/sum_mgc.lh5",
        sim_path="input/sim/sim_all_muons_with_hpge.hf5"
    output:
        "gen/muon_saturation_sim_vs_data.yaml"
    run:
        import lgdo.lh5 as lh5
        import lgdo.types as types
        import h5py
        import numpy as np
        import awkward as ak
        from ge77m_search_workflow import utils as ut
        import yaml
        # Read the data from the LH5 file
        data = lh5.read_as("mgc", str(input.data_path), "ak")
        data_energies = ak.ravel(data["geds"]["energy"])

        with h5py.File(str(input.sim_path), "r") as f:
            group = f["awkward"]
            reconstituted = ak.from_buffers(
                ak.forms.from_json(group.attrs["form"]),
                group.attrs["length"],
                {k: np.asarray(v) for k, v in group.items()},
            )
        sim = ak.Array(reconstituted)

        mask_time = sim["ged_w_t"] < 3e6
        sim_energies = np.array(ak.ravel(sim["ged_etot_a"][mask_time])) * 1e3

        physics_energy_range = [3000,6000]
        saturated_energy_range = [6500, float(np.max([np.max(data_energies), np.max(sim_energies)]))]

        n_d_physics = np.sum((data_energies > physics_energy_range[0]) & (data_energies < physics_energy_range[1]))
        n_s_physics = np.sum((sim_energies > physics_energy_range[0]) & (sim_energies < physics_energy_range[1]))

        n_d_saturated = np.sum((data_energies > saturated_energy_range[0]) & (data_energies < saturated_energy_range[1]))
        n_s_saturated = np.sum((sim_energies > saturated_energy_range[0]) & (sim_energies < saturated_energy_range[1]))

        r_d = n_d_physics / n_d_saturated
        unc_r_d = np.sqrt( n_d_physics / n_d_saturated**2 + n_d_saturated * n_d_physics**2 / n_d_saturated**4 )
        r_s = n_s_physics / n_s_saturated
        unc_r_s = r_s * np.sqrt( 1 / n_s_physics + 1 / n_s_saturated * r_s )
        print(n_s_physics,n_s_saturated,unc_r_s,n_s_physics / n_s_saturated**2 ,n_s_saturated * n_s_physics**2 / n_s_saturated**4 )

        output_data = {
            "data": {
                "n_physics": int(n_d_physics),
                "n_saturated": int(n_d_saturated),
                "r": float(r_d),
                "unc_r": float(unc_r_d),
            },
            "sim": {
                "n_physics": int(n_s_physics),
                "n_saturated": int(n_s_saturated),
                "r": float(r_s),
                "unc_r": float(unc_r_s),
            },
            "physics_energy_range": physics_energy_range,
            "saturated_energy_range": saturated_energy_range
        }

        with open(str(output[0]), "w") as f:
            yaml.dump(output_data, f)

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(2,1,figsize=(8, 6))

        data_energies_squashed = np.array(data_energies)
        data_energies_squashed[data_energies_squashed > 6500] = 7000

        sim_energies_squashed = np.array(sim_energies)
        sim_energies_squashed[sim_energies_squashed > 6500] = 7000

        ax[0].hist(data_energies_squashed, bins=np.linspace(3000, 7500, 101), label="Data", histtype="step", lw=2)

        ax[1].hist(sim_energies_squashed, bins=np.linspace(3000, 7500, 101), label="Sim", histtype="step", lw=2)

        fig.legend()
        ax[1].set_xlabel("Energy [keV]")

        ax[0].set_ylabel("Counts")
        ax[1].set_ylabel("Counts")

        ax[0].set_yscale("log")
        ax[1].set_yscale("log")

        plt.savefig("./gen/figs/muon_saturation_sim_vs_data.pdf", bbox_inches="tight")