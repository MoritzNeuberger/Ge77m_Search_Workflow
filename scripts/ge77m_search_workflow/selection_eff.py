import yaml
import h5py
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
from tqdm import tqdm
import scipy as sp

def load_sim_file(sim_file):
    with h5py.File(sim_file, "r") as f:
        group = f["awkward"]
        reconstituted = ak.from_buffers(
            ak.forms.from_json(group.attrs["form"]),
            group.attrs["length"],
            {k: np.asarray(v) for k, v in group.items()},
        )
    return ak.Array(reconstituted)

def energy_mask(e, e_range,scaling = 1e-3):
    tmp_arrays = None
    for e_min, e_max in e_range:
        if tmp_arrays is None:
            tmp_arrays = (e > e_min*scaling) & (e < e_max*scaling)
        else:
            tmp_arrays = tmp_arrays | ((e > e_min*scaling) & (e < e_max*scaling))
    return tmp_arrays

def calculate_selection_eff(
            sim_file_prompt,
            sim_file_delayed,
            output_file,
            selection
            ):
    
    data_sim_prompt = load_sim_file(sim_file_prompt)
    data_sim_delayed = load_sim_file(sim_file_delayed)

    # prepare prompt data
    mask_time = data_sim_prompt["ged_w_t"] < 3e6
    mask_status = data_sim_prompt["ged_status"] == "on" 
    mask_ge77 = data_sim_prompt["ged_contains_Ge77_gs"] | data_sim_prompt["ged_contains_Ge77_is"]
    mask_wf = data_sim_prompt["ged_w_t"] < 54e3

    energies_prompt = ak.ravel(data_sim_prompt["ged_etot_a"][mask_wf & mask_status & mask_ge77]) * 1e3

    mask_energy = (energies_prompt > selection["main"]["prompt_e_range"][0]) & (energies_prompt < selection["main"]["prompt_e_range"][1])

    n_ge77_selected = ak.sum(mask_energy)
    n_ge77_total = ak.sum(mask_ge77 & mask_status & mask_time)

    eps_prompt = float(n_ge77_selected / n_ge77_total) if n_ge77_total > 0 else 0
    unc_prompt = float(eps_prompt * np.sqrt(1 / n_ge77_selected + 1 / n_ge77_total * eps_prompt)) if n_ge77_selected > 0 else 0

    # draw efficiency dependence on lower energy threshold
    x_p = np.linspace(0,selection["main"]["prompt_e_range"][1],selection["main"]["prompt_e_range"][1] // 50 + 1)
    y_p = np.zeros_like(x_p)
    for i, e in enumerate(x_p):
        mask_energy = (energies_prompt > e) & (energies_prompt < selection["main"]["prompt_e_range"][1])
        n_ge77_selected = ak.sum(mask_energy)
        y_p[i] = n_ge77_selected / n_ge77_total if n_ge77_total > 0 else 0
    plt.plot(x_p, y_p, label="Prompt Efficiency")
    plt.axvline(selection["main"]["prompt_e_range"][0], color='red', linestyle='--', label="Lower Energy Threshold")
    plt.axvline(selection["main"]["prompt_e_range"][1], color='blue', linestyle='--', label="Upper Energy Threshold")
    plt.xlabel("Lower Energy Threshold (keV)")
    plt.ylabel("Efficiency")
    plt.savefig("./gen/figs/efficiency_prompt_over_lower_thr.pdf", bbox_inches='tight')
    plt.close()


    # prepare delayed data
    mask_ge77 = (data_sim_delayed["ged_w_t"] == 0)
    n_ge77 = len(data_sim_delayed)

    masks = {
        "energy": energy_mask(data_sim_delayed["ged_etot_a"], selection["main"]["delayed_e_range"]),
        "multiplicity": data_sim_delayed["ged_multiplicity_veto"],
        "lar": data_sim_delayed["lar_veto"],
        "psd": data_sim_delayed["MSE_veto"]
    }

    combined_masks = {
        "all": mask_ge77,
        "mult": mask_ge77 & masks["multiplicity"],
        "mult_lar": mask_ge77  & masks["multiplicity"] & masks["lar"],
        "mult_lar_psd": mask_ge77 & masks["multiplicity"] & masks["lar"] & masks["psd"]
    }

    n_ge77_selected = {
        key: ak.sum(masks["energy"] & combined_masks[key]) for key in combined_masks
    }

    eps_delayed_energs = {
        "val": float(n_ge77_selected["all"] / n_ge77) if n_ge77 > 0 else 0,
        "unc": float(np.sqrt(n_ge77_selected["all"] / n_ge77**2 + n_ge77_selected["all"]**2 / n_ge77**3) if n_ge77_selected["all"] > 0 else 0)
    }
    eps_delayed_analysis_cuts = {
        key: float(n_ge77_selected[key] / n_ge77_selected["all"]) for key in n_ge77_selected
    }

    additional_unc = {
        "all": 0,
        "mult": 0,
        "mult_lar": 0.05,
        "mult_lar_psd": 0.07
    }
    unc_delayed_analysis_cuts = {
        key: float(np.sqrt(n_ge77_selected[key]/n_ge77_selected["all"]**2 + n_ge77_selected["all"]**2/n_ge77_selected["all"]**3 + additional_unc[key]**2) if n_ge77_selected[key] > 0 else 0)
        for key in n_ge77_selected
    }

    # draw efficiency dependence on upper and lower energy threshold
    e_low = np.linspace(1500, 2700, (2700 - 1500) // 100 + 1)
    e_high = np.linspace(1500, 2700, (2700 - 1500) // 100 + 1)
    x, y = np.meshgrid(e_low, e_high)
    scans = {}
    for key, mask in combined_masks.items():
        z = np.zeros_like(x)
        for i in tqdm(range(len(e_low))):
            for j in range(len(e_high)):
                if e_low[i] < e_high[j]:
                    mask_energy = energy_mask(data_sim_delayed["ged_etot_a"], [(e_low[i], e_high[j])])
                    n_ge77_selected = ak.sum(mask_energy & combined_masks[key])
                    z[i, j] = float(n_ge77_selected / n_ge77) if n_ge77 > 0 else 0
        scans[key] = z.tolist()
    scans["x"] = x.tolist()
    scans["y"] = y.tolist()

    plt.pcolormesh(x, y, z, shading='auto')
    plt.colorbar(label="Efficiency")
    plt.xlabel("Lower Energy Threshold (keV)")
    plt.ylabel("Upper Energy Threshold (keV)")
    plt.title("Efficiency Dependence on Energy Thresholds")
    plt.savefig("./gen/figs/efficiency_delayed_over_energy_thr.pdf", bbox_inches='tight')
    plt.close()

    # time selection
    half_life = 53.7 # sec
    lifetime = half_life / np.log(2)

    eps_time = float(sp.stats.expon.cdf(selection["main"]["on_dataset"][1],0,lifetime))

    # save results
    results = {
        "eps_prompt": {
            "val": eps_prompt,
            "unc": unc_prompt,
            "scan": {
                "x": x_p.tolist(),
                "y": y_p.tolist()
            }
        },
        "eps_delayed": {
            "energy": eps_delayed_energs,
            "selection": {
                key: {
                    "val": eps_delayed_analysis_cuts[key],
                    "unc": unc_delayed_analysis_cuts[key]
                } for key in eps_delayed_analysis_cuts.keys()
            },
            "scan": scans
        },
        "eps_time": {
            "val": eps_time,
            "unc": 0  # no uncertainty for time selection
        },
        "selection": selection
    }

    with open(output_file, "w") as f:
        yaml.dump(results, f)