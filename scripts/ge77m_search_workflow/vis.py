import numpy as np
import matplotlib.pyplot as plt
import awkward as ak
import lgdo.lh5 as lh5 

def draw_waveform(hpge, muon, energy, hpge_rawid, mu_diff, ax=None):

    wf_p_ged = hpge["waveform_presummed"]
    wf_p_muon = muon["waveform_presummed"]

    x = np.arange(len(wf_p_ged["values"])) * wf_p_ged["dt"] * 1e-3  # Assuming 2 ns sampling rate
    y_ged = wf_p_ged["values"]
    y_mu = wf_p_muon["values"]
    
    if ax is None:
        fig, ax = plt.subplots()
        
    ax.plot(x, y_ged, label="ch{}".format(hpge_rawid) + " ({:.1f} keV, {:.1f}".format(energy, mu_diff) + r"$\mu$s)", color="black")
    ax.plot(x, y_mu, label="mu gate", color="tab:orange")
    ax.set_xlabel(r"Time ($\mu$s)")
    ax.set_ylabel("ADC")
    ax.set_xlim(50-10, 50+54)  # Adjust the x-axis limits as needed
    ax.legend()

def draw_event(entry, paths, ax=None):
    """
    Draws the event for a single mgc entry.
    Loads HPGe and muon waveforms using the read function.
    Assumes entry contains:
      - 'raw_path'
      - 'hpge_hit_table', 'hpge_hit_idx'
      - 'muon_hit_table', 'muon_hit_idx'
    """
    # Use the same read function as in draw_waveform
    hpge = lh5.read("ch{}/raw/".format(entry["geds"]["id"]['hit_table']), paths["raw"], idx=[entry["geds"]["id"]['hit_idx']]).view_as("ak")[0]
    muon = lh5.read("ch{}/raw/".format(entry["mu"]["id"]['hit_table']), paths["raw"], idx=[entry["mu"]["id"]['hit_idx']]).view_as("ak")[0]
    return draw_waveform(hpge, muon, entry["geds"]["energy"], entry["geds"]["id"]['hit_table'], entry["coinc"]["mu_diff"], ax=ax)

