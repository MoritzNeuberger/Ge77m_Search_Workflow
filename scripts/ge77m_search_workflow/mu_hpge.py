import os
import numpy as np
from lgdo.lh5 import LH5Store, write
import lgdo.types as types 
import awkward as ak

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


import utils as ut

acc_range = [-2000, 5000]


def get_pet_data(store, paths):
    tmp = store.read("/evt/geds/", paths["pet"])
    if isinstance(tmp, tuple):
        pet_data_geds = store.read("/evt/geds/", paths["pet"])[0].view_as("ak")
        pet_data_coinc = store.read("/evt/coincident/", paths["pet"])[0].view_as("ak")
        pet_data_trigger = store.read("/evt/trigger/", paths["pet"])[0].view_as("ak")
    else:
        pet_data_geds = store.read("/evt/geds/", paths["pet"]).view_as("ak")
        pet_data_coinc = store.read("/evt/coincident/", paths["pet"]).view_as("ak")
        pet_data_trigger = store.read("/evt/trigger/", paths["pet"]).view_as("ak")
    return pet_data_geds, pet_data_coinc, pet_data_trigger

def get_pht_hpge(store, paths, selected_id, selected_idx):
    tmp = store.read("ch{}/hit/".format(selected_id), paths["pht"], idx=[selected_idx])
    if isinstance(tmp, tuple):
        return tmp[0].view_as("ak")
    else:
        return tmp.view_as("ak")

def get_psp_hpge(store, paths, selected_id, selected_idx):
    tmp = store.read("ch{}/dsp/".format(selected_id), paths["psp"], idx=[selected_idx])
    if isinstance(tmp, tuple):
        return tmp[0].view_as("ak")
    else:
        return tmp.view_as("ak")

def get_psp_muon(store, paths, selected_id, selected_idx, chmap):
    tmp = store.read("ch{}/dsp/".format(chmap.map("name")["MUON01"].daq.rawid), paths["psp"], idx=[selected_idx])
    if isinstance(tmp, tuple):
        return tmp[0].view_as("ak")
    else:
        return tmp.view_as("ak")

def fill_entry(output_data, selected_id, selected_idx, chmap, data_pht_hpge, data_psp_hpge, data_psp_muon, pet_data_geds, pet_data_trigger, pet_data_coinc, paths):
    evt_id = np.where(pet_data_geds["rawid"][selected_idx] == selected_id)[0][0]

    output_data["geds"]["energy"].append(float(data_pht_hpge["cuspEmax_ctc_cal"][0]))
    output_data["geds"]["other_energy_estimators"]["trapEmax_ctc_cal"].append(float(data_pht_hpge["trapEmax_ctc_cal"][0]))
    output_data["geds"]["other_energy_estimators"]["cuspEmax_ctc_cal"].append(float(data_pht_hpge["cuspEmax_ctc_cal"][0]))
    output_data["geds"]["other_energy_estimators"]["zacEmax_ctc_cal"].append(float(data_pht_hpge["zacEmax_ctc_cal"][0]))
    output_data["geds"]["tp_01"].append(float(data_psp_hpge["tp_01"][0]))
    output_data["geds"]["tp_max"].append(float(data_psp_hpge["tp_max"][0]))
    output_data["geds"]["wf_max"].append(float(data_psp_hpge["wf_max"][0]))
    output_data["geds"]["tailEmax"].append(float(data_psp_hpge["tailEmax"][0]))
    output_data["geds"]["quality"]["is_bb_like"].append(bool(pet_data_geds[selected_idx]["quality"]["is_bb_like"]))
    output_data["geds"]["quality"]["is_good_channel"].append(bool(pet_data_geds[selected_idx]["quality"]["is_good_channel"][evt_id]))
    output_data["geds"]["quality"]["is_saturated"].append(bool(data_pht_hpge["is_saturated"][0]))
    if selected_id in pet_data_geds[selected_idx]["quality"]["is_not_bb_like"]["rawid"]:
        quality_idx = np.where(selected_id == pet_data_geds[selected_idx]["quality"]["is_not_bb_like"]["rawid"])[0][0]
        output_data["geds"]["quality"]["channel_is_positive_polarity"].append(bool(pet_data_geds[selected_idx]["quality"]["is_not_bb_like"]["is_pos_polarity_bits"][quality_idx] == 127))
    else:
        output_data["geds"]["quality"]["channel_is_positive_polarity"].append(True)
    output_data["geds"]["id"]["hit_table"].append(int(selected_id))
    output_data["geds"]["id"]["hit_idx"].append(int(selected_idx))

    output_data["mu"]["tp_max"].append(float(data_psp_muon["tp_max"][0]))
    output_data["mu"]["id"]["hit_table"].append(int(chmap.map("name")["MUON01"].daq.rawid))
    output_data["mu"]["id"]["hit_idx"].append(int(selected_idx))
    output_data["mu"]["coinc_flags"]["muon"].append(bool(pet_data_coinc[selected_idx]["muon"]))
    output_data["mu"]["coinc_flags"]["muon_offline"].append(bool(pet_data_coinc[selected_idx]["muon_offline"]))


    output_data["coinc"]["mu_diff"].append(float(data_psp_hpge["tp_01"][0] - data_psp_muon["tp_max"][0]))
    output_data["coinc"]["id"]["evt_idx"].append(int(selected_idx))
    output_data["coinc"]["id"]["timestamp"].append(float(pet_data_trigger[selected_idx]["timestamp"]))


def process_one_entry(selected_id, selected_idx, store, paths, chmap, pet_data_geds, pet_data_trigger, pet_data_coinc, output_data, min_cuspEmax):
    if not selected_id in pet_data_geds["rawid"][selected_idx]:
        print(f"Warning: Selected ID {selected_id} not found in pet_data_geds at index {selected_idx}. Skipping this entry.")
        return

    data_pht_hpge = get_pht_hpge(store, paths, selected_id, selected_idx)
    if data_pht_hpge["cuspEmax_ctc_cal"] < min_cuspEmax or np.isnan(data_pht_hpge["cuspEmax_ctc_cal"]):
        return
    
    data_psp_hpge = get_psp_hpge(store, paths, selected_id, selected_idx)
    data_psp_muon = get_psp_muon(store, paths, selected_id, selected_idx, chmap)

    fill_entry(output_data, selected_id, selected_idx, chmap, data_pht_hpge, data_psp_hpge, data_psp_muon, pet_data_geds, pet_data_trigger, pet_data_coinc, paths)

def make_selection(pet_data_coinc,pet_data_geds):
    mask = (pet_data_coinc["muon"] | pet_data_coinc["muon_offline"]) & ~pet_data_coinc["puls"]
    
    selection_idx = ak.to_numpy(ak.flatten(pet_data_geds["hit_idx"][mask]))
    selection_id = ak.to_numpy(ak.flatten(pet_data_geds["rawid"][mask]))
    return selection_idx, selection_id


def process_mu_hpge_coinc(input, output, default_ref_version="ref-v2.1.0", fallback_ref_version="ref-v2.0.0", metadata=None, input_raw=None, min_cuspEmax=25):
    """
    Process coincidences between muon channel and HPGe channel.

    Args:
        input (list): List of input files.
        output (list): List of output files.

    This function reads the necessary data from the input files, processes it to find coincidences
    between muon and HPGe channels, and writes the results to the output files.
    
    """
    paths = ut.generate_paths_of_different_tiers_from_pht(str(input),default_ref_version=default_ref_version,fallback_ref_version=fallback_ref_version, input_raw=input_raw)

    store = LH5Store()
    pet_data_geds, pet_data_coinc, pet_data_trigger = get_pet_data(store, paths)
    
    timestamp = ut.extract_timestamp_raw(paths["pet"])
    chmap = ut.generate_channel_map(timestamp, metadata=metadata)
    data_streams_hpge = ut.select_datastreams(chmap, "HPGE")

    selected_idx, selected_id = make_selection(pet_data_coinc, pet_data_geds)

    output_data = {
        "geds": {
            "energy": [],
            "tp_01": [],
            "quality": {
                "is_bb_like": [],
                "is_good_channel": [],
                "is_saturated": [],
                "channel_is_positive_polarity": []
            },
            "id": {
                "hit_table": [],
                "hit_idx": []
            },
            "other_energy_estimators": {
                "trapEmax_ctc_cal": [],
                "cuspEmax_ctc_cal": [],
                "zacEmax_ctc_cal": []
            },
            "tp_max": [],
            "wf_max": [],
            "tailEmax": []
        },
        "mu": {
            "tp_max": [],
            "id": {
                "hit_table": [],
                "hit_idx": [],
            },
            "coinc_flags":{
                "muon": [],
                "muon_offline": []
            }

        },
        "coinc": {
            "mu_diff": [],
            "id":{
                "evt_idx": [],
                "timestamp": []
            }
        }
    }

      
    for i in range(len(selected_idx)):
        if selected_id[i] not in data_streams_hpge:
            continue
        process_one_entry(selected_id[i], selected_idx[i], store, paths, chmap, pet_data_geds, pet_data_trigger, pet_data_coinc, output_data, min_cuspEmax)

    output_lh5 = types.Table(ut.enforce_type_prompt(ak.Array(output_data)))
    write(output_lh5, name="mgc", lh5_file=str(output))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process muon HPGe coincidences.")
    parser.add_argument("input", type=str, help="Input file path.")
    parser.add_argument("output", type=str, help="Output file path.")
    args = parser.parse_args()
    print("Processing muon HPGe coincidences...")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    process_mu_hpge_coinc(args.input, args.output)
    print("Processing completed.")