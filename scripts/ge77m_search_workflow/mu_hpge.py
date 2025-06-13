import os
import numpy as np
from lgdo.lh5 import LH5Store, write
import lgdo.types as types 
import awkward as ak

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


import utils as ut


#columns = ["hit_table", "hit_idx", "evt_idx", "timestamp", "cuspEmax_ctc_cal", "is_good_hit","mu_diff","tp_01_hpge","tp_max_muon","is_in_coincidence_with_mu","is_saturated"]
acc_range = [-2000, 5000]
min_cuspEmax = 25

def get_pht_hpge(store, paths, selected_id, selected_idx):
    return store.read("ch{}/hit/".format(selected_id), paths["pht"], idx=[selected_idx])[0].view_as("ak")

def get_psp_hpge(store, paths, selected_id, selected_idx):
    return store.read("ch{}/dsp/".format(selected_id), paths["psp"], idx=[selected_idx])[0].view_as("ak")

def get_psp_muon(store, paths, selected_id, selected_idx):
    return store.read("ch{}/dsp/".format(selected_id), paths["psp"], idx=[selected_idx])[0].view_as("ak")

def fill_entry(output_data, selected_id, selected_idx, chmap, data_pht_hpge, data_psp_hpge, data_psp_muon, pet_data_geds, pet_data_trigger):
    evt_id = np.where(pet_data_geds["rawid"][selected_idx] == selected_id)[0][0]

    output_data["geds"]["energy"].append(data_pht_hpge["cuspEmax_ctc_cal"][0])
    output_data["geds"]["tp_01"].append(data_psp_hpge["tp_01"][0])
    output_data["geds"]["quality"]["is_bb_like"].append(pet_data_geds[selected_idx]["quality"]["is_bb_like"])
    output_data["geds"]["quality"]["is_good_channel"].append(pet_data_geds[selected_idx]["quality"]["is_good_channel"][evt_id])
    output_data["geds"]["quality"]["is_saturated"].append(data_pht_hpge["is_saturated"][0])
    output_data["geds"]["quality"]["channel_is_positive_polarity"].append(pet_data_geds[selected_idx]["quality"]["is_not_bb_like"]["is_pos_polarity_bits"] == 127)
    output_data["geds"]["id"]["hit_table"].append(selected_id)
    output_data["geds"]["id"]["hit_idx"].append(selected_idx)

    output_data["mu"]["tp_01"].append(data_psp_muon["tp_01"][0])
    output_data["mu"]["id"]["hit_table"].append(chmap.map("name")["MUON01"].daq.rawid)
    output_data["mu"]["id"]["hit_idx"].append(selected_idx)

    output_data["coinc"]["mu_diff"].append((data_psp_hpge["tp_01"] - data_psp_muon["tp_01"])[0])
    output_data["coinc"]["is_in_coincidence_with_mu"].append(
        (acc_range[0] < output_data["coinc"]["mu_diff"][-1] < acc_range[1])
    )
    output_data["coinc"]["id"]["evt_idx"].append(selected_idx)
    output_data["coinc"]["id"]["timestamp"].append(pet_data_trigger[selected_idx]["timestamp"])


def process_one_entry(selected_id, selected_idx, store, paths, chmap, pet_data_geds, pet_data_trigger, output_data):
    if not selected_id in pet_data_geds["rawid"][selected_idx]:
        print(f"Warning: Selected ID {selected_id} not found in pet_data_geds at index {selected_idx}. Skipping this entry.")
        return

    data_pht_hpge = get_pht_hpge(store, paths, selected_id, selected_idx)
    if data_pht_hpge["cuspEmax_ctc_cal"] < min_cuspEmax:
        return
    
    data_psp_hpge = get_psp_hpge(store, paths, selected_id, selected_idx)
    data_psp_muon = get_psp_muon(store, paths, selected_id, selected_idx)

    fill_entry(output_data, selected_id, selected_idx, chmap, data_pht_hpge, data_psp_hpge, data_psp_muon, pet_data_geds, pet_data_trigger)


def process_mu_hpge_coinc(input, output, default_ref_version="ref-v2.1.0", fallback_ref_version="ref-v2.0.0", metadata=None, raw_path=None):
    """
    Process coincidences between muon channel and HPGe channel.

    Args:
        input (list): List of input files.
        output (list): List of output files.

    This function reads the necessary data from the input files, processes it to find coincidences
    between muon and HPGe channels, and writes the results to the output files.
    
    """
    paths = ut.generate_paths_of_different_tiers_from_pht(str(input),default_ref_version=default_ref_version,fallback_ref_version=fallback_ref_version, raw_path=raw_path)

    store = LH5Store()
    pet_data_geds = store.read("/evt/geds/", paths["pet"])[0].view_as("ak")
    pet_data_coinc = store.read("/evt/coincident/", paths["pet"])[0].view_as("ak")
    pet_data_trigger = store.read("/evt/trigger/", paths["pet"])[0].view_as("ak")
    
    timestamp = ut.extract_timestamp_raw(paths["pet"])
    chmap = ut.generate_channel_map(timestamp, metadata=metadata)

    mask_muon_coinc = pet_data_coinc["muon"] & ~pet_data_coinc["puls"]  # & pet_data["geds"]["is_good_hit"]
    selected_idx = ak.to_numpy(ak.flatten(pet_data_geds["hit_idx"][mask_muon_coinc]))
    selected_id = ak.to_numpy(ak.flatten(pet_data_geds["rawid"][mask_muon_coinc]))


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
                "hit_idx": [],
            }
        },
        "mu": {
            "tp_01": [],
            "id": {
                "hit_table": [],
                "hit_idx": [],
            }
        },
        "coinc": {
            "mu_diff": [],
            "is_in_coincidence_with_mu": [],
            "id":{
                "evt_idx": [],
                "timestamp": []
            }
        }
    }

      
    for i in range(len(selected_idx)):
        process_one_entry(selected_id[i], selected_idx[i], store, paths, chmap, pet_data_geds, pet_data_trigger, output_data)

    output_lh5 = types.Table(col_dict=ut.dict_to_lgdo(output_data))
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