import glob
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils as ut
import pandas as pd
import pygama
from packaging.version import Version
from pygama import logging
from pygama.flow import DataLoader
from legendmeta import LegendMetadata
import lgdo.lh5 as lh5
import awkward as ak
import numpy as np
from tqdm import tqdm
import lgdo.types as types 

def load_mgc_data_of_run(mgc_sum_file):
    return lh5.read_as("mgc", mgc_sum_file,"ak")

def get_pht_hpge(store, paths, selected_id, selected_idx):
    if isinstance(selected_idx, int):
        tmp = store.read("ch{}/hit/".format(selected_id), paths["pht"], idx=[selected_idx])
        if isinstance(tmp, tuple):
            return tmp[0].view_as("ak")
        else:
            return tmp.view_as("ak")
    elif isinstance(selected_idx, np.ndarray):
        tmp = store.read("ch{}/hit/".format(selected_id), paths["pht"], idx=selected_idx)
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

def get_pet_data(store, paths):
    tmp = store.read("/evt/geds/", paths["pet"])
    if isinstance(tmp, tuple):
        pet_data_geds = store.read("/evt/geds/", paths["pet"])[0].view_as("ak")
        pet_data_coinc = store.read("/evt/coincident/", paths["pet"])[0].view_as("ak")
    else:
        pet_data_geds = store.read("/evt/geds/", paths["pet"]).view_as("ak")
        pet_data_coinc = store.read("/evt/coincident/", paths["pet"]).view_as("ak")
    return pet_data_geds, pet_data_coinc

def generate_hpge_data(pht_file_name, default_ref_version, fallback_ref_version, metadata=None, min_cuspEmax=25):
    paths = ut.generate_paths_of_different_tiers_from_pht(
        pht_file_name, default_ref_version=default_ref_version, fallback_ref_version=fallback_ref_version
    )
    store = lh5.LH5Store()
    pet_data_geds, pet_data_coinc = get_pet_data(store, paths)
    timestamp = ut.extract_timestamp_raw(paths["pet"])
    chmap = ut.generate_channel_map(timestamp, metadata=metadata)
    data_streams_hpge = list(ut.select_datastreams(chmap, "HPGE"))

    mask = (
        pet_data_coinc["geds"] &
        ~pet_data_coinc["muon_offline"] &
        ~pet_data_coinc["puls"] &
        (pet_data_geds["energy_sum"] > 25)
    )


    selected_idx = ak.to_numpy(ak.flatten(pet_data_geds["hit_idx"][mask]))
    selected_id = ak.to_numpy(ak.flatten(pet_data_geds["rawid"][mask]))

    mask_selected_id = np.isin(selected_id, data_streams_hpge)
    selected_idx = selected_idx[mask_selected_id]
    selected_id = selected_id[mask_selected_id]

    # Preallocate output arrays
    unique_ids = np.unique(selected_id)
    n = np.max(selected_idx) + 1
    m = np.max(unique_ids)  + 1
    energy = np.empty((m,n), dtype=np.float64)
    quality_is_bb_like = np.empty((m,n), dtype=bool)
    psd_is_good = np.empty((m,n), dtype=bool)
    psd_is_bb_like = np.empty((m,n), dtype=bool)
    hit_table = np.empty((m,n), dtype=np.int64)
    hit_idx = np.empty((m,n), dtype=np.int64)
    timestamp_arr = np.empty((m,n), dtype=np.float64)
    multiplicity = np.empty((m,n), dtype=np.int64)
    spm_coinc = np.empty((m,n), dtype=bool)

    # Vectorized or batch processing here if possible
    valid_indices = np.zeros((m,n), dtype=bool)

    counter = -1
    for id in unique_ids:
        counter += 1
        mask_id = selected_id == id
        downselected_idx = selected_idx[mask_id]
        data_pht_hpge = get_pht_hpge(store, paths, id, downselected_idx)
        mask_cuspEmax = (
            (data_pht_hpge["cuspEmax_ctc_cal"] >= min_cuspEmax) &
            (~np.isnan(data_pht_hpge["cuspEmax_ctc_cal"]))
        )
        if not np.any(mask_cuspEmax):
            continue

        # Get the indices in the original arrays to update
        valid_idx = downselected_idx[mask_cuspEmax]

        
        evt_id = pet_data_geds[valid_idx]["rawid"] == id
        hit_table[counter][valid_idx] = id
        hit_idx[counter][valid_idx] = valid_idx
        timestamp_arr[counter][valid_idx] = data_pht_hpge["timestamp"][mask_cuspEmax]
        energy[counter][valid_idx] = data_pht_hpge["cuspEmax_ctc_cal"][mask_cuspEmax]
        quality_is_bb_like[counter][valid_idx] = pet_data_geds[valid_idx]["quality"]["is_bb_like"]
        psd_is_good[counter][valid_idx] = ak.ravel(pet_data_geds[valid_idx]["psd"]["is_good"])[ak.ravel(evt_id)]
        psd_is_bb_like[counter][valid_idx] = ak.ravel(pet_data_geds[valid_idx]["psd"]["is_bb_like"])[ak.ravel(evt_id)]
        multiplicity[counter][valid_idx] = pet_data_geds[valid_idx]["multiplicity"]
        spm_coinc[counter][valid_idx] = pet_data_coinc[valid_idx]["spms"]
        valid_indices[counter][valid_idx] = True

    # After the loop, slice arrays to only valid entries
    output_data = {
        "geds": {
            "energy": ak.ravel(energy[valid_indices]),
            "quality": {
                "quality_is_bb_like": ak.ravel(quality_is_bb_like[valid_indices]),
                "psd_is_good": ak.ravel(psd_is_good[valid_indices]),
                "psd_is_bb_like": ak.ravel(psd_is_bb_like[valid_indices]),
            },
            "id": {
                "hit_table": ak.ravel(hit_table[valid_indices]),
                "hit_idx": ak.ravel(hit_idx[valid_indices]),
                "timestamp": ak.ravel(timestamp_arr[valid_indices])
            }
        },
        "coinc": {
            "multiplicity": ak.ravel(multiplicity[valid_indices]),
            "spm_coinc":ak.ravel(spm_coinc[valid_indices])
        }
    }
    ak_output = ak.Array(output_data)
    argsort_id = ak.argsort(ak_output["geds"]["id"]["hit_table"])
    ak_output = ak_output[argsort_id]
    argsort_ts = ak.argsort(ak_output["geds"]["id"]["timestamp"])
    # sort by geds/id/timestamp
    return ak_output[argsort_ts]


def enforce_type(array):
    
    prompt = ut.enforce_type_prompt(array["prompt"])
    delayed = ut.enforce_type_delayed(array["delayed"])

    return  ak.Array({
        "dT": np.array(array["dT"]).astype("float64"),
        "prompt": prompt,
        "delayed": delayed
    })

def find_delayed_coincicence_candidates(hpge_data, mgc_data, max_delta_sec):
    # Merge the two dataframes on the "hit_table" column

    output_data = {
        "dT": [],
        "delayed": [],
        "prompt": []
    }


    for i in range(len(hpge_data)):
        #print(hpge_data[i])
        ts_hpge = hpge_data[i]["geds"]["id"]["timestamp"]
        dTs = ts_hpge - mgc_data["coinc"]["id"]["timestamp"]
        mask_dT = (dTs > 0) & (dTs <= max_delta_sec)
        mask_hit_table = (mgc_data["geds"]["id"]["hit_table"] == hpge_data[i]["geds"]["id"]["hit_table"])
        mask = mask_dT & mask_hit_table

        for j in range(len(mgc_data[mask])):

            dT = hpge_data[i]["geds"]["id"]["timestamp"] - mgc_data[mask][j]["coinc"]["id"]["timestamp"]
            output_data["prompt"].append(mgc_data[mask][j].to_list())
            output_data["delayed"].append(hpge_data[i].to_list())
            output_data["dT"].append(float(dT))

    return ak.Array(output_data)


def save_output(mdc_data,output):
    lh5.write(types.Table(enforce_type(mdc_data)), name="mdc", lh5_file=str(output))


def process_mu_delayed_coinc(input, output,default_ref_version="ref-v2.1.0", fallback_ref_version="ref-v2.0.0", max_delta_sec=100*53.7, metadata=None):
    import time

    start = time.time()
    mgc_data = load_mgc_data_of_run(input.mgc_sum_file)
    end = time.time()
    print(f"Loaded MGC data in {end - start:.2f} seconds")
    start = time.time()
    hpge_data = generate_hpge_data( input.pht_files, default_ref_version=default_ref_version, fallback_ref_version=fallback_ref_version, metadata=metadata)
    end = time.time()
    print(f"Generated HPGe data in {end - start:.2f} seconds")
    start = time.time()
    mdc_data = find_delayed_coincicence_candidates(hpge_data,mgc_data, max_delta_sec=max_delta_sec)
    end = time.time()
    print(f"Found delayed coincidence candidates in {end - start:.2f} seconds")
    save_output(mdc_data,output)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process muon delayed coincidence candidates.")
    parser.add_argument("mgc_sum_file", type=str, help="Path to the MGC sum file.")
    parser.add_argument("pht_files", type=str, help="Path to the PHT files.")
    parser.add_argument("output", type=str, help="Output file path.")
    args = parser.parse_args()

    process_mu_delayed_coinc(args, args.output)