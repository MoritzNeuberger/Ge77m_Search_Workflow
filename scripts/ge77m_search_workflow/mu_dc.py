
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

min_cuspEmax = 25

def load_mgc_data_of_run(mgc_skm_file):
    return lh5.read_as("mgc", mgc_skm_file,"ak")

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

def generate_hpge_data(pht_file_name, default_ref_version, fallback_ref_version, metadata=None):

    paths = ut.generate_paths_of_different_tiers_from_pht(pht_file_name,default_ref_version=default_ref_version,fallback_ref_version=fallback_ref_version)

    store = lh5.LH5Store()

    pet_data_geds, pet_data_coinc, _ = get_pet_data(store, paths)


    timestamp = ut.extract_timestamp_raw(paths["pet"])
    chmap = ut.generate_channel_map(timestamp, metadata=metadata)
    data_streams_hpge = ut.select_datastreams(chmap, "HPGE")

    mask = ~pet_data_coinc["muon_offline"] & ~pet_data_coinc["puls"] & (pet_data_geds["energy_sum"] > 50)
    
    selected_idx = ak.to_numpy(ak.flatten(pet_data_geds["hit_idx"][mask]))
    selected_id = ak.to_numpy(ak.flatten(pet_data_geds["rawid"][mask]))

    output_data = {
        "geds": {
            "energy": [],
            "quality": {
                "quality_is_bb_like": [],
                "psd_is_good": [],
                "psd_is_bb_like": [],
            },
            "id": {
                "hit_table": [],
                "hit_idx": [],
                "timestamp": []
            }
        },
        "coinc":{
            "multiplicity": [],
            "spm_coinc": []
        }
    }

    def _process_one_entry(idx, id, store, pet_data_geds, pet_data_coinc, output_data):
        if not id in pet_data_geds["rawid"][idx]:
            logging.warning("Selected ID {} not found in pet_data_geds".format(id))
            return

        data_pht_hpge = get_pht_hpge(store, paths, id, idx)
        if data_pht_hpge["cuspEmax_ctc_cal"] < min_cuspEmax or np.isnan(data_pht_hpge["cuspEmax_ctc_cal"]):
            return
        
        evt_id = np.where(pet_data_geds["rawid"][idx] == id)[0][0]
        
        output_data["geds"]["id"]["hit_table"].append(int(id))
        output_data["geds"]["id"]["hit_idx"].append(int(idx))
        output_data["geds"]["id"]["timestamp"].append(float(data_pht_hpge["timestamp"][0]))
        
        output_data["geds"]["energy"].append(float(data_pht_hpge["cuspEmax_ctc_cal"][0]))        
        output_data["geds"]["quality"]["quality_is_bb_like"].append(bool(pet_data_geds[idx]["quality"]["is_bb_like"]))
        output_data["geds"]["quality"]["psd_is_good"].append(bool(pet_data_geds[idx]["psd"]["is_good"][evt_id]))
        output_data["geds"]["quality"]["psd_is_bb_like"].append(bool(pet_data_geds[idx]["psd"]["is_bb_like"][evt_id]))
        
        output_data["coinc"]["multiplicity"].append(int(pet_data_geds[idx]["multiplicity"]))
        output_data["coinc"]["spm_coinc"].append(bool(pet_data_coinc[idx]["spms"]))


    for i in range(len(selected_idx)):
        if not selected_id[i] in data_streams_hpge:
            continue
        _process_one_entry(selected_idx[i], selected_id[i], store, pet_data_geds, pet_data_coinc, output_data)
        
    return ak.Array(output_data) #ut.dict_to_lgdo(output_data)


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

            dT = mgc_data[mask][j]["coinc"]["id"]["timestamp"] - hpge_data[i]["geds"]["id"]["timestamp"]
            output_data["prompt"].append(mgc_data[mask][j].to_list())
            output_data["delayed"].append(hpge_data[i].to_list())
            output_data["dT"].append(float(dT))

    return ak.Array(output_data)


def save_output(mdc_data,output):
    lh5.write(types.Table(enforce_type(mdc_data)), name="mdc", lh5_file=str(output))


def process_mu_delayed_coinc(input, output,default_ref_version="ref-v2.1.0", fallback_ref_version="ref-v2.0.0", max_delta_sec=100*53.7, metadata=None):
    import time

    start = time.time()
    mgc_data = load_mgc_data_of_run(input.mgc_skm_file)
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
