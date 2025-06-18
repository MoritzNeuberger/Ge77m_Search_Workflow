
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

def load_mgc_data_of_run(mgc_file_name):

    mgc_files = np.sort(glob.glob(mgc_file_name.rsplit("/",1)[0] + "*.lh5"))

    aggregated_data = []
    for path in mgc_files:
        data = lh5.read_as("mgc", path,"ak")
        aggregated_data.append(data)

    return ak.concatenate(aggregated_data)


def generate_hpge_data(pht_file_name, default_ref_version, fallback_ref_version, metadata=None):

    paths = ut.generate_paths_of_different_tiers_from_pht(pht_file_name,default_ref_version=default_ref_version,fallback_ref_version=fallback_ref_version)

    store = lh5.LH5Store()

    pet_data_geds = store.read("/evt/geds/", paths["pet"])[0].view_as("ak")
    pet_data_coinc = store.read("/evt/coincident/", paths["pet"])[0].view_as("ak")

    tcm_id = store.read("/hardware_tcm_1/array_id", paths["tcm"])[0].view_as("np")
    tcm_idx = store.read("/hardware_tcm_1/array_idx", paths["tcm"])[0].view_as("np")

    timestamp = ut.extract_timestamp_raw(paths["pet"])
    chmap = ut.generate_channel_map(timestamp, metadata=metadata)
    data_streams_hpge = ut.select_datastreams(chmap, "HPGE")

    mask = ~pet_data_coinc["muon_offline"] & ~pet_data_coinc["puls"] & (pet_data_geds["energy_sum"] > 50)
    
    selected_idx = ak.to_numpy(ak.flatten(pet_data_geds["hit_idx"][mask]))
    selected_id = ak.to_numpy(ak.flatten(pet_data_geds["rawid"][mask]))

    output_data = {
        "geds": {
            "energy": [],
            "tp_01": [],
            "quality": {
                "quality_is_bb_like": [],
                "psd_is_good": [],
                "psd_is_bb_like": [],
                "is_good_hit": [],
                "is_saturated": []
            },
            "id": {
                "hit_table": [],
                "hit_idx": [],
                "evt_idx": [],
                "timestamp": []
            }
        },
        "coinc":{
            "multiplicity": [],
            "spm_coinc": []
        }
    }
    
    for i in range(len(selected_idx)):
        
        hpge_idx = np.where((tcm_idx == selected_idx[i]) & (tcm_id == selected_id[i]))[0][0] #selected_idx[i]

        if tcm_id[hpge_idx] in data_streams_hpge:
        
            data_pht_hpge = store.read("ch{}/hit/".format(tcm_id[hpge_idx]),
                   paths["pht"],
                   idx=[tcm_idx[hpge_idx]])[0].view_as("ak")
            
            data_dsp_hpge = store.read("ch{}/dsp/".format(tcm_id[hpge_idx]),
                   paths["psp"],
                   idx=[tcm_idx[hpge_idx]])[0].view_as("ak")
            
            evt_idx = tcm_idx[hpge_idx]
            evt_id = np.where(pet_data_geds["rawid"][evt_idx] == selected_id[i])[0][0]
            
            output_data["geds"]["id"]["hit_table"].append(tcm_id[hpge_idx])
            output_data["geds"]["id"]["hit_idx"].append(hpge_idx)
            output_data["geds"]["id"]["evt_idx"].append(tcm_idx[hpge_idx])
            output_data["geds"]["id"]["timestamp"].append(data_pht_hpge["timestamp"][0])
            
            output_data["geds"]["energy"].append(data_pht_hpge["cuspEmax_ctc_cal"][0])
            output_data["geds"]["tp_01"].append(data_dsp_hpge["tp_01"][0])
            
            output_data["geds"]["quality"]["quality_is_bb_like"].append(pet_data_geds[evt_idx]["quality"]["is_bb_like"])
            output_data["geds"]["quality"]["psd_is_good"].append(pet_data_geds[evt_idx]["psd"]["is_good"][evt_id])
            output_data["geds"]["quality"]["psd_is_bb_like"].append(pet_data_geds[evt_idx]["psd"]["is_bb_like"][evt_id])
            
            output_data["coinc"]["multiplicity"].append(pet_data_geds[evt_idx]["multiplicity"])
            output_data["coinc"]["spm_coinc"].append(pet_data_coinc[evt_idx]["spms"])

    return ut.dict_to_lgdo(output_data)


def enforce_type(array):
    
    dT_type = ak.types.NumpyType("float64")  # dT

    prompt_type = ak.Array(array)["prompt"].type

    delayed_type = ak.types.RecordType(
        [
            ak.types.NumpyType("float64"),  # energy
            ak.types.NumpyType("float64"),  # tp_01
            ak.types.RecordType(
                [
                    ak.types.NumpyType("bool"),    # quality_is_bb_like
                    ak.types.NumpyType("bool"),    # psd_is_good
                    ak.types.NumpyType("bool"),    # psd_is_bb_like
                    ak.types.NumpyType("bool"),    # is_good_hit
                    ak.types.NumpyType("bool"),    # is_saturated
                ],
                ["quality_is_bb_like", "psd_is_good", "psd_is_bb_like", "is_good_hit", "is_saturated"]
            ),  # quality
            ak.types.RecordType(
                [
                    ak.types.NumpyType("int64"),   # hit_table
                    ak.types.NumpyType("int64"),   # hit_idx
                    ak.types.NumpyType("int64"),   # evt_idx
                    ak.types.NumpyType("float64"), # timestamp
                ],
                ["hit_table", "hit_idx", "evt_idx", "timestamp"]
            ),  # id
        ],
        ["energy", "tp_01", "quality", "id"]
    )

    output_type = ak.types.RecordType(
        [
            dT_type,
            prompt_type,
            delayed_type,
        ],
        ["dT", "prompt", "delayed"]
    )

    return ak.enforce_type(array, output_type)


def find_delayed_coincicence_candidates(hpge_data, mgc_data, max_delta_sec):
    # Merge the two dataframes on the "hit_table" column

    output_data = {
        "dT": [],
        "delayed": [],
        "prompt": []
    }

    for i in range(len(hpge_data)):
        for j in range(len(mgc_data)):

            dT = mgc_data[j]["coinc"]["id"]["timestamp"] - hpge_data[i]["geds"]["id"]["timestamp"]
            print(dT)
            if dT < 0:
                continue
            if dT > max_delta_sec:
                break

            output_data["prompt"].append(mgc_data[j].to_list())
            output_data["delayed"].append(hpge_data[i])
            output_data["dT"].append(dT)

    return ut.dict_to_lgdo(output_data)


def save_output(mdc_data,output):
    lh5.write(enforce_type(mdc_data), name="mdc", lh5_file=str(output))


def process_mu_delayed_coinc(input, output,default_ref_version="ref-v2.1.0", fallback_ref_version="ref-v2.0.0", max_delta_sec=100*53.7, metadata=None):

    mgc_data = load_mgc_data_of_run(input.mgc_files)
    hpge_data = generate_hpge_data( input.pht_files, default_ref_version=default_ref_version, fallback_ref_version=fallback_ref_version, metadata=metadata)
    mdc_data = find_delayed_coincicence_candidates(mgc_data,hpge_data, max_delta_sec=max_delta_sec)
    save_output(mdc_data,output)
