import os
import numpy as np
import pygama
from packaging.version import Version
from pygama import logging
from pygama.flow import DataLoader
from lgdo.lh5 import LH5Store, write
from lgdo.types import Table
from legendmeta import LegendMetadata
import pandas as pd
import json
import awkward as ak
import utils as ut

def generate_paths_of_different_tiers(input_path, fallback_defult="ref-v2.0.0"):
    """
    Generate paths for different tiers based on the input path.
    Args:
        input_path (str): Path to the input file.
    Returns:
        dict: Dictionary containing paths for pht, pet, dsp, and tcm.
    """
    paths = {}
    paths['pht'] = input_path
    paths['pet'] = input_path.replace("pht", "pet")
    paths["dsp"] = input_path.replace("pht", "dsp")
    if not os.path.exists(paths['dsp']):
        paths['dsp'] = input_path.replace("pht", "dsp").replace("ref-v2.1.0", fallback_defult)
    paths['tcm'] = input_path.replace("pht", "tcm")
    if not os.path.exists(paths['tcm']):
        paths['tcm'] = input_path.replace("pht", "tcm").replace("ref-v2.1.0", fallback_defult)
    return paths


def process_mu_hpge_coinc(input, output):
    """
    Process coincidences between muon channel and HPGe channel.
    Args:
        input (list): List of input files.
        output (list): List of output files.
    """

    paths = generate_paths_of_different_tiers(input)
    
    # Get tcm file
    store = LH5Store()
    pet_data_geds = store.read("/evt/geds/", paths["pet"])[0].view_as("ak")
    pet_data_coinc = store.read("/evt/coincident/", paths["pet"])[0].view_as("ak")

    timestamp = ut.extract_timestamp_raw(paths["pet"])
    chmap = ut.generate_channel_map(timestamp)
    data_streams_hpge = ut.select_datastreams(chmap, "HPGE")

    mask_muon_coinc = pet_data_coinc["muon_offline"] & ~pet_data_coinc["puls"]  # & pet_data["geds"]["is_good_hit"]

    tcm_id = store.read("/hardware_tcm_1/array_id", paths["tcm"])[0].view_as("np")
    tcm_idx = store.read("/hardware_tcm_1/array_idx", paths["tcm"])[0].view_as("np")

    # Load dsp and hit data
    columns = ["hit_table", "hit_idx", "evt_idx", "timestamp", "cuspEmax_ctc_cal", "is_good_hit","mu_diff","is_in_coincidence_with_mu","is_saturated"]
    acc_range = [-2000, 5000]

    selected_idx = ak.to_numpy(ak.flatten(pet_data_geds["hit_idx"][mask_muon_coinc]))

    output_data = {}
    for col in columns:
        output_data[col] = []
        
    for i in range(len(selected_idx)):
        hpge_idx = selected_idx[i]
        
        if not tcm_id[hpge_idx] in data_streams_hpge:
            continue
        
        muon_idx = np.where([(tcm_idx == tcm_idx[hpge_idx]) & (tcm_id == chmap.map("name")["MUON01"].daq.rawid)])[1][0]
            
        data_pht_hpge = store.read("ch{}/hit/".format(tcm_id[hpge_idx]),
                    paths["pht"],
                    idx=[tcm_idx[hpge_idx]])[0].view_as("ak")
        
        data_dsp_hpge = store.read("ch{}/dsp/".format(tcm_id[hpge_idx]),
                    paths["dsp"],
                    idx=[tcm_idx[hpge_idx]])[0].view_as("ak")

        data_dsp_muon = store.read("ch{}/dsp/".format(tcm_id[muon_idx]),
                    paths["dsp"],
                    idx=[tcm_idx[muon_idx]])[0].view_as("ak")

        if data_pht_hpge["cuspEmax_ctc_cal"] > 25:
            evt_idx = tcm_idx[hpge_idx]
            evt_id = np.where(pet_data_geds[evt_idx]["hit_idx"] == hpge_idx)[0][0]
            output_data["hit_table"].append(tcm_id[hpge_idx])
            output_data["hit_idx"].append(hpge_idx)
            output_data["evt_idx"].append(tcm_idx[hpge_idx])
            output_data["timestamp"].append(data_pht_hpge["timestamp"][0])
            output_data["cuspEmax_ctc_cal"].append(data_pht_hpge["cuspEmax_ctc_cal"][0])
            output_data["is_good_hit"].append(pet_data_geds[evt_idx]["is_good_hit"][evt_id])
            output_data["mu_diff"].append((data_dsp_hpge["tp_01"] - data_dsp_muon["tp_max"])[0])
            output_data["is_in_coincidence_with_mu"].append((acc_range[0] < output["mu_diff"][-1] < acc_range[1]))
            output_data["is_saturated"].append(data_pht_hpge["is_saturated"][0])

    # write to lh5 file
    output_lh5 = Table(col_dict=output_data)
    write(output_lh5, name="mgc", lh5_file=output)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process muon HPGe coincidences.")
    parser.add_argument("input", type=str, help="Input file path.")
    parser.add_argument("output", type=str, help="Output file path.")
    args = parser.parse_args()
    logging.info("Processing muon HPGe coincidences...")
    logging.info(f"Input: {args.input}")
    logging.info(f"Output: {args.output}")
    process_mu_hpge_coinc(args.input, args.output)
    logging.info("Processing completed.")