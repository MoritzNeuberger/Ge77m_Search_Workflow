
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
from lgdo.lh5 import LH5Store
import awkward as ak
import numpy as np
from tqdm import tqdm

"""


import glob
import utils as ut
import os
import pandas as pd
import pygama
from packaging.version import Version
from pygama import logging
from pygama.flow import DataLoader
from legendmeta import LegendMetadata
from lgdo.lh5 import LH5Store
import awkward as ak
import numpy as np
from tqdm import tqdm

def generate_channel_map(ts):
    lmeta = LegendMetadata()
    chmap = lmeta.channelmap(on=ts)
    return chmap

# Check pygama version
assert Version(pygama.__version__) >= Version("1.3.1")

# Set up logging level
logging.setup(level=logging.INFO)  # or level=logging.INFO, logging.DEBUG for more verbose output

# Set the PRODENV environment variable
os.environ['PRODENV'] = '/mnt/atlas01/projects/legend/data/prodenv/prod-blind/ref-v2.0.0'

# Define functions for data loading and processing

def get_muon_tagged_data_of_run(run_info,muon_tagged_folder):
    run_folder = ut.check_folder_path(os.path.abspath(ut.check_folder_path(muon_tagged_folder) + ut.generate_folder_path_from_run_info(run_info)))
    file_list = glob.glob(run_folder + "*.lh5")
    sorted_file_list = sorted(file_list, key=ut.extract_timestamp)
    output = []
    for path in sorted_file_list:
        df = pd.read_hdf(path)
        df["quality_is_bb_like"] = df["quality_is_bb_like"].astype(int)
        df["is_in_coincidence_with_mu"] = df["is_in_coincidence_with_mu"].astype(int)
        df["is_saturated"] = df["is_saturated"].astype(int)
        output.append(df)
    return pd.concat(output,ignore_index=True)
    
    
def get_hpge_data(run_info,config_file):

    
    # Get tcm file
    store = LH5Store()
    path_pet = "/mnt/atlas01/projects/legend/data/prodenv/prod-blind/ref-v2.0.0/generated/tier/pet/{}/{}/{}/l200-{}-{}-{}-{}-tier_pet.lh5".format(run_info["datatype"],run_info["period"],run_info["run"],run_info["period"],run_info["run"],run_info["datatype"],run_info["timestamp"])
 
    pet_data_geds= store.read("/evt/geds/",path_pet)[0].view_as("ak")
    pet_data_coinc= store.read("/evt/coincident/",path_pet)[0].view_as("ak")
    pet_data_trigger= store.read("/evt/trigger/",path_pet)[0].view_as("ak")

    chmap = generate_channel_map(run_info["timestamp"])

    mask = ~pet_data_coinc["muon_offline"] & ~pet_data_coinc["puls"] & (pet_data_geds["energy_sum"] > 50)

    path_tcm = "/mnt/atlas01/projects/legend/data/prodenv/prod-blind/ref-v2.0.0/generated/tier/tcm/{}/{}/{}/l200-{}-{}-{}-{}-tier_tcm.lh5".format(run_info["datatype"],run_info["period"],run_info["run"],run_info["period"],run_info["run"],run_info["datatype"],run_info["timestamp"])
    tcm_id = store.read("/hardware_tcm_1/array_id",path_tcm)[0].view_as("np")
    tcm_idx = store.read("/hardware_tcm_1/array_idx",path_tcm)[0].view_as("np")

    selected_idx = ak.to_numpy(ak.flatten(pet_data_geds["hit_idx"][mask]))
    selected_id = ak.to_numpy(ak.flatten(pet_data_geds["rawid"][mask]))
    
    columns = [
        "hit_table", 
        "hit_idx",
        "evt_idx",
        "cuspEmax_ctc_cal", 
        "timestamp", 
        "tp_01", 
        "tp_max", 
        "quality_is_bb_like",
        "psd_is_good",
        "psd_is_bb_like",
        "multiplicity",
        "spm_coinc"
    ]

    output = {}
    for col in columns:
        output[col] = []
    
    data_streams = ut.select_datastreams(chmap,"HPGE")
    for i in range(len(selected_idx)):
        #hpge_idx = selected_idx[i]
        
        hpge_idx = np.where((tcm_idx == selected_idx[i]) & (tcm_id == selected_id[i]))[0][0] #selected_idx[i]

        
        if tcm_id[hpge_idx] in data_streams:
        
            data_pht_hpge = store.read("ch{}/hit/".format(tcm_id[hpge_idx]),
                   "/mnt/atlas01/projects/legend/data/prodenv/prod-blind/ref-v2.0.0/generated/tier/pht/{}/{}/{}/l200-{}-{}-{}-{}-tier_pht.lh5".format(run_info["datatype"],run_info["period"],run_info["run"],run_info["period"],run_info["run"],run_info["datatype"],run_info["timestamp"]),
                   idx=[tcm_idx[hpge_idx]])[0].view_as("ak")
            
            if not "AoE_Double_Sided_Cut" in data_pht_hpge.fields: # and not "LQ_Cut" in data_pht_hpge.fields:
                continue
            
            data_dsp_hpge = store.read("ch{}/dsp/".format(tcm_id[hpge_idx]),
                   "/mnt/atlas01/projects/legend/data/prodenv/prod-blind/ref-v2.0.0/generated/tier/psp/{}/{}/{}/l200-{}-{}-{}-{}-tier_psp.lh5".format(run_info["datatype"],run_info["period"],run_info["run"],run_info["period"],run_info["run"],run_info["datatype"],run_info["timestamp"]),
                   idx=[tcm_idx[hpge_idx]])[0].view_as("ak")
            
            evt_idx = tcm_idx[hpge_idx]
            #print(pet_data_geds["hit_idx"][evt_idx],hpge_idx)
            evt_id = np.where(pet_data_geds["rawid"][evt_idx] == selected_id[i])[0][0]
            
            output["hit_table"].append(tcm_id[hpge_idx])
            output["hit_idx"].append(hpge_idx)
            output["evt_idx"].append(tcm_idx[hpge_idx])
            output["timestamp"].append(data_pht_hpge["timestamp"][0])
            output["cuspEmax_ctc_cal"].append(data_pht_hpge["cuspEmax_ctc_cal"][0])
            output["quality_is_bb_like"].append(pet_data_geds[evt_idx]["quality"]["is_bb_like"])
            output["psd_is_good"].append(pet_data_geds[evt_idx]["psd"]["is_good"][evt_id])
            output["tp_01"].append(data_dsp_hpge["tp_01"][0])
            output["tp_max"].append(data_dsp_hpge["tp_max"][0])
            output["psd_is_bb_like"].append(pet_data_geds[evt_idx]["psd"]["is_bb_like"][evt_id])
            output["multiplicity"].append(pet_data_geds[evt_idx]["multiplicity"])
            output["spm_coinc"].append(pet_data_coinc[evt_idx]["spms"])

    
    return pd.DataFrame(output)


def match_muon_HPGe(hpge_data, muon_tagged_data, max_delta_sec):
    # Merge the two dataframes on the "hit_table" column
    merged_data = pd.merge(hpge_data, muon_tagged_data, on="hit_table", suffixes=("_hpge", "_muon"))

    # Calculate the time difference
    merged_data["dt"] = merged_data["timestamp_hpge"] - merged_data["timestamp_muon"]

    # Filter based on time difference and maximum delta
    filtered_data = merged_data[(merged_data["dt"] > 0) & (merged_data["dt"] < max_delta_sec) & merged_data["is_in_coincidence_with_mu"]]
    
    output_data = pd.DataFrame(data={
        "hit_table": filtered_data["hit_table"].astype(int),
        "evt_idx_hpge": filtered_data["evt_idx_hpge"].astype(int),
        "timestamp_hpge": filtered_data["timestamp_hpge"].astype(float),
        "cuspEmax_ctc_cal_hpge": filtered_data["cuspEmax_ctc_cal_hpge"].astype(float),
        "quality_is_bb_like_hpge": filtered_data["quality_is_bb_like_hpge"].astype(int),
        "psd_is_good": filtered_data["psd_is_good"].astype(int),
        "psd_is_bb_like": filtered_data["psd_is_bb_like"].astype(int),
        "multiplicity": filtered_data["multiplicity"].astype(int),
        "spm_coinc": filtered_data["spm_coinc"].astype(int),
        "evt_idx_muon": filtered_data["evt_idx_muon"].astype(int),
        "timestamp_muon": filtered_data["timestamp_muon"].astype(float),
        "mu_diff": filtered_data["mu_diff"].astype(float),
        "dt": filtered_data["dt"].astype(float),
        "cuspEmax_ctc_cal_muon": filtered_data["cuspEmax_ctc_cal_muon"].astype(float),
        "quality_is_bb_like_muon": filtered_data["quality_is_bb_like_muon"].astype(int),
        "is_in_coincidence_with_mu": filtered_data["is_in_coincidence_with_mu"].astype(int),
        "is_saturated_muon": filtered_data["is_saturated"].astype(int)
    })
    
    output_data = output_data.sort_values("timestamp_hpge").reset_index(drop=True)

    return output_data

def save_output(run_info,output,output_folder):
    file_info = run_info
    if not output_folder[-1] == "/":
        output_folder += "/"
    output_folder = (output_folder 
                        + file_info["datatype"] 
                        + "/" + file_info["period"] 
                        + "/" + file_info["run"] 
                        + "/" )
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    output_file_name = (file_info["run"] 
                        + "-" + file_info["period"] 
                        + "-" + file_info["run"] 
                        + "-" + file_info["datatype"] 

                        + "-" + file_info["timestamp"] 
                        + "-" + "tier_mu_hpge_mdc.lh5")
    #ls = LH5Store()
    #ls.write_object(output,"muon_HPGe_mdc",output_file_name)
    #print(output)
    output.to_hdf((output_folder + output_file_name),"df")

muon_tagged_data = get_muon_tagged_data_of_run(run_info,args["muon_tagged_folder"])

hpge_data = get_hpge_data(run_info,args["config_file"])

merged_data = match_muon_HPGe(hpge_data, muon_tagged_data,args["max_time_gate"])

save_output(run_info,merged_data,args["output_folder"])

"""



def get_muon_tagged_data_of_run(run_info,muon_tagged_folder):
    run_folder = ut.check_folder_path(os.path.abspath(ut.check_folder_path(muon_tagged_folder) + ut.generate_folder_path_from_run_info(run_info)))
    file_list = glob.glob(run_folder + "*.lh5")
    sorted_file_list = sorted(file_list, key=ut.extract_timestamp)
    output = []
    for path in sorted_file_list:
        df = pd.read_hdf(path)
        df["quality_is_bb_like"] = df["quality_is_bb_like"].astype(int)
        df["is_in_coincidence_with_mu"] = df["is_in_coincidence_with_mu"].astype(int)
        df["is_saturated"] = df["is_saturated"].astype(int)
        output.append(df)
    return pd.concat(output,ignore_index=True)
    
    
def get_hpge_data(run_info,config_file):

    
    # Get tcm file
    store = LH5Store()
    path_pet = "/mnt/atlas01/projects/legend/data/prodenv/prod-blind/ref-v2.0.0/generated/tier/pet/{}/{}/{}/l200-{}-{}-{}-{}-tier_pet.lh5".format(run_info["datatype"],run_info["period"],run_info["run"],run_info["period"],run_info["run"],run_info["datatype"],run_info["timestamp"])
 
    pet_data_geds= store.read("/evt/geds/",path_pet)[0].view_as("ak")
    pet_data_coinc= store.read("/evt/coincident/",path_pet)[0].view_as("ak")
    pet_data_trigger= store.read("/evt/trigger/",path_pet)[0].view_as("ak")

    chmap = generate_channel_map(run_info["timestamp"])

    mask = ~pet_data_coinc["muon_offline"] & ~pet_data_coinc["puls"] & (pet_data_geds["energy_sum"] > 50)

    path_tcm = "/mnt/atlas01/projects/legend/data/prodenv/prod-blind/ref-v2.0.0/generated/tier/tcm/{}/{}/{}/l200-{}-{}-{}-{}-tier_tcm.lh5".format(run_info["datatype"],run_info["period"],run_info["run"],run_info["period"],run_info["run"],run_info["datatype"],run_info["timestamp"])
    tcm_id = store.read("/hardware_tcm_1/array_id",path_tcm)[0].view_as("np")
    tcm_idx = store.read("/hardware_tcm_1/array_idx",path_tcm)[0].view_as("np")

    selected_idx = ak.to_numpy(ak.flatten(pet_data_geds["hit_idx"][mask]))
    selected_id = ak.to_numpy(ak.flatten(pet_data_geds["rawid"][mask]))
    
    columns = [
        "hit_table", 
        "hit_idx",
        "evt_idx",
        "cuspEmax_ctc_cal", 
        "timestamp", 
        "tp_01", 
        "tp_max", 
        "quality_is_bb_like",
        "psd_is_good",
        "psd_is_bb_like",
        "multiplicity",
        "spm_coinc"
    ]

    output = {}
    for col in columns:
        output[col] = []
    
    data_streams = ut.select_datastreams(chmap,"HPGE")
    for i in range(len(selected_idx)):
        #hpge_idx = selected_idx[i]
        
        hpge_idx = np.where((tcm_idx == selected_idx[i]) & (tcm_id == selected_id[i]))[0][0] #selected_idx[i]

        
        if tcm_id[hpge_idx] in data_streams:
        
            data_pht_hpge = store.read("ch{}/hit/".format(tcm_id[hpge_idx]),
                   "/mnt/atlas01/projects/legend/data/prodenv/prod-blind/ref-v2.0.0/generated/tier/pht/{}/{}/{}/l200-{}-{}-{}-{}-tier_pht.lh5".format(run_info["datatype"],run_info["period"],run_info["run"],run_info["period"],run_info["run"],run_info["datatype"],run_info["timestamp"]),
                   idx=[tcm_idx[hpge_idx]])[0].view_as("ak")
            
            if not "AoE_Double_Sided_Cut" in data_pht_hpge.fields: # and not "LQ_Cut" in data_pht_hpge.fields:
                continue
            
            data_dsp_hpge = store.read("ch{}/dsp/".format(tcm_id[hpge_idx]),
                   "/mnt/atlas01/projects/legend/data/prodenv/prod-blind/ref-v2.0.0/generated/tier/psp/{}/{}/{}/l200-{}-{}-{}-{}-tier_psp.lh5".format(run_info["datatype"],run_info["period"],run_info["run"],run_info["period"],run_info["run"],run_info["datatype"],run_info["timestamp"]),
                   idx=[tcm_idx[hpge_idx]])[0].view_as("ak")
            
            evt_idx = tcm_idx[hpge_idx]
            #print(pet_data_geds["hit_idx"][evt_idx],hpge_idx)
            evt_id = np.where(pet_data_geds["rawid"][evt_idx] == selected_id[i])[0][0]
            
            output["hit_table"].append(tcm_id[hpge_idx])
            output["hit_idx"].append(hpge_idx)
            output["evt_idx"].append(tcm_idx[hpge_idx])
            output["timestamp"].append(data_pht_hpge["timestamp"][0])
            output["cuspEmax_ctc_cal"].append(data_pht_hpge["cuspEmax_ctc_cal"][0])
            output["quality_is_bb_like"].append(pet_data_geds[evt_idx]["quality"]["is_bb_like"])
            output["psd_is_good"].append(pet_data_geds[evt_idx]["psd"]["is_good"][evt_id])
            output["tp_01"].append(data_dsp_hpge["tp_01"][0])
            output["tp_max"].append(data_dsp_hpge["tp_max"][0])
            output["psd_is_bb_like"].append(pet_data_geds[evt_idx]["psd"]["is_bb_like"][evt_id])
            output["multiplicity"].append(pet_data_geds[evt_idx]["multiplicity"])
            output["spm_coinc"].append(pet_data_coinc[evt_idx]["spms"])

    
    return pd.DataFrame(output)


def match_muon_HPGe(hpge_data, muon_tagged_data, max_delta_sec):
    # Merge the two dataframes on the "hit_table" column
    merged_data = pd.merge(hpge_data, muon_tagged_data, on="hit_table", suffixes=("_hpge", "_muon"))

    # Calculate the time difference
    merged_data["dt"] = merged_data["timestamp_hpge"] - merged_data["timestamp_muon"]

    # Filter based on time difference and maximum delta
    filtered_data = merged_data[(merged_data["dt"] > 0) & (merged_data["dt"] < max_delta_sec) & merged_data["is_in_coincidence_with_mu"]]
    
    output_data = pd.DataFrame(data={
        "hit_table": filtered_data["hit_table"].astype(int),
        "evt_idx_hpge": filtered_data["evt_idx_hpge"].astype(int),
        "timestamp_hpge": filtered_data["timestamp_hpge"].astype(float),
        "cuspEmax_ctc_cal_hpge": filtered_data["cuspEmax_ctc_cal_hpge"].astype(float),
        "quality_is_bb_like_hpge": filtered_data["quality_is_bb_like_hpge"].astype(int),
        "psd_is_good": filtered_data["psd_is_good"].astype(int),
        "psd_is_bb_like": filtered_data["psd_is_bb_like"].astype(int),
        "multiplicity": filtered_data["multiplicity"].astype(int),
        "spm_coinc": filtered_data["spm_coinc"].astype(int),
        "evt_idx_muon": filtered_data["evt_idx_muon"].astype(int),
        "timestamp_muon": filtered_data["timestamp_muon"].astype(float),
        "mu_diff": filtered_data["mu_diff"].astype(float),
        "dt": filtered_data["dt"].astype(float),
        "cuspEmax_ctc_cal_muon": filtered_data["cuspEmax_ctc_cal_muon"].astype(float),
        "quality_is_bb_like_muon": filtered_data["quality_is_bb_like_muon"].astype(int),
        "is_in_coincidence_with_mu": filtered_data["is_in_coincidence_with_mu"].astype(int),
        "is_saturated_muon": filtered_data["is_saturated"].astype(int)
    })
    
    output_data = output_data.sort_values("timestamp_hpge").reset_index(drop=True)

    return output_data

def save_output(run_info,output,output_folder):
    file_info = run_info
    if not output_folder[-1] == "/":
        output_folder += "/"
    output_folder = (output_folder 
                        + file_info["datatype"] 
                        + "/" + file_info["period"] 
                        + "/" + file_info["run"] 
                        + "/" )
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    output_file_name = (file_info["run"] 
                        + "-" + file_info["period"] 
                        + "-" + file_info["run"] 
                        + "-" + file_info["datatype"] 

                        + "-" + file_info["timestamp"] 
                        + "-" + "tier_mu_hpge_mdc.lh5")
    #ls = LH5Store()
    #ls.write_object(output,"muon_HPGe_mdc",output_file_name)
    #print(output)
    output.to_hdf((output_folder + output_file_name),"df")


def process_mu_delayed_coinc(input, output):

    print("Input: ", input, type(input), input.mgc_files, input.pht_files)
    print("Output: ", output)

#    
#    muon_tagged_data = get_muon_tagged_data_of_run(run_info,args["muon_tagged_folder"])
#
#    hpge_data = get_hpge_data(run_info,args["config_file"])
#
#    merged_data = match_muon_HPGe(hpge_data, muon_tagged_data,args["max_time_gate"])
#
#    save_output(run_info,merged_data,args["output_folder"])
