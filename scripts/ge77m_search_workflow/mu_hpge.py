import os
#import numpy as np
#import pygama
#from packaging.version import Version
#from pygama import logging
#from pygama.flow import DataLoader
#from lgdo.lh5 import LH5Store
#from legendmeta import LegendMetadata
#import pandas as pd
#import json
#import awkward as ak

"""

def generate_channel_map(ts):
    lmeta = LegendMetadata()
    chmap = lmeta.channelmap(on=ts)
    return chmap


def get_HPGe_channels(chmap):
    channels = chmap.map("system", unique=False)["geds"]
    channels = channels.map("analysis.usability", unique=False)["on"]
    return channels


def select_datastreams(chmap,stream):
    if stream == "HPGE":
        datastreams = get_HPGe_channels(chmap).map("daq.rawid").keys()
    else:
        datastreams = [chmap.map("name")[stream].daq.rawid]
    return datastreams


def get_data(dl, chmap, stream, entry_list, columns):
    if len(entry_list) == 0:
        return pd.DataFrame(columns=[*columns,"evt_idx_abs",'hit_table', 'hit_idx', 'evt_idx', 'file'])
    dl.set_datastreams(select_datastreams(chmap,stream), "ch")
    dl.set_output(columns=columns, fmt="pd.DataFrame")
    return dl.load(entry_list)

# Define the is_in_coincidence_with_mu function
def is_in_coincidence_with_mu(hpge_data, mu_data, acc_range):
    hpge_data_np = hpge_data[['tp_01', 'evt_idx_abs']].to_numpy()
    mu_data_np = mu_data[['evt_idx_abs', 'tp_max']].to_numpy()

    def nb_func(hpge_data_np,mu_data_np,acc_range):
        output_bool = []
        output_diff = []
        for i in range(len(hpge_data_np)):
            sig_tp_01 = hpge_data_np[i, 0]
            index = hpge_data_np[i, 1]
            mu_index = np.where(mu_data_np[:, 0] == index)[0]
            if len(mu_index) > 0:
                mu_tp_max = mu_data_np[mu_index[0], 1]
                diff = sig_tp_01 - mu_tp_max 
                if acc_range[0] < diff < acc_range[1]:
                    output_bool.append(True)
                    output_diff.append(diff)
                else:
                    output_bool.append(False)
                    output_diff.append(np.nan)
            else:
                output_bool.append(False)
                output_diff.append(np.nan)
        return output_diff, output_bool
    return nb_func(hpge_data_np,mu_data_np,acc_range)


def check_folder_path(path):
    if not path[-1] == "/":
        path += "/"
    return path

def save_output(run_info,output,output_folder):
    if not output_folder[-1] == "/":
        output_folder += "/"
    output_folder = (output_folder 
                        + run_info["datatype"] 
                        + "/" + run_info["period"] 
                        + "/" + run_info["run"] 
                        + "/" )
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    output_file_name = ("l200" 
                        + "-" + run_info["period"] 
                        + "-" + run_info["run"] 
                        + "-" + run_info["datatype"] 
                        + "-" + run_info["timestamp"] 
                        + "-" + "tier_mu_coinc.lh5")
    output.to_hdf((output_folder + output_file_name),"df")




# Parse command-line arguments
parser = argparse.ArgumentParser(description='Generate muon tagged HPGe files')
parser.add_argument('-p', '--period', help='period', required=True)
parser.add_argument('-d', '--datatype', help='datatype', required=True)
parser.add_argument('-r', '--run', help='run', required=True)
parser.add_argument('-t', '--timestamp', help='timestamp', required=True)
parser.add_argument('-o', '--output_folder', help='output_folder', required=True)
parser.add_argument('-c', '--config_file', help='config_file', required=True)
args = vars(parser.parse_args())

run_info = {
    "period": args["period"],
    "datatype": args["datatype"],
    "run": args["run"],
    "timestamp": args["timestamp"]
}

# Get tcm filestore =LH5Store()
path_pet = "/global/cfs/projectdirs/m2676/data/lngs/l200/public/prodenv-new/prod-blind/ref-v1.1.0/generated/tier/pet/{}/{}/{}/l200-{}-{}-{}-{}-tier_pet.lh5".format(run_info["datatype"],run_info["period"],run_info["run"],run_info["period"],run_info["run"],run_info["datatype"],run_info["timestamp"])
pet_data_geds= store.read("/evt/geds/",path_pet)[0].view_as("ak")
pet_data_coinc= store.read("/evt/coincident/",path_pet)[0].view_as("ak")
pet_data_trigger= store.read("/evt/trigger/",path_pet)[0].view_as("ak")


chmap = generate_channel_map(run_info["timestamp"])
data_streams_hpge = ut.select_datastreams(chmap,"HPGE")

mask_muon_coinc = pet_data_coinc["muon_offline"] & ~pet_data_coinc["puls"]# & pet_data["geds"]["is_good_hit"]

path_tcm = "/global/cfs/projectdirs/m2676/data/lngs/l200/public/prodenv-new/prod-blind/ref-v1.0.0/generated/tier/tcm/{}/{}/{}/l200-{}-{}-{}-{}-tier_tcm.lh5".format(run_info["datatype"],run_info["period"],run_info["run"],run_info["period"],run_info["run"],run_info["datatype"],run_info["timestamp"])
tcm_id = store.read("/hardware_tcm_1/array_id",path_tcm)[0].view_as("np")
tcm_idx = store.read("/hardware_tcm_1/array_idx",path_tcm)[0].view_as("np")

# Load dsp and hit data
columns = ["hit_table", "hit_idx", "evt_idx", "timestamp", "cuspEmax_ctc_cal", "is_good_hit","mu_diff","is_in_coincidence_with_mu","is_saturated"]
acc_range = [-2000, 5000]

selected_idx = ak.to_numpy(ak.flatten(pet_data_geds["hit_idx"][mask_muon_coinc]))

output = {}
for col in columns:
    output[col] = []
    
print(data_streams_hpge)
    
for i in range(len(selected_idx)):
    hpge_idx = selected_idx[i]
    
    print(tcm_id[hpge_idx],tcm_id[hpge_idx] in data_streams_hpge)
    if not tcm_id[hpge_idx] in data_streams_hpge:
        continue
    
    muon_idx = np.where([(tcm_idx == tcm_idx[hpge_idx]) & (tcm_id == chmap.map("name")["MUON01"].daq.rawid)])[1][0]
        
    data_pht_hpge = store.read("ch{}/hit/".format(tcm_id[hpge_idx]),
                   "/global/cfs/projectdirs/m2676/data/lngs/l200/public/prodenv-new/prod-blind/ref-v1.1.0/generated/tier/pht/{}/{}/{}/l200-{}-{}-{}-{}-tier_pht.lh5".format(run_info["datatype"],run_info["period"],run_info["run"],run_info["period"],run_info["run"],run_info["datatype"],run_info["timestamp"]),
                   idx=[tcm_idx[hpge_idx]])[0].view_as("ak")
    
    data_dsp_hpge = store.read("ch{}/dsp/".format(tcm_id[hpge_idx]),
                   "/global/cfs/projectdirs/m2676/data/lngs/l200/public/prodenv-new/prod-blind/ref-v1.0.0/generated/tier/dsp/{}/{}/{}/l200-{}-{}-{}-{}-tier_dsp.lh5".format(run_info["datatype"],run_info["period"],run_info["run"],run_info["period"],run_info["run"],run_info["datatype"],run_info["timestamp"]),
                   idx=[tcm_idx[hpge_idx]])[0].view_as("ak")
    
    data_dsp_muon = store.read("ch{}/dsp/".format(tcm_id[muon_idx]),
                   "/global/cfs/projectdirs/m2676/data/lngs/l200/public/prodenv-new/prod-blind/ref-v1.0.0/generated/tier/dsp/{}/{}/{}/l200-{}-{}-{}-{}-tier_dsp.lh5".format(run_info["datatype"],run_info["period"],run_info["run"],run_info["period"],run_info["run"],run_info["datatype"],run_info["timestamp"]),
                   idx=[tcm_idx[muon_idx]])[0].view_as("ak")

    if data_pht_hpge["cuspEmax_ctc_cal"] > 10:
        evt_idx = tcm_idx[hpge_idx]
        evt_id = np.where(pet_data_geds[evt_idx]["hit_idx"] == hpge_idx)[0][0]
        output["hit_table"].append(tcm_id[hpge_idx])
        output["hit_idx"].append(hpge_idx)
        output["evt_idx"].append(tcm_idx[hpge_idx])
        output["timestamp"].append(data_pht_hpge["timestamp"][0])
        output["cuspEmax_ctc_cal"].append(data_pht_hpge["cuspEmax_ctc_cal"][0])
        output["is_good_hit"].append(pet_data_geds[evt_idx]["is_good_hit"][evt_id])
        output["mu_diff"].append((data_dsp_hpge["tp_01"] - data_dsp_muon["tp_max"])[0])
        output["is_in_coincidence_with_mu"].append((acc_range[0] < output["mu_diff"][-1] < acc_range[1]))
        output["is_saturated"].append(data_pht_hpge["is_saturated"][0])

output = pd.DataFrame(output)
#output = hpge_data[["hit_table", "hit_idx", "evt_idx", "file", "evt_idx_abs", "timestamp", "trapEmax_ctc_cal", "is_physical", "mu_diff", "is_in_coincidence_with_mu"]]

# Save the output
save_output(run_info, output, args["output_folder"])

"""


def process_mu_hpge_coinc(input, output):
    """
    Process coincidences between muon channel and HPGe channel.
    Args:
        input (list): List of input files.
        output (list): List of output files.
    """










    print("Processing mu_hpge_coinc...")
    print("input:", len(input))
    print("output:", len(output))    
    os.system(f"touch {output}")