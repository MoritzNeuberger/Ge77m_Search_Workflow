def check_folder_path(path):
    if not path[-1] == "/":
        path += "/"
    return path

def generate_folder_path_from_run_info(run_info):
    return "{}/{}/{}/".format(run_info["datatype"],run_info["period"],run_info["run"])

from lgdo import types
import os

def dict_to_lgdo(in_dict):
    """
    A helper function that turns a nested dict containing equal sized arrays to a lgdo.type.Table, 
    with each sub dict also being turn into a lgdo.type.Table and each array into a lgdo.type.Array.
    """
    tables = {}
    for key, value in in_dict.items():
        if isinstance(value, dict):
            tables[key] = dict_to_lgdo(value)
        else:
            print(value)
            tables[key] = types.Array(value)
    return types.Table(tables)

def generate_paths_of_different_tiers_from_pht(input_path, default_ref_version, fallback_ref_version, input_raw=None):
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
    paths["psp"] = input_path.replace("pht", "psp")
    if not os.path.exists(paths['psp']):
        paths['psp'] = input_path.replace("pht", "psp").replace(default_ref_version, fallback_ref_version)
    paths['tcm'] = input_path.replace("pht", "tcm")
    if not os.path.exists(paths['tcm']):
        paths['tcm'] = input_path.replace("pht", "tcm").replace(default_ref_version, fallback_ref_version)
    if input_raw is not None:
        paths['raw'] = (input_raw + "/" + "/".join(input_path.rsplit("/",4)[1:])).replace("pht", "raw")
    return paths

#/dvs_ro/cfs/cdirs/m2676/users/pertoldi/legend-prodenv/prod-blind/ref-v2.1.0/generated/tier/pht/phy/p03/r001/l200-p03-r001-phy-20230320T025504Z-tier_pht.lh5


from datetime import datetime 
import re
# Define a custom key function to extract and convert the timestamp
def extract_timestamp(filename):
    # Use regular expression to find the timestamp in the filename
    match = re.search(r'\d{8}T\d{6}Z', filename)
    if match:
        timestamp_str = match.group()
        timestamp = datetime.strptime(timestamp_str, '%Y%m%dT%H%M%SZ')
        return timestamp
    else:
        # Handle cases where the timestamp is not found
        return datetime.min
    
def extract_timestamp_raw(filename):
    # Use regular expression to find the timestamp in the filename
    match = re.search(r'\d{8}T\d{6}Z', filename)
    if match:
        timestamp_str = match.group()
        return timestamp_str
    else:
        # Handle cases where the timestamp is not found
        return None
    
import glob
import numpy as np

def generate_timestamp_to_path_dict(input_path):
    input_paths = glob.glob(input_path + "/*/*/*.lh5")
    timestamps = []
    paths = []
    for path in input_paths:
        timestamps.append(extract_timestamp(path).timestamp())
        paths.append(path)
    arg_sort_mask = np.argsort(timestamps)
    sorted_timestamps = np.array(timestamps)[arg_sort_mask]
    sorted_paths = np.array(paths)[arg_sort_mask]
    return dict(zip(sorted_timestamps, sorted_paths))

def find_file_from_closest_larger_timestamp(timestamp, timestamp_path_dict):
    """
    Find the file with the closest larger timestamp in the timestamp_path_dict.
    
    Args:
        timestamp (float): The timestamp to compare against.
        timestamp_path_dict (dict): A dictionary mapping timestamps to file paths.
        
    Returns:
        str: The path of the file with the closest larger timestamp, or None if not found.
    """
    for ts in sorted(timestamp_path_dict.keys()):
        if ts > timestamp:
            return timestamp_path_dict[ts]
    return None

from pygama.flow import DataLoader
from legendmeta import LegendMetadata

def generate_data_loader(run_info,config_file):
    with open(config_file,"r") as f:
        filedb = json.load(f)["setups"]["l200"]["db"]
    dl = DataLoader(config_file + "[setups/l200/dataloader]", filedb=filedb)
    file_query = "period == '{}' and datatype == '{}' and run == '{}' and timestamp == '{}'".format(run_info['period'],run_info['datatype'],run_info['run'],run_info['timestamp'])
    dl.set_files(file_query)
    return dl

def generate_channel_map(ts,metadata=None):
    if metadata is not None:
        lmeta = LegendMetadata(metadata)
    else:
        lmeta = LegendMetadata()
    chmap = lmeta.channelmap(on=ts)
    return chmap

def get_HPGe_channels(chmap,usability_condition='on',psd_type="low_aoe"):
    channels = chmap.map("system", unique=False)["geds"]
    channels = channels.map("analysis.usability", unique=False)[usability_condition]
    #channels = channels.map("analysis.psd.status."+psd_type, unique=False)['valid']
    channels = channels.map("daq.rawid", unique=False)
    return channels

def select_datastreams(chmap,stream):
    if stream == "HPGE":
        datastreams = get_HPGe_channels(chmap).map("daq.rawid").keys()
    else:
        datastreams = [chmap.map("name")[stream].daq.rawid]
    return datastreams


import awkward as ak


def enforce_type_prompt(array):
    
    geds_type = ak.types.RecordType(
        [
            ak.types.NumpyType("float64"),  # energy
            ak.types.NumpyType("float64"),  # tp_01
            ak.types.RecordType(
                [
                    ak.types.NumpyType("bool"),  # is_bb_like
                    ak.types.NumpyType("bool"),  # is_good_channel
                    ak.types.NumpyType("bool"),  # is_saturated
                    ak.types.NumpyType("bool"),  # channel_is_positive_polarity
                ],
                ["is_bb_like", "is_good_channel", "is_saturated", "channel_is_positive_polarity"]
            ),  # quality
            ak.types.RecordType(
                [
                    ak.types.NumpyType("int64"),  # hit_table
                    ak.types.NumpyType("int64"),  # hit_idx
                ],
                ["hit_table", "hit_idx"]
            ),  # id
        ],
        ["energy", "tp_01", "quality", "id"]
    )

    mu_type = ak.types.RecordType(
        [
            ak.types.NumpyType("float64"),  # tp_max
            ak.types.RecordType(
                [
                    ak.types.NumpyType("int64"),  # hit_table
                    ak.types.NumpyType("int64"),  # hit_idx
                ],
                ["hit_table", "hit_idx"]
            ),  # id
            ak.types.RecordType(
                [
                    ak.types.NumpyType("bool"),  # evt_idx
                    ak.types.NumpyType("bool"),  # timestamp
                ],
                ["muon","muon_offline"]
            )
        ],
        ["tp_max", "id", "coinc_flags"]
    )

    coinc_type = ak.types.RecordType(
        [
            ak.types.NumpyType("float64"),  # mu_diff
            ak.types.NumpyType("bool"),     # is_in_coincidence_with_mu
            ak.types.RecordType(
                [
                    ak.types.NumpyType("int64"),   # evt_idx
                    ak.types.NumpyType("float64"), # timestamp
                ],
                ["evt_idx", "timestamp"]
            ),  # id
        ],
        ["mu_diff", "is_in_coincidence_with_mu", "id"]
    )

    output_type = ak.types.RecordType(
        [
            geds_type,
            mu_type,
            coinc_type,
        ],
        ["geds", "mu", "coinc"]
    )

    return ak.enforce_type(array, output_type)



def enforce_type_delayed(array):
    """
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
    """

    
    delayed_type = ak.types.RecordType([

        ak.types.RecordType(
            [
                ak.types.NumpyType("float64"),  # energy
                ak.types.RecordType(
                    [
                        ak.types.NumpyType("bool"),    # quality_is_bb_like
                        ak.types.NumpyType("bool"),    # psd_is_good
                        ak.types.NumpyType("bool"),    # psd_is_bb_like
                    ],
                    ["quality_is_bb_like", "psd_is_good", "psd_is_bb_like"]
                ),  # quality
                ak.types.RecordType(
                    [
                        ak.types.NumpyType("int64"),   # hit_table
                        ak.types.NumpyType("int64"),   # hit_idx
                        ak.types.NumpyType("float64"), # timestamp
                    ],
                    ["hit_table", "hit_idx", "timestamp"]
                ),  # id
            ],
            ["energy", "quality", "id"]
        ),
        ak.types.RecordType(
            [
                ak.types.NumpyType("int64"),  # multiplicity
                ak.types.NumpyType("bool"),   # spm_coinc
            ],
            ["multiplicity", "spm_coinc"]
        )  # coinc
        ],
        ["geds", "coinc"]
    )

    print(array, array.fields)

    return ak.enforce_type(array, delayed_type)
