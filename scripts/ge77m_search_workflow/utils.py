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
    Find the file that actually contains the target timestamp by checking timestamp ranges.
    
    Args:
        timestamp (float): The timestamp to compare against.
        timestamp_path_dict (dict): A dictionary mapping timestamps to file paths.
        
    Returns:
        str: The path of the file that contains the timestamp, or None if not found.
    """
    import lgdo.lh5 as lh5
    
    # First, try the original approach - find files with larger timestamps
    candidate_files = []
    for ts in sorted(timestamp_path_dict.keys()):
        if ts > timestamp:
            candidate_files.append(timestamp_path_dict[ts])
    
    # If no files found with larger timestamps, also check smaller ones
    if not candidate_files:
        for ts in sorted(timestamp_path_dict.keys(), reverse=True):
            if ts <= timestamp:
                candidate_files.append(timestamp_path_dict[ts])
    
    # Now check each candidate file to see if it actually contains the timestamp
    for file_path in candidate_files:
        try:
            # Get list of available groups in the file
            groups = lh5.ls(file_path)
            
            # Look for channel groups (ch + number pattern)
            channels = [group for group in groups if group.startswith('ch') and group[2:].isdigit()]
            
            if not channels:
                continue
            
            # Try the first available channel
            test_channel = channels[0]
            timestamps_path = f"{test_channel}/hit/timestamp"  # Changed from raw to hit
            
            # Check if timestamp field exists
            try:
                # Read first few timestamps
                first_timestamps = lh5.read(timestamps_path, file_path, n_rows=5)
                
                # Get the total length by reading the shape
                full_timestamps = lh5.read(timestamps_path, file_path)
                if hasattr(full_timestamps, 'nda'):
                    total_length = len(full_timestamps.nda)
                elif hasattr(full_timestamps, 'value'):
                    total_length = len(full_timestamps.value)
                else:
                    total_length = len(full_timestamps)
                
                # Read last few timestamps
                start_row = max(0, total_length - 5)
                n_rows = min(5, total_length - start_row)
                last_timestamps = lh5.read(timestamps_path, file_path, start_row=start_row, n_rows=n_rows)
                
                # Extract timestamp values from LGDO objects
                if hasattr(first_timestamps, 'nda'):
                    first_ts = first_timestamps.nda[0]
                elif hasattr(first_timestamps, 'value'):
                    first_ts = first_timestamps.value[0] if hasattr(first_timestamps.value, '__len__') else first_timestamps.value
                else:
                    first_ts = first_timestamps[0]
                
                if hasattr(last_timestamps, 'nda'):
                    last_ts = last_timestamps.nda[-1]
                elif hasattr(last_timestamps, 'value'):
                    last_ts = last_timestamps.value[-1] if hasattr(last_timestamps.value, '__len__') else last_timestamps.value
                else:
                    last_ts = last_timestamps[-1]
                
                # Check if our target timestamp falls within this range
                if first_ts <= timestamp <= last_ts:
                    return file_path
                    
            except Exception as inner_e:
                print(f"Warning: Could not read timestamps from {timestamps_path} in {file_path}: {inner_e}")
                continue
                        
        except Exception as e:
            # If we can't read this file, continue to the next one
            print(f"Warning: Could not check timestamp range for {file_path}: {e}")
            continue
    
    # If no file contains the timestamp, return None
    return None

from pygama.flow import DataLoader
from legendmeta import LegendMetadata
import json

def generate_data_loader(run_info,config_file):
    with open(config_file,"r") as f:
        filedb = json.load(f)["setups"]["l200"]["db"]
    dl = DataLoader(config_file + "[setups/l200/dataloader]", filedb=filedb)
    file_query = "period == '{}' and datatype == '{}' and run == '{}' and timestamp == '{}'".format(run_info['period'],run_info['datatype'],run_info['run'],run_info['timestamp'])
    dl.set_files(file_query)
    return dl

def generate_channel_map(ts,metadata=None):
    if metadata is not None:
        if isinstance(metadata, str):
            lmeta = LegendMetadata(metadata)
        elif isinstance(metadata, LegendMetadata):
            lmeta = metadata
    else:
        lmeta = LegendMetadata()
    chmap = lmeta.channelmap(on=ts)
    return chmap

def get_HPGe_channels(chmap,usability_condition='on',psd_ms_type=None):
    channels = chmap.map("system", unique=False)["geds"]
    channels = channels.map("analysis.usability", unique=False)[usability_condition]
    if isinstance(psd_ms_type, str):
        if psd_ms_type == "bb_like":
            channels = channels.map("daq.rawid", unique=False)  # return all channels
            output = []
            for ch in channels:
                is_bb_like = channels[ch]["analysis"]["psd"]["is_bb_like"].split(" & ")
                is_valid = True
                for psd_type in is_bb_like:
                    if psd_type == "missing" or channels[ch]["analysis"]["psd"]["status"][psd_type] != "valid":
                        is_valid = False
                        break
                if is_valid:
                    output.append(channels[ch]["daq"]["rawid"])
            return output
        else:
            try:
                channels = channels.map("analysis.psd.status."+psd_ms_type, unique=False)
                if "valid" in channels:
                    channels = channels['valid']
                else:
                    return None
            except KeyError:
                print(channels)
                raise
            channels = channels.map("daq.rawid", unique=False) 
            return channels
    elif isinstance(psd_ms_type, list):
        tmp_channels =  []
        for psd_type in psd_ms_type:
            try:
                tmp = channels.map("analysis.psd.status."+psd_type, unique=False)
                if "valid" in tmp:
                    tmp = tmp["valid"]
                    tmp = tmp.map("daq.rawid", unique=False)
                    tmp_channels.append(tmp)
            except KeyError:
                print(channels)
                raise
        return tmp_channels
    else:
        channels = channels.map("daq.rawid", unique=False)
        return channels

def select_datastreams(chmap,stream,psd_ms_type=None):

    psd_ms_type_options = ["low_aoe", "coax_rt", "bb_like"]

    if psd_ms_type is not None and isinstance(psd_ms_type,str) and psd_ms_type not in psd_ms_type_options:
        if psd_ms_type == "any":
            psd_ms_type = None
        else:
            raise ValueError(f"Invalid psd_ms_type: {psd_ms_type}. Must be one of {psd_ms_type_options} or None.")

    if psd_ms_type is not None and isinstance(psd_ms_type, list):
        for psd_type in psd_ms_type:
            if psd_type not in psd_ms_type_options:
                raise ValueError(f"Invalid psd_ms_type: {psd_type}. Must be one of {psd_ms_type_options} or None.")
    if stream == "HPGE":
        HPGe_channels = get_HPGe_channels(chmap, psd_ms_type=psd_ms_type)
        if HPGe_channels is None:
            return []
        if isinstance(HPGe_channels, list):
            datastreams = []
            if isinstance(HPGe_channels[0], dict):
                for ch in HPGe_channels:
                    datastreams.extend(ch.keys())
            else:
                for ch in HPGe_channels:
                    datastreams.append(ch)
        else:
            datastreams = HPGe_channels.keys()
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
            ak.types.RecordType(
                [
                    ak.types.NumpyType("float64"),  # trapEmax_ctc_cal
                    ak.types.NumpyType("float64"),  # cuspEmax_ctc_cal
                    ak.types.NumpyType("float64"),  # zacEmax_ctc_cal
                ],
                ["trapEmax_ctc_cal", "cuspEmax_ctc_cal", "zacEmax_ctc_cal"]
            ),  # other_energy_estimators
            ak.types.NumpyType("float64"),  # tp_max
            ak.types.NumpyType("float64"),  # wf_max
            ak.types.NumpyType("float64"),  # tailEmax
        ],
        ["energy", "tp_01", "quality", "id", "other_energy_estimators", "tp_max", "wf_max", "tailEmax"]
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
            ak.types.RecordType(
                [
                    ak.types.NumpyType("int64"),   # evt_idx
                    ak.types.NumpyType("float64"), # timestamp
                ],
                ["evt_idx", "timestamp"]
            ),  # id
        ],
        ["mu_diff", "id"]
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

    return ak.enforce_type(array, delayed_type)


import re
def expand_range_token(token: list[str]) -> list[str]:
    """
    Given a token like 'r000..r005', returns ['r000','r001',...,'r005'].
    """

    output = []
    for t in token:
        # 1) Parse out the two endpoints
        left, right = t.split("..")
        # 2) Split each endpoint into its alphabetic prefix and numeric suffix
        m1 = re.match(r"^([A-Za-z_]*)(\d+)$", left)
        m2 = re.match(r"^([A-Za-z_]*)(\d+)$", right)
        if not (m1 and m2):
            raise ValueError(f"Token {t!r} is not in expected format prefix+num..prefix+num")
        prefix1, num1 = m1.groups()
        prefix2, num2 = m2.groups()
        # 3) They should share the same prefix and same width
        if prefix1 != prefix2:
            raise ValueError("Mismatched prefixes: "
                            f"{prefix1!r} vs {prefix2!r}")
        width = max(len(num1), len(num2))
        start, end = int(num1), int(num2)
        step = 1 if start <= end else -1
        output.append([
            f"{prefix1}{i:0{width}d}"
            for i in range(start, end + step, step)
        ])
    # 4) Generate the list
    return np.concatenate(output).tolist()