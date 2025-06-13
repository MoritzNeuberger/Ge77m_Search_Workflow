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
            tables[key] = types.Array(value)
    return types.Table(tables)

def generate_paths_of_different_tiers_from_pht(input_path, default_ref_version, fallback_ref_version, raw_path=None):
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
    if raw_path is not None:
        paths['raw'] = raw_path + "/".join(input_path.rsplit("/",4)[1:])
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
    channels = channels.map("analysis.psd.status."+psd_type, unique=False)['valid']
    channels = channels.map("daq.rawid", unique=False)
    return channels

def select_datastreams(chmap,stream):
    if stream == "HPGE":
        datastreams = get_HPGe_channels(chmap).map("daq.rawid").keys()
    else:
        datastreams = [chmap.map("name")[stream].daq.rawid]
    return datastreams
