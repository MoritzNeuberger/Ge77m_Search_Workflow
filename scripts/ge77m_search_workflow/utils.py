def check_folder_path(path):
    if not path[-1] == "/":
        path += "/"
    return path

def generate_folder_path_from_run_info(run_info):
    return "{}/{}/{}/".format(run_info["datatype"],run_info["period"],run_info["run"])

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

def generate_channel_map(dl):
    file_list = dl.get_file_list()[["experiment", "period", "run", "timestamp", "datatype"]]
    first_key = file_list.iloc[0].timestamp
    lmeta = LegendMetadata()
    chmap = lmeta.channelmap(on=first_key)
    return chmap

def get_HPGe_channels(chmap):
    channels = chmap.map("system", unique=False)["geds"]
    channels = channels.map("analysis.usability", unique=False)['on']
    return channels

def select_datastreams(chmap,stream):
    if stream == "HPGE":
        datastreams = get_HPGe_channels(chmap).map("daq.rawid").keys()
    else:
        datastreams = [chmap.map("name")[stream].daq.rawid]
    return datastreams


def get_entry_list(dl, chmap, stream, cuts):
    dl.set_datastreams(select_datastreams(chmap,stream), "ch")
    dl.set_cuts(cuts)
    dl.set_output(columns=["trapEmax_ctc_cal","trapTmax", "timestamp","tp_01", "tp_max", "wf_max", "waveform"], fmt="pd.DataFrame")
    el = dl.build_entry_list(tcm_level="tcm")
    add_evt_idx_abs(el)
    return el

def get_data(dl, chmap, stream, entry_list, columns):
    dl.set_datastreams(select_datastreams(chmap,stream), "ch")
    dl.set_output(columns=columns, fmt="pd.DataFrame")
    return dl.load(entry_list)

def get_data_new(dl, chmap, stream, entry_list, columns):
    dl.set_datastreams(select_datastreams(chmap,stream), "ch")
    dl.set_output(columns=columns, fmt="pd.DataFrame")
    return dl.load(entry_list)


def add_evt_idx_abs(el):
    el["evt_idx_abs"] = el.file * 10000 + el.evt_idx
