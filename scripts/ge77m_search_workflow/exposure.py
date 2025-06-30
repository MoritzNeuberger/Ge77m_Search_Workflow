from legendmeta import LegendMetadata
from ge77m_search_workflow import utils as ut
import yaml
from tqdm import tqdm

def calculate_exposure(output_path, metadata_path, psd_ms_type, ignore_periods=None, physics_dataset="ovbb"):

    exposures = {}

    lmeta = LegendMetadata(metadata_path)
    p_dataset = lmeta.datasets.runlists[physics_dataset]["phy"]

    for period in tqdm(lmeta.datasets.runinfo.keys()):

        if not ignore_periods is None and period in ignore_periods:
            continue

        if period not in exposures:
            exposures[period] = {}

        if not period in p_dataset:
            continue
        
        runs_in_physics_period = ut.expand_range_token(p_dataset[period])

        for run in lmeta.datasets.runinfo[period].keys():

            if run not in runs_in_physics_period:
                continue

            if "phy" not in lmeta.datasets.runinfo[period][run]:
                continue

            runinfo = lmeta.datasets.runinfo[period][run].phy

            if "livetime_in_s" not in runinfo:
                continue

            if run not in exposures[period]:
                exposures[period][run] = {}

            chmap = ut.generate_channel_map(runinfo.start_key, metadata=lmeta)
            data_streams_hpge = ut.select_datastreams(chmap, "HPGE",psd_ms_type=psd_ms_type)
            chmap = chmap.map("daq.rawid", unique=False)
            for det_id in data_streams_hpge:
                gedet = chmap[det_id]
                exposures[period][run][gedet.name] = (
                    gedet.production.mass_in_g
                    / 1000
                    * runinfo.livetime_in_s
                    / 60
                    / 60
                    / 24
                    / 365.25 
                )

    exposures_per_run = {}
    for period, runs in exposures.items():
        exposures_per_run[period] = {}
        for run, detectors in runs.items():
            exposures_per_run[period][run] = 0
            for det_name, exposure in detectors.items():
                exposures_per_run[period][run] += exposure
    
    exposures_per_period = {}
    for period, runs in exposures.items():
        exposures_per_period[period] = 0
        for run, detectors in runs.items():
            for det_name, exposure in detectors.items():
                exposures_per_period[period] += exposure
    
    total_exposure = 0
    for period, exposure in exposures_per_period.items():
        total_exposure += exposure

    output = {
        "exposures": exposures,
        "exposures_per_run": exposures_per_run,
        "exposures_per_period": exposures_per_period,
        "total_exposure": total_exposure,
        "psd_ms_type": psd_ms_type,
        "ignore_periods": ignore_periods
    }
    
    with open(output_path, "w") as f:
        yaml.dump(output, f, default_flow_style=False)
        
    

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Calculate exposure for HPGe detectors.")
    parser.add_argument("output", type=str, help="Output YAML file path.")
    parser.add_argument("metadata", type=str, help="Path to metadata file.")
    parser.add_argument("--psd_ms_type", type=str, default=None, help="PSD MS type to consider (default: None).")
    parser.add_argument("--ignore_periods", nargs='*', default=None, help="Periods to ignore (default: None).")
    args = parser.parse_args()

    print("Calculating exposure for HPGe detectors...")
    print(f"Output: {args.output}")
    print(f"Metadata: {args.metadata}")
    
    calculate_exposure(args.output, args.metadata, args.psd_ms_type, args.ignore_periods)