# Snakefile

import os
import glob
import re

configfile: "config.yaml"

# include all rule files
include: "rules/mu_hpge_coinc.smk"
include: "rules/delayed_coinc.smk"
include: "rules/make_skm.smk"

# discover all initial inputs
# They live under config["input_root"]/*/*/*.lh5
initial = []
for lvl1 in os.listdir(config["input_root"]):
    p1 = os.path.join(config["input_root"], lvl1)
    if os.path.isdir(p1):
        for lvl2 in os.listdir(p1):
            p2 = os.path.join(p1, lvl2)
            if os.path.isdir(p2):
                for f in glob.glob(os.path.join(p2, "*.lh5")):
                    base = os.path.basename(f).rsplit(".lh5", 1)[0]
                    # Filter if target_file is set
                    if config.get("target_file") and base != config["target_file"].rsplit(".lh5", 1)[0]:
                        continue
                    initial.append(dict(lvl1=lvl1, lvl2=lvl2, base=base))

# Collect all output files in a single list
all_outputs = []

for d in initial:
    # mu_hpge_coinc output
    all_outputs.append(
        os.path.join(
            config["out_root"],
            config["workflow"][0],
            d["lvl1"],
            d["lvl2"],
            d["base"].replace("tier_pht", "tier_mgc") + ".lh5"
        )
    )
    # delayed_coinc output
    all_outputs.append(
        os.path.join(
            config["out_root"],
            config["workflow"][1],
            d["lvl1"],
            d["lvl2"],
            d["base"].replace("tier_pht", "tier_dc") + ".lh5"
        )
    )
    # make_skm output
    all_outputs.append(
        os.path.join(
            config["out_root"],
            "make_skm",
            d["lvl1"],
            f"{d['lvl1']}-{d['lvl2']}.lh5"
        )
    )

# rule all: build everything in the order given by config["workflow"]
rule all:
    input:
        all_outputs