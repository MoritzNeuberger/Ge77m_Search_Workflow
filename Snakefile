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
print("Number of initial inputs found:", len(initial))

# rule all: build everything in the order given by config["workflow"]
rule all:
    input:
        expand(
            os.path.join(
                config["out_root"],
                config["workflow"][0],
                "{lvl1}",
                "{lvl2}",
                "{base}".replace("tier_pht", "tier_mgc") + ".lh5",
            ),
            **{**{k: [d[k] for d in initial] for k in ["lvl1","lvl2","base"]} }
        ),
        expand(
            os.path.join(
                config["out_root"],
                config["workflow"][1],
                "{lvl1}",
                "{lvl2}",
                "{base}".replace("tier_pht", "tier_dc") + ".lh5",
            ),
            **{**{k: [d[k] for d in initial] for k in ["lvl1","lvl2","base"]} }
        ),
        expand(
            os.path.join(
                config["out_root"],
                "make_skm",
                "{lvl1}",
                "{lvl1}-{lvl2}.lh5"
            ),
            lvl1=[x[0] for x in lvl1_lvl2_pairs],
            lvl2=[x[1] for x in lvl1_lvl2_pairs]
        )