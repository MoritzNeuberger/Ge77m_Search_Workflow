# Snakefile

import os
import glob
import re

configfile: "config.yaml"

input_root = config["input_root"].format(default_ref_version=config["default_ref_version"])

# include all rule files
include: "rules/mu_hpge_coinc.smk"
include: "rules/make_mgc_skm.smk"
include: "rules/mu_delayed_coinc.smk"
include: "rules/make_skm.smk"
include: "rules/draw_mu_diff_vs_energy.smk"
include: "rules/draw_mgc_candidates.smk"
include: "rules/check_muon_tagging.smk"

# discover all initial inputs
# They live under config["input_root"]/*/*/*.lh5

initial = []
for lvl1 in os.listdir(input_root):
    if lvl1 in config["ignore_periods"]:
        continue  # Skip ignored periods
    p1 = os.path.join(input_root, lvl1)
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
all_inputs = []

muon_tagging = ["gen/muon_tagging_stats.txt"]
mgc_outputs = []
mdc_outputs = []

for d in initial:
    mgc_outputs.append(
        os.path.join(
            config["out_root"],
            config["workflow"][0],
            d["lvl1"],
            d["lvl2"],
            d["base"].replace("tier_pht", "tier_mgc") + ".lh5"
        )
    )
    mdc_outputs.append(
        os.path.join(
            config["out_root"],
            config["workflow"][1],
            d["lvl1"],
            d["lvl2"],
            d["base"].replace("tier_pht", "tier_mdc") + ".lh5"
        )
    )

skm_mgc_output = [
    os.path.join(
        config["out_root"],
        config["workflow"][0],
        "skm_mgc.lh5"
    )
]

waveform_block = "/tmp/mgc_plots_all.done"

draw_plots = [
    "gen/figs/mu_diff_vs_energy.pdf",
    waveform_block
]

os.makedirs(os.path.dirname(draw_plots[0]), exist_ok=True)

all_inputs = sorted(set(muon_tagging + mgc_outputs + skm_mgc_output + draw_plots + mdc_outputs))  # Remove duplicates and sort

# rule all: build everything in the order given by config["workflow"]
rule all:
    input:
        all_inputs

rule mgc_all:
    input:
        mgc_outputs

rule mdc_all:
    input:
        mgc_outputs,
        mdc_outputs


barrier_file = "/tmp/mgc_all.done"

rule mgc_barrier:
    input:
        mgc_outputs
    output:
        barrier_file
    shell:
        "touch {output}"
        
    
    # This ensures mdc_all only runs after mgc_all is complete
    # Use 'run' or 'shell' as appropriate for your workflow
    # Here, we just use input/output to enforce order

# Make mdc_all depend on mgc_all
ruleorder: mgc_all > mdc_all
