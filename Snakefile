# Snakefile

import os
import glob

configfile: "config.yaml"

input_root = config["input_root"].format(default_ref_version=config["default_ref_version"])

# --- Include all rule files ---
include: "rules/prod/build_mgc.smk"
include: "rules/prod/make_mgc_sum.smk"
include: "rules/prod/build_mdc.smk"
include: "rules/prod/make_mdc_sum.smk"
include: "rules/prod/build_skm.smk"
include: "rules/prod/calculate_exposure.smk"
include: "rules/prod/calculate_rc_bkg.smk"
include: "rules/prod/calculate_selection_eff.smk"
include: "rules/prod/calculate_ideal_main_selection.smk"
include: "rules/prod/select_candidates.smk"
include: "rules/prod/run_frequentist_calculator.smk"

include: "rules/xchecks_and_drawers/draw_mu_diff_vs_energy.smk"
#include: "rules/xchecks_and_drawers/draw_mgc_candidates.smk"
include: "rules/xchecks_and_drawers/check_muon_tagging.smk"
include: "rules/xchecks_and_drawers/draw_mgc_rate.smk"
include: "rules/xchecks_and_drawers/check_saturation_sim_vs_data.smk"
include: "rules/xchecks_and_drawers/draw_mgc_non_physical_vs_physical_distr.smk"
include: "rules/xchecks_and_drawers/draw_mgc_data_vs_sim.smk"
include: "rules/xchecks_and_drawers/check_frequentist_calculator.smk"

# --- Helper functions for output paths ---
def mgc_output_path(d):
    return os.path.join(
        config["out_root"],
        config["workflow"][0],
        d["lvl1"],
        d["lvl2"],
        d["base"].replace("tier_pht", "tier_mgc") + ".lh5"
    )

def mdc_output_path(d):
    tmp = os.path.join(
        config["out_root"],
        config["workflow"][1],
        d["lvl1"],
        d["lvl2"],
        d["base"].replace("tier_pht", "tier_mdc") + ".lh5"
    )
    # Always return the path, regardless of whether file exists
    return tmp 

# --- Discover all initial inputs ---
# They live under config["input_root"]/*/*/*.lh5
initial = []
for lvl1 in os.listdir(input_root):
    if lvl1 in config["ignore_periods"]:
        continue
    p1 = os.path.join(input_root, lvl1)
    if os.path.isdir(p1):
        for lvl2 in os.listdir(p1):
            p2 = os.path.join(p1, lvl2)
            if os.path.isdir(p2):
                for f in glob.glob(os.path.join(p2, "*.lh5")):
                    base = os.path.basename(f).rsplit(".lh5", 1)[0]
                    if config.get("target_file") and base != config["target_file"].rsplit(".lh5", 1)[0]:
                        continue
                    initial.append(dict(lvl1=lvl1, lvl2=lvl2, base=base))

# --- Output file lists ---
muon_tagging = ["gen/muon_tagging_stats.txt"]
mgc_outputs = [mgc_output_path(d) for d in initial]
mdc_outputs = [mdc_output_path(d) for d in initial]

# Generate sum_mgc files for each period/run combination found in the initial data
sum_mgc_per_run_outputs = []
periods_runs = set()
for d in initial:
    periods_runs.add((d["lvl1"], d["lvl2"]))

for period, run in periods_runs:
    sum_mgc_per_run_outputs.append(f"gen/mu_hpge_coinc/{period}/sum_mgc_{run}.lh5")

sum_mgc_output = [os.path.join(config["out_root"], config["workflow"][0], "sum_mgc.lh5")]
sum_mdc_output = [os.path.join(config["out_root"], config["workflow"][1], "sum_mdc.lh5")]

draw_plots = [
    "gen/figs/mu_diff_vs_energy.pdf",
    "gen/figs/mgc_rate.pdf",
#    "/tmp/mgc_plots_all.done",
    "gen/muon_saturation_sim_vs_data.yaml",
    "gen/figs/mu_non_physical_vs_physical_distr_energy_lin.pdf",
    "gen/figs/mgc_data_vs_sim_energy_log.pdf",
    "input/freq_tests/only_signal_signal_eff_one/output.lh5",
]
os.makedirs(os.path.dirname(draw_plots[0]), exist_ok=True)

exposure_files = [
    f"gen/exposure/exposures_{'-'.join(psd) if isinstance(psd, list) else psd}.yaml"
    for psd in config["psd_ms_type"]
]

skm_files = [
    f"gen/skm/skm_{'-'.join(psd) if isinstance(psd, list) else psd}.lh5"
    for psd in config["psd_ms_type"]
]


rc_bkg_files = [
    f"gen/rc_bkg/rates_{'-'.join(psd) if isinstance(psd, list) else psd}.yaml"
    for psd in config["psd_ms_type"]
]

efficiency_file = [ "gen/sim/selection_eff.yaml", "gen/ideal_selection.yaml" ]

select_candidates_file = ["gen/selected_candidates.lh5"]

freq_files = [
        f"gen/freq_calc/{state}-{sys}.lh5"
    for state in ["Ge77m", "Ge77_and_Ge77m"]
    for sys in ["w_sys", "wo_sys"]
]

# --- Collect all output files in a single list ---
all_inputs = sorted(set(filter(None, 
    muon_tagging +
    mgc_outputs +
    sum_mgc_per_run_outputs +
    sum_mgc_output +
    draw_plots +
    mdc_outputs +
    sum_mdc_output +
    exposure_files + 
    skm_files + 
    rc_bkg_files + 
    efficiency_file + 
    select_candidates_file + 
    freq_files
)))

# --- Rules ---
rule all:
    input:
        all_inputs

#rule mgc_all:
#    input:
#        mgc_outputs
#
#rule mdc_all:
#    input:
#        mgc_outputs,
#        mdc_outputs

rule mgc_barrier:
    input:
        mgc_outputs
    output:
        "/tmp/mgc_all.done"
    shell:
        "touch {output}"

rule mdc_barrier:
    input:
        mdc_outputs
    output:
        "/tmp/mdc_all.done"
    shell:
        "touch {output}"

#ruleorder: mgc_all > mdc_all
