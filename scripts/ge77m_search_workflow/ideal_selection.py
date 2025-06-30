import yaml
import numpy as np
import awkward as ak

def fom(s,b):
    b[b == 0] = 1e-10  # Avoid division by zero
    return np.sqrt(2 * ((s + b) * np.log(1 + s / b) - s))

def find_ideal_given_rates(rc_rate_scans, s_eff_scans, other_eff, upper, lower, prod_rate=0.21, rel_prod_ratio=0.5, final_exp=1e3):
    signal_rate_scans = prod_rate * rel_prod_ratio * other_eff * np.array(s_eff_scans)

    x = np.array(upper)
    y = np.array(lower)

    x_mask = x > 2100
    y_mask = y <= 2000

    pre_select_mask = x_mask & y_mask

    foms = fom(signal_rate_scans[pre_select_mask] * final_exp, rc_rate_scans[pre_select_mask] * final_exp)
    args_max = np.unravel_index(np.argmax(foms, axis=None), foms.shape)

    x_max = x[pre_select_mask][args_max]
    y_max = y[pre_select_mask][args_max]

    return {
        "emax": float(x_max),
        "emin": float(y_max),
        "fom_max": float(foms[args_max]),
        "signal_rate": float(signal_rate_scans[pre_select_mask][args_max]),
        "rc_rate": float(rc_rate_scans[pre_select_mask][args_max]),
        "signal": float(signal_rate_scans[pre_select_mask][args_max] * final_exp),
        "rc": float(rc_rate_scans[pre_select_mask][args_max] * final_exp),
    }

def ideal_selection(rc_rate_file, selection_eff_file, output_file):
    
    with open(rc_rate_file, "r") as f:
        rc_rates = yaml.safe_load(f)

    with open(selection_eff_file, "r") as f:
        selection_eff = yaml.safe_load(f)

    result = {}

    for key in rc_rates["scans"].keys():
        if key == "x" or key == "y":
            continue

        rc_rate_scans = np.array(rc_rates["scans"][key])
        s_eff_scans = np.array(selection_eff["eps_delayed"]["scan"][key])
        other_eff = selection_eff["eps_prompt"]["val"] * selection_eff["eps_time"]["val"]

        result[key]= find_ideal_given_rates(
            rc_rate_scans=rc_rate_scans,
            s_eff_scans=s_eff_scans,
            other_eff=other_eff,
            upper=np.array(rc_rates["scans"]["x"]),
            lower=np.array(rc_rates["scans"]["y"])
        )

    with open(output_file, "w") as f:
        yaml.dump(result, f)
    print(f"Results written to {output_file}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Calculate ideal selection based on RC rates and selection efficiency.")
    parser.add_argument("--rc_rate_file", type=str, required=True, help="Path to the RC rate file (YAML format).")
    parser.add_argument("--selection_eff_file", type=str, required=True, help="Path to the selection efficiency file (YAML format).")
    parser.add_argument("--output_file", type=str, required=True, help="Path to save the output results (YAML format).")

    args = parser.parse_args()

    ideal_selection(args.rc_rate_file, args.selection_eff_file, args.output_file)