import numpy as np
import awkward as ak
import scipy as sp
from numba_stats import poisson, norm
from tqdm import tqdm
from iminuit import Minuit
import argparse
import yaml
import lgdo.lh5 as lh5
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from functools import partial
import multiprocessing as mp
import lgdo.types as types

# Precision control: Set numpy to use high precision
np.seterr(all='raise')  # Raise exceptions on numerical errors
np.set_printoptions(precision=17)  # Display more decimal places

# Precision strategy:
# - Use longdouble (extended precision) for critical calculations like nllr
# - Use float64 for storage and I/O to maintain compatibility
# - Use numba_stats for statistical functions (faster than scipy.stats and multiprocessing-compatible)
# - Convert inputs to float64 for numba_stats functions, then results back to longdouble
# - Hybrid approach: high precision arithmetic, standard precision for statistical functions

# Use extended precision for better numerical stability
DEFAULT_FLOAT_TYPE = np.longdouble
STORAGE_FLOAT_TYPE = np.float64

print("Using longdouble precision for calculations with numba_stats (hybrid approach)")

def ensure_precision(value, dtype=DEFAULT_FLOAT_TYPE):
    """Ensure a value uses the specified precision (default for calculations)"""
    if isinstance(value, (int, float)):
        return dtype(value)
    elif isinstance(value, np.ndarray):
        return value.astype(dtype)
    return value


def nll(data, model):
        # Convert inputs to dictionaries with consistent precision
        data_dict = to_dict(data)
        model_dict = to_dict(model)
        
        # Extract values with consistent precision
        n_sb = data_dict["n_sb"]
        n_b_ext = data_dict["n_b_ext"]
        eps_sel_data = data_dict["eps_sel"]
        
        eps_sel = model_dict["eps_sel"]
        r_Ge77m = model_dict["r_Ge77m"]
        r_b = model_dict["r_b"]
        ratio = model_dict["ratio"]
        exp = model_dict["exp"]
        unc_sel = model_dict["unc_sel"]
        
        # Calculate with consistent precision
        mu_sb = (eps_sel * r_Ge77m + r_b) * exp
        mu_b_ext = (r_b / ratio) * exp
        
        # Convert to float64 for numba_stats functions, then back to high precision
        # numba_stats is faster and multiprocessing-compatible
        ll_sb = DEFAULT_FLOAT_TYPE(poisson.logpmf(float(n_sb), float(mu_sb)))
        ll_b_ext = DEFAULT_FLOAT_TYPE(poisson.logpmf(float(n_b_ext), float(mu_b_ext)))
        ll_eps_sel = DEFAULT_FLOAT_TYPE(norm.logpdf(float(eps_sel_data), loc=float(eps_sel), scale=float(unc_sel)))

        result = DEFAULT_FLOAT_TYPE(-2.0) * (ll_sb + ll_b_ext + ll_eps_sel)
        return result


def estimate_MLE(dat,set_r_Ge77m=False,verbose=False):

        def nll_wrapper(n_sb,n_b_ext,r_Ge77m,r_b,eps_sel):
                tmp_data = {
                        "n_sb": n_sb,
                        "n_b_ext": n_b_ext,
                        "eps_sel": dat["eps_sel"]
                }
                tmp_model = {
                        "r_Ge77m": r_Ge77m,
                        "r_b": r_b,
                        "eps_sel": eps_sel,
                }
                # Handle both dict and awkward array cases
                if hasattr(dat, 'keys'):
                        keys_to_check = dat.keys()
                else:
                        # For awkward arrays, get the field names
                        keys_to_check = dat.fields if hasattr(dat, 'fields') else ['ratio', 'exp', 'eps_sel', 'unc_sel']
                
                for key in keys_to_check:
                        if "ratio" in key or "exp" in key or "unc_" in key:
                                tmp_model[key] = dat[key]

                tmp = nll(tmp_data,tmp_model)

                if verbose:
                        print(f"n_sb: {n_sb}, n_b_ext: {n_b_ext}, r_Ge77m: {r_Ge77m}, r_b: {r_b}, eps_sel: {eps_sel} -> NLL: {tmp}")

                return tmp

        def manual_bf(dat):
                output = {}
                n_sb = DEFAULT_FLOAT_TYPE(dat["n_sb"])
                n_b_ext = DEFAULT_FLOAT_TYPE(dat["n_b_ext"])
                ratio = DEFAULT_FLOAT_TYPE(dat["ratio"])
                exp = DEFAULT_FLOAT_TYPE(dat["exp"])
                eps_sel = DEFAULT_FLOAT_TYPE(dat["eps_sel"])
                unc_sel = DEFAULT_FLOAT_TYPE(dat["unc_sel"])
                
                # Calculate with high precision using DEFAULT_FLOAT_TYPE
                r_Ge77m_calc = (n_sb - n_b_ext * ratio) / exp / eps_sel
                output["r_Ge77m"] = DEFAULT_FLOAT_TYPE(np.max([r_Ge77m_calc, DEFAULT_FLOAT_TYPE(0.0)]))
                
                r_b_calc = n_sb / exp - output["r_Ge77m"] * eps_sel
                output["r_b"] = DEFAULT_FLOAT_TYPE(np.max([r_b_calc, DEFAULT_FLOAT_TYPE(0.0)]))
                
                output["ratio"] = DEFAULT_FLOAT_TYPE(ratio)
                output["exp"] = DEFAULT_FLOAT_TYPE(exp)
                output["eps_sel"] = DEFAULT_FLOAT_TYPE(eps_sel)
                output["unc_sel"] = DEFAULT_FLOAT_TYPE(unc_sel)
                return output

        output = manual_bf(dat)

        if verbose:
                print("Manual BF output:", output)

        m = Minuit(nll_wrapper, 
                        n_sb = dat["n_sb"], 
                        n_b_ext = dat["n_b_ext"], 
                        r_Ge77m = output["r_Ge77m"], 
                        r_b = output["r_b"],
                        eps_sel = dat["eps_sel"]
                )

        dat_ak = ak.Array([dat])

        for key in dat_ak.fields:
                if "n_" in key:
                        m.fixed[key] = True
        for key in output.keys():
                if "r_" in key:
                        m.limits[key] = [0,100]
        for key in output.keys():
                if "eps_" in key:
                        m.limits[key] = [0,1]
        if set_r_Ge77m:
                m.values["r_Ge77m"] = set_r_Ge77m
                m.fixed["r_Ge77m"] = True

        m.errordef = Minuit.LIKELIHOOD
        
        # Set tighter tolerances for better precision
        m.tol = 1e-10  # Tighter convergence tolerance
        m.precision = 1e-12  # Higher precision for internal calculations
        
        m.migrad()

        if verbose:
                print("Minuit output:", m.params)

        output["r_Ge77m"] = DEFAULT_FLOAT_TYPE(m.params["r_Ge77m"].value)
        output["r_b"] = DEFAULT_FLOAT_TYPE(m.params["r_b"].value)
        output["eps_sel"] = DEFAULT_FLOAT_TYPE(m.params["eps_sel"].value)

        if hasattr(dat, 'keys'):
                keys_to_check = dat.keys()
        else:
                # For awkward arrays, get the field names
                keys_to_check = dat.fields if hasattr(dat, 'fields') else ['ratio', 'exp']
        
        for key in keys_to_check:
                if "ratio" in key or "exp" in key or "unc_" in key:
                        output[key] = dat[key]

        return output

def nllr(data, model):
        # Enhanced precision approach for likelihood ratio calculation
        # The key is to minimize precision loss from subtracting similar values
        
        # Convert inputs to dictionaries with consistent precision
        data_hp = to_dict(data)
        model_hp = to_dict(model)
        
        # Get the fixed r_Ge77m value with high precision
        r_Ge77m_fixed = model_hp["r_Ge77m"]
        
        # Compute MLEs with enhanced precision settings
        model_MLE_r_Ge77m_fixed = estimate_MLE(data_hp, set_r_Ge77m=r_Ge77m_fixed)
        model_MLE = estimate_MLE(data_hp)
        
        # Convert MLE results to consistent precision
        model_MLE_r_Ge77m_fixed_hp = to_dict(model_MLE_r_Ge77m_fixed)
        model_MLE_hp = to_dict(model_MLE)
        
        # Compute NLLs with careful precision handling
        nll_fixed = nll(data_hp, model_MLE_r_Ge77m_fixed_hp)
        nll_free = nll(data_hp, model_MLE_hp)
        
        # Return the difference with careful precision handling
        result = DEFAULT_FLOAT_TYPE(nll_fixed) - DEFAULT_FLOAT_TYPE(nll_free)
        
        # Ensure result is non-negative (likelihood ratio should be >= 0)
        return DEFAULT_FLOAT_TYPE(max(result, 0.0))

def generate_toy_data(model,sample_size=1000):
        def sample(input_data):
                # Ensure all inputs use highest precision
                eps_sel = DEFAULT_FLOAT_TYPE(input_data["eps_sel"])
                unc_sel = DEFAULT_FLOAT_TYPE(input_data["unc_sel"])
                r_b = DEFAULT_FLOAT_TYPE(input_data["r_b"])
                r_Ge77m = DEFAULT_FLOAT_TYPE(input_data["r_Ge77m"])
                exp = DEFAULT_FLOAT_TYPE(model["exp"])
                ratio = DEFAULT_FLOAT_TYPE(input_data["ratio"])
                
                true_eps_sel = DEFAULT_FLOAT_TYPE(np.random.normal(eps_sel, unc_sel))
                true = {
                        "n_b": np.random.poisson(r_b * exp),
                        "n_s": np.random.poisson(true_eps_sel * r_Ge77m * exp),
                        "n_b_ext": np.random.poisson(r_b / ratio * exp),
                        "eps_sel": true_eps_sel
                }
                output = {
                        "n_sb": DEFAULT_FLOAT_TYPE(true["n_s"] + true["n_b"]),
                        "n_b_ext": DEFAULT_FLOAT_TYPE(true["n_b_ext"]),
                        "ratio": ratio,
                        "exp": exp,
                        "eps_sel": eps_sel, # Use nominal value, not the true one
                        "unc_sel": unc_sel
                }
                
                return output

        return ak.Array([sample(model) for _ in range(sample_size)])

def generate_model(template_model,r_Ge77m):
        output = template_model.copy()
        output["r_Ge77m"] = r_Ge77m
        return output

def get_interval(para_range,p_value_range,alpha=0.1):
        f = sp.interpolate.interp1d(para_range, p_value_range - alpha, kind='linear')
        # Find approximate roots by evaluating on a fine grid
        x_fine = np.linspace(para_range.min(), para_range.max(), 10000)
        y_fine = f(x_fine)
        
        # Find sign changes
        roots = []
        for i in range(len(y_fine)-1):
            if y_fine[i] * y_fine[i+1] < 0:
                # Linear interpolation to get more precise root
                x1, x2 = x_fine[i], x_fine[i+1]
                y1, y2 = y_fine[i], y_fine[i+1]
                root = x1 - y1 * (x2 - x1) / (y2 - y1)
                roots.append(root)
        return np.array(roots)

def calc_p_val(thr,distr):
        return np.count_nonzero(distr >= thr)/len(distr)


def generate_toy_mcs_batch(args):
    """Helper function to compute NLL H0 for a specific model"""
    model, sample_size = args
    return generate_toy_data(model,sample_size=sample_size)

def compute_nll_for_toy_batch(args):
    """Helper function to compute NLL for a batch of toy data"""
    toy_batch, model = args
    return [nllr(data_test, model) for data_test in toy_batch]


def compute_nll_h0_for_model(args):
    """Helper function to compute NLL H0 for a specific model"""
    toy_h0, model = args
    return [nllr(toy_h0_entry, model) for toy_h0_entry in toy_h0]


def compute_p_val_for_batch(args):
    """Helper function to compute p-values for a batch"""
    nll_h0_batch, nll_batch = args
    return [calc_p_val(nll_h0_entry, nll_batch) for nll_h0_entry in nll_h0_batch]



def calculate_mle_and_ci(
        output_path,
        rc_rates_path,
        sel_eff_path,
        selected_candidates_path,
        exposure_path,
        n_threads=None,
        sample_size=10000,
        range_max=2.0,
        range_points=40
        ):

        if n_threads is None:
                n_threads = mp.cpu_count()

        with open(rc_rates_path, 'r') as f:
                rc_rates = yaml.safe_load(f)
        with open(sel_eff_path, 'r') as f:
                sel_eff = yaml.safe_load(f)
        selected_candidates = lh5.read_as("mdc",selected_candidates_path,"ak")
        with open(exposure_path, 'r') as f:
                exposure = yaml.safe_load(f)

        n_obs = len(selected_candidates)
        n_b_ext = rc_rates["selection"]["mult_lar_psd"]["cts"]
        ratio = rc_rates["scaling"]
        exp = exposure["total_exposure"]

        selection_effs = [sel_eff["eps_time"], sel_eff["eps_prompt"], sel_eff["eps_delayed"]["selection"]["mult_lar_psd"]]

        # Calculate with highest precision
        eps_values = [DEFAULT_FLOAT_TYPE(eps["val"]) for eps in selection_effs]
        eps_sel = DEFAULT_FLOAT_TYPE(np.prod(eps_values))
        
        # Uncertainty propagation with highest precision
        unc_terms = []
        for i in range(len(selection_effs)):
                unc_i = DEFAULT_FLOAT_TYPE(selection_effs[i]["unc"])
                prod_others = DEFAULT_FLOAT_TYPE(np.prod([eps_values[j] for j in range(len(eps_values)) if j != i]))
                unc_terms.append((unc_i * prod_others) ** 2)
        
        unc_sel = DEFAULT_FLOAT_TYPE(np.sqrt(np.sum(unc_terms)))

        data = {
                "n_sb": DEFAULT_FLOAT_TYPE(n_obs),
                "n_b_ext": DEFAULT_FLOAT_TYPE(n_b_ext),
                "ratio": DEFAULT_FLOAT_TYPE(ratio),
                "exp": DEFAULT_FLOAT_TYPE(exp),
                "eps_sel": eps_sel,
                "unc_sel": unc_sel,
        }
        print("Data: ", data)

        model_MLE = estimate_MLE(data)
        print("Estimated MLE parameters:", model_MLE)

        para_range = np.linspace(0, range_max, range_points, dtype=DEFAULT_FLOAT_TYPE)
        bin_width = para_range[1] - para_range[0]
        min_val = np.min([model_MLE["r_Ge77m"],0.001* bin_width])
        para_range[para_range<min_val] = min_val  # Avoid numerical issues with very small values

        model_range = ak.Array([generate_model(model_MLE,x) for x in para_range])

        print("Generating toy data...")
        with ProcessPoolExecutor(max_workers=n_threads) as executor:
                # Prepare arguments for parallel processing
                args_list = [(model_range[i], sample_size) for i in range(len(model_range))]
                
                # Use list comprehension with executor.map for parallel processing
                toy_range = list(tqdm(executor.map(generate_toy_mcs_batch, args_list), 
                                    total=len(args_list), desc="toy generation"))
        

        print("Computing NLL range using multithreading...")
        with ProcessPoolExecutor(max_workers=n_threads) as executor:
                # Prepare arguments for parallel processing
                args_list = [(toy_range[i], model_range[i]) for i in range(len(toy_range))]
                
                # Use list comprehension with executor.map for parallel processing
                nll_range = list(tqdm(executor.map(compute_nll_for_toy_batch, args_list), 
                                    total=len(args_list), desc="NLL calculations"))
        
        nll_range = np.array(nll_range)
        data_nll = np.array([nllr(data,model_range[i]) for i in range(len(toy_range))])
        p_value_range = np.array([np.count_nonzero(nll_range[i] >= data_nll[i])/len(nll_range[i]) for i in range(len(toy_range))])

        model_H0 = model_MLE.copy()
        model_H0["r_Ge77m"] = 0

        toy_H0 = generate_toy_data(model_H0,sample_size=sample_size)

        print("Computing NLL H0 range using multithreading...")
        with ProcessPoolExecutor(max_workers=n_threads) as executor:
                args_list = [(toy_H0, model_range[i]) for i in range(len(toy_range))]
                nll_H0_range = list(tqdm(executor.map(compute_nll_h0_for_model, args_list), 
                                       total=len(args_list), desc="NLL H0 calculations"))
        
        nll_H0_range = np.array(nll_H0_range)

        print("Computing p-value H0 using multithreading...")
        with ProcessPoolExecutor(max_workers=n_threads) as executor:
                args_list = [(nll_H0_range[i], nll_range[i]) for i in range(len(nll_range))]
                p_val_H0 = list(tqdm(executor.map(compute_p_val_for_batch, args_list), 
                                   total=len(args_list), desc="p-value H0 calculations"))
        
        p_val_H0 = np.array(p_val_H0)

        median_H0 = np.median(p_val_H0,axis=1)
        bound_68_lower = np.array([np.sort(p_val_H0_entry)[int(len(p_val_H0_entry)*((1 - 0.68)/2.))] for p_val_H0_entry in p_val_H0])
        bound_68_higher = np.array([np.sort(p_val_H0_entry)[int(len(p_val_H0_entry)*(1 - (1 - 0.68)/2.))] for p_val_H0_entry in p_val_H0])
        bound_90_lower = np.array([np.sort(p_val_H0_entry)[int(len(p_val_H0_entry)*(0.10/2.))] for p_val_H0_entry in p_val_H0])
        bound_90_higher = np.array([np.sort(p_val_H0_entry)[int(len(p_val_H0_entry)*(1 - 0.10/2.))] for p_val_H0_entry in p_val_H0])

        import json
        
        # Convert all output values to JSON-serializable format (float64)
        output = ensure_precision_dict({
                "MLE": model_MLE,
                "para": list(para_range),
                "obs":  list(p_value_range),
                "exp":  list(median_H0),
                "m_exp": list(bound_68_lower),
                "p_exp": list(bound_68_higher),
                "m_exp_90": list(bound_90_lower),
                "p_exp_90": list(bound_90_higher)
        }, dtype=STORAGE_FLOAT_TYPE)  # Ensure float64 for JSON compatibility
        
        # Use high precision when writing to JSON
        class HighPrecisionEncoder(json.JSONEncoder):
            def encode(self, obj):
                if isinstance(obj, np.floating):
                    return format(float(obj), '.17g')  # Convert to Python float first
                return super().encode(obj)
                
            def default(self, obj):
                if isinstance(obj, np.floating):
                    return float(obj)  # Convert numpy floats to Python floats
                elif isinstance(obj, np.integer):
                    return int(obj)    # Convert numpy ints to Python ints
                elif isinstance(obj, np.ndarray):
                    return obj.tolist()  # Convert arrays to lists
                return super().default(obj)
        
        with open(output_path,"w") as f:
                json.dump(output, f, indent=4, cls=HighPrecisionEncoder)

        
        # write nll_range, data_nll and nll_H0 to lh5 file
        import os
        lh5_path = output_path.replace(".json", ".lh5")
        
        # Remove existing file to avoid shape conflicts
        if os.path.exists(lh5_path):
                os.remove(lh5_path)
        
        # Ensure arrays have consistent dtypes for storage (convert from longdouble to float64)
        nll_range = np.array(nll_range, dtype=STORAGE_FLOAT_TYPE)
        data_nll = np.array(data_nll, dtype=STORAGE_FLOAT_TYPE)
        nll_H0_range = np.array(nll_H0_range, dtype=STORAGE_FLOAT_TYPE)
        para_range = np.array(para_range, dtype=STORAGE_FLOAT_TYPE)
        
        output_lh5 = {
                "nll_range": types.Array(nll_range),
                "data_nll": types.Array(data_nll),
                "nll_H0_range": types.Array(nll_H0_range),
                "para_range": types.Array(para_range)
        }
        
        print("LH5 output shapes:")
        print(f"nll_range: {nll_range.shape}")
        print(f"data_nll: {data_nll.shape}")  
        print(f"nll_H0_range: {nll_H0_range.shape}")
        print(f"para_range: {para_range.shape}")

        lh5.write(types.Struct(output_lh5), name="nll", lh5_file=lh5_path)


def ensure_precision_dict(data_dict, dtype=STORAGE_FLOAT_TYPE):
    """Recursively ensure all numeric values in a dictionary use specified precision for storage"""
    result = {}
    for key, value in data_dict.items():
        if isinstance(value, dict):
            result[key] = ensure_precision_dict(value, dtype)
        elif isinstance(value, list):
            converted_list = []
            for v in value:
                if isinstance(v, (int, float, np.number)):
                    converted = ensure_precision(v, dtype)
                    if dtype == STORAGE_FLOAT_TYPE:
                        converted_list.append(float(converted))  # Convert to Python float for JSON
                    else:
                        converted_list.append(converted)
                else:
                    converted_list.append(v)
            result[key] = converted_list
        elif isinstance(value, (int, float, np.number)):
            # Convert to the specified dtype, then to Python native type for JSON compatibility
            converted = ensure_precision(value, dtype)
            if dtype == STORAGE_FLOAT_TYPE:
                result[key] = float(converted)  # Convert to Python float for JSON
            else:
                result[key] = converted
        else:
            result[key] = value
    return result

def to_dict(obj):
    """Convert awkward array or dict to a regular dictionary with high precision values"""
    result = {}
    if hasattr(obj, 'items'):
        # It's already a dictionary
        for key, value in obj.items():
            result[key] = DEFAULT_FLOAT_TYPE(value)
    elif hasattr(obj, 'fields'):
        # It's an awkward array
        for field in obj.fields:
            result[field] = DEFAULT_FLOAT_TYPE(obj[field])
    else:
        # Try to iterate over it as if it's a mapping
        try:
            for key in obj:
                result[key] = DEFAULT_FLOAT_TYPE(obj[key])
        except:
            raise ValueError(f"Cannot convert object of type {type(obj)} to dictionary")
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate MLE and confidence intervals for a given dataset.")
    parser.add_argument("output_path", type=str, help="Path to the output file where results will be saved.")
    parser.add_argument("rc_rates_path", type=str, help="Path to the YAML file containing RC rates.")
    parser.add_argument("sel_eff_path", type=str, help="Path to the YAML file containing selection efficiency.")
    parser.add_argument("selected_candidates_path", type=str, help="Path to the LH5 file containing selected candidates.")
    parser.add_argument("exposure_path", type=str, help="Path to the YAML file containing exposure data.")
    parser.add_argument("--n_threads", type=int, default=None, help="Number of threads to use for parallel processing (default: all available cores).")
    parser.add_argument("--sample_size", type=int, default=10000, help="Number of samples for each model (default: 1000).")
    parser.add_argument("--range_max", type=float, default=2.0, help="Maximum range for the parameter (default: 2.0).")
    parser.add_argument("--range_points", type=int, default=40, help="Number of points in the parameter range (default: 40).")

    args = parser.parse_args()

    print("Starting MLE and CI calculation with the following parameters:")
    print(f"Output Path: {args.output_path}")
    print(f"RC Rates Path: {args.rc_rates_path}")
    print(f"Selection Efficiency Path: {args.sel_eff_path}")
    print(f"Selected Candidates Path: {args.selected_candidates_path}")
    print(f"Exposure Path: {args.exposure_path}")
    print(f"Number of Threads: {args.n_threads}")
    print(f"Sample Size: {args.sample_size}")

    calculate_mle_and_ci(
        args.output_path,
        args.rc_rates_path,
        args.sel_eff_path,
        args.selected_candidates_path,
        args.exposure_path,
        n_threads=args.n_threads,
        sample_size=args.sample_size,
        range_max=args.range_max,
        range_points=args.range_points
    )