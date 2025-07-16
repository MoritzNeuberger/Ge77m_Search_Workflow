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


def nll(data, model, verbose=False):
        # Convert inputs to dictionaries with consistent precision
        data_dict = to_dict(data)
        model_dict = to_dict(model)
        
        # Extract values with consistent precision
        n_sb = data_dict["n_sb"]
        n_b_ext = data_dict["n_b_ext"]
        eps_sel_data = data_dict["eps_sel"]
        
        eps_sel = model_dict["eps_sel"]
        r_b = model_dict["r_b"]
        ratio = model_dict["ratio"]
        exp = model_dict["exp"]
        unc_sel = model_dict["unc_sel"]
        
        # Handle flexible parameter choice
        param_type = model_dict["param_type"]  # Default to original behavior
        scaling_factor = model_dict["scaling_factor"]
        unc_scaling = model_dict["unc_scaling"]
        scaling_factor_data = data_dict["scaling_factor"]
        
        if param_type == "r_Ge77m":
            # Original behavior: r_Ge77m is the primary parameter
            r_Ge77m = model_dict["r_Ge77m"]
            # Calculate derived r_Ge77_and_Ge77m: r_Ge77_and_Ge77m = r_Ge77m / scaling_factor
            r_Ge77_and_Ge77m = model_dict["r_Ge77_and_Ge77m"]
        else:  # param_type == "r_Ge77_and_Ge77m"
            # New behavior: r_Ge77_and_Ge77m is the primary parameter
            r_Ge77_and_Ge77m = model_dict["r_Ge77_and_Ge77m"]
            # Calculate r_Ge77m from r_Ge77_and_Ge77m: r_Ge77m = r_Ge77_and_Ge77m * scaling_factor
            r_Ge77m = r_Ge77_and_Ge77m * scaling_factor
        
        # Calculate with consistent precision
        mu_sb = (eps_sel * r_Ge77m + r_b) * exp
        mu_b_ext = (r_b / ratio) * exp
        
        # Convert to float64 for numba_stats functions, then back to high precision
        ll_sb = DEFAULT_FLOAT_TYPE(poisson.logpmf(float(n_sb), float(mu_sb)))
        ll_b_ext = DEFAULT_FLOAT_TYPE(poisson.logpmf(float(n_b_ext), float(mu_b_ext)))
        ll_eps_sel = DEFAULT_FLOAT_TYPE(norm.logpdf(float(eps_sel_data), loc=float(eps_sel), scale=float(unc_sel)))
        
        # Add constraint for scaling factor when using r_Ge77_and_Ge77m as primary parameter
        ll_scaling_factor = DEFAULT_FLOAT_TYPE(0.0)
        if param_type == "r_Ge77_and_Ge77m":
            ll_scaling_factor = DEFAULT_FLOAT_TYPE(norm.logpdf(float(scaling_factor_data), loc=float(scaling_factor), scale=float(unc_scaling)))

        if verbose:
                print(f"nll>: param_type: {param_type}")
                print(f"nll>: n_sb: {n_sb}, n_b_ext: {n_b_ext}, eps_sel_data: {eps_sel_data}")
                print(f"nll>: r_Ge77m: {r_Ge77m}, r_Ge77_and_Ge77m: {r_Ge77_and_Ge77m}")
                print(f"nll>: scaling_factor: {scaling_factor}, scaling_factor_data: {scaling_factor_data}")
                print(f"nll>: mu_sb: {mu_sb}, mu_b_ext: {mu_b_ext}, eps_sel: {eps_sel}, unc_sel: {unc_sel}")
                print(f"nll>: ll_sb: {ll_sb}, ll_b_ext: {ll_b_ext}, ll_eps_sel: {ll_eps_sel}, ll_scaling_factor: {ll_scaling_factor}")

        result = DEFAULT_FLOAT_TYPE(-2.0) * (ll_sb + ll_b_ext + ll_eps_sel + ll_scaling_factor)
        return result


def estimate_MLE(dat, set_param_value=False, verbose=False, fix_systematics=False, param_type="r_Ge77m"):

        # Flexible wrapper for different parameter types
        if param_type == "r_Ge77m":
            def nll_wrapper(n_sb, n_b_ext, r_Ge77m, r_b, eps_sel, scaling_factor=0.5):
                tmp_data = {
                        "n_sb": n_sb,
                        "n_b_ext": n_b_ext,
                        "eps_sel": dat["eps_sel"],
                        "scaling_factor": dat["scaling_factor"]
                }
                tmp_model = {
                        "r_Ge77m": r_Ge77m,
                        "r_Ge77_and_Ge77m": r_Ge77m / scaling_factor,  # Derived parameter
                        "r_b": r_b,
                        "eps_sel": eps_sel,
                        "scaling_factor": scaling_factor,
                        "param_type": param_type
                }
                # Handle both dict and awkward array cases
                if hasattr(dat, 'keys'):
                        keys_to_check = dat.keys()
                else:
                        keys_to_check = dat.fields if hasattr(dat, 'fields') else ['ratio', 'exp', 'eps_sel', 'unc_sel', 'unc_scaling']
                
                for key in keys_to_check:
                        if "ratio" in key or "exp" in key or "unc_" in key:
                                tmp_model[key] = dat[key]

                tmp = nll(tmp_data, tmp_model, verbose=verbose)
                if verbose:
                        print(f"nll_wrapper (r_Ge77m): r_Ge77m={r_Ge77m}, r_b={r_b}, eps_sel={eps_sel}, scaling_factor={scaling_factor} -> NLL: {tmp}")
                return tmp
                
        else:  # param_type == "r_Ge77_and_Ge77m"
            def nll_wrapper(n_sb, n_b_ext, r_Ge77_and_Ge77m, r_b, eps_sel, scaling_factor=0.5):
                tmp_data = {
                        "n_sb": n_sb,
                        "n_b_ext": n_b_ext,
                        "eps_sel": dat["eps_sel"],
                        "scaling_factor": dat["scaling_factor"]
                }
                tmp_model = {
                        "r_Ge77_and_Ge77m": r_Ge77_and_Ge77m,
                        "r_Ge77m": r_Ge77_and_Ge77m * scaling_factor,  # Derived parameter
                        "r_b": r_b,
                        "eps_sel": eps_sel,
                        "scaling_factor": scaling_factor,
                        "param_type": param_type
                }
                # Handle both dict and awkward array cases
                if hasattr(dat, 'keys'):
                        keys_to_check = dat.keys()
                else:
                        keys_to_check = dat.fields if hasattr(dat, 'fields') else ['ratio', 'exp', 'eps_sel', 'unc_sel', 'unc_scaling']
                
                for key in keys_to_check:
                        if "ratio" in key or "exp" in key or "unc_" in key:
                                tmp_model[key] = dat[key]

                tmp = nll(tmp_data, tmp_model, verbose=verbose)
                if verbose:
                        print(f"nll_wrapper (r_Ge77_and_Ge77m): r_Ge77_and_Ge77m={r_Ge77_and_Ge77m}, r_b={r_b}, eps_sel={eps_sel}, scaling_factor={scaling_factor} -> NLL: {tmp}")
                return tmp

        def manual_bf(dat, param_type):
                output = {}
                n_sb = DEFAULT_FLOAT_TYPE(dat["n_sb"])
                n_b_ext = DEFAULT_FLOAT_TYPE(dat["n_b_ext"])
                ratio = DEFAULT_FLOAT_TYPE(dat["ratio"])
                exp = DEFAULT_FLOAT_TYPE(dat["exp"])
                eps_sel = DEFAULT_FLOAT_TYPE(dat["eps_sel"])
                unc_sel = DEFAULT_FLOAT_TYPE(dat["unc_sel"])
                
                # Get scaling factor parameters
                scaling_factor = DEFAULT_FLOAT_TYPE(dat["scaling_factor"])
                unc_scaling = DEFAULT_FLOAT_TYPE(dat["unc_scaling"])

                # Calculate total expected rate from observed data
                total_rate_calc = (n_sb - n_b_ext * ratio) / exp / eps_sel
                total_rate_calc = DEFAULT_FLOAT_TYPE(np.max([total_rate_calc, DEFAULT_FLOAT_TYPE(0.0)]))
                
                output["r_Ge77m"] = total_rate_calc
                output["r_Ge77_and_Ge77m"] = total_rate_calc / scaling_factor
                
                r_b_calc = n_sb / exp - total_rate_calc * eps_sel
                output["r_b"] = DEFAULT_FLOAT_TYPE(np.max([r_b_calc, DEFAULT_FLOAT_TYPE(0.0)]))
                
                output["ratio"] = DEFAULT_FLOAT_TYPE(ratio)
                output["exp"] = DEFAULT_FLOAT_TYPE(exp)
                output["eps_sel"] = DEFAULT_FLOAT_TYPE(eps_sel)
                output["unc_sel"] = DEFAULT_FLOAT_TYPE(unc_sel)
                output["scaling_factor"] = scaling_factor
                output["unc_scaling"] = unc_scaling
                output["param_type"] = param_type
                return output

        output = manual_bf(dat, param_type)

        if verbose:
                print("estimate_MLE> Manual BF output:", output)

        # Initialize Minuit with appropriate parameters based on param_type
        if param_type == "r_Ge77m":
                m = Minuit(nll_wrapper, 
                        n_sb = dat["n_sb"], 
                        n_b_ext = dat["n_b_ext"], 
                        r_Ge77m = output["r_Ge77m"], 
                        r_b = output["r_b"],
                        eps_sel = dat["eps_sel"],
                        scaling_factor = output["scaling_factor"]
                )
                primary_param = "r_Ge77m"
        else:  # param_type == "r_Ge77_and_Ge77m"
                m = Minuit(nll_wrapper, 
                        n_sb = dat["n_sb"], 
                        n_b_ext = dat["n_b_ext"], 
                        r_Ge77_and_Ge77m = output["r_Ge77_and_Ge77m"], 
                        r_b = output["r_b"],
                        eps_sel = dat["eps_sel"],
                        scaling_factor = output["scaling_factor"]
                )
                primary_param = "r_Ge77_and_Ge77m"

        dat_ak = ak.Array([dat])

        for key in dat_ak.fields:
                if "n_" in key:
                        m.fixed[key] = True
        
        # Set limits for rate parameters

        m.limits[primary_param] = [0, 100]
        m.limits["r_b"] = [0, 100]
        
        # Set limits for efficiency parameters
        m.limits["eps_sel"] = [0, 1]
        if fix_systematics:
                m.fixed["eps_sel"] = True
        
        # Set limits and constraints for scaling factor
        m.limits["scaling_factor"] = [0.1, 1.0]  # Reasonable range around 0.5
        if fix_systematics:
                m.fixed["scaling_factor"] = True

        # Handle parameter fixing
        if set_param_value:
                m.values[primary_param] = set_param_value
                m.fixed[primary_param] = True

        m.errordef = Minuit.LIKELIHOOD
        
        # Set tighter tolerances for better precision
        m.tol = 1e-10  # Tighter convergence tolerance
        m.precision = 1e-12  # Higher precision for internal calculations
        
        m.migrad()

        if verbose:
                print("estimate_MLE> Minuit output:", m.params)

        # Extract fitted parameters
        if param_type == "r_Ge77m":
                output["r_Ge77m"] = DEFAULT_FLOAT_TYPE(m.params["r_Ge77m"].value)
                output["r_Ge77_and_Ge77m"] = output["r_Ge77m"] / DEFAULT_FLOAT_TYPE(m.params["scaling_factor"].value)
        else:
                output["r_Ge77_and_Ge77m"] = DEFAULT_FLOAT_TYPE(m.params["r_Ge77_and_Ge77m"].value)
                # Calculate r_Ge77m from fitted values
                output["r_Ge77m"] = output["r_Ge77_and_Ge77m"] * DEFAULT_FLOAT_TYPE(m.params["scaling_factor"].value)
        
        output["r_b"] = DEFAULT_FLOAT_TYPE(m.params["r_b"].value)
        output["eps_sel"] = DEFAULT_FLOAT_TYPE(m.params["eps_sel"].value)
        output["scaling_factor"] = DEFAULT_FLOAT_TYPE(m.params["scaling_factor"].value)

        if hasattr(dat, 'keys'):
                keys_to_check = dat.keys()
        else:
                # For awkward arrays, get the field names
                keys_to_check = dat.fields if hasattr(dat, 'fields') else ['ratio', 'exp']
        
        for key in keys_to_check:
                if "ratio" in key or "exp" in key or "unc_" in key:
                        output[key] = dat[key]

        return output

def nllr(data, model, verbose=False, fix_systematics=False, param_type="r_Ge77m"):
        # Enhanced precision approach for likelihood ratio calculation
        # The key is to minimize precision loss from subtracting similar values
        
        # Convert inputs to dictionaries with consistent precision
        data_hp = to_dict(data)
        model_hp = to_dict(model)
        
        # Get the fixed parameter value based on param_type
        if param_type == "r_Ge77m":
            param_fixed = model_hp["r_Ge77m"]
        else:  # param_type == "r_Ge77_and_Ge77m"
            param_fixed = model_hp["r_Ge77_and_Ge77m"]
        
        # Compute MLEs with enhanced precision settings
        model_MLE_param_fixed = estimate_MLE(data_hp, set_param_value=param_fixed, verbose=verbose, fix_systematics=fix_systematics, param_type=param_type)
        model_MLE = estimate_MLE(data_hp, verbose=verbose, fix_systematics=fix_systematics, param_type=param_type)
        
        # Convert MLE results to consistent precision
        model_MLE_param_fixed_hp = to_dict(model_MLE_param_fixed)
        model_MLE_hp = to_dict(model_MLE)
        
        # Compute NLLs with careful precision handling
        nll_fixed = nll(data_hp, model_MLE_param_fixed_hp)
        nll_free = nll(data_hp, model_MLE_hp)
        
        # Return the difference with careful precision handling
        result = DEFAULT_FLOAT_TYPE(nll_fixed) - DEFAULT_FLOAT_TYPE(nll_free)
        
        if verbose:
                print(f"Data: {data_hp}")
                print(f"Model (fixed {param_type}): {model_MLE_param_fixed_hp}")
                print(f"Model (free {param_type}): {model_MLE_hp}")
                print(f"NLL (fixed {param_type}): {nll_fixed}")
                print(f"NLL (free {param_type}): {nll_free}")
                print(f"NLL Ratio Result: {result}")

        # Ensure result is non-negative (likelihood ratio should be >= 0)
        return DEFAULT_FLOAT_TYPE(max(result, 0.0))

def generate_toy_data(model, sample_size=1000, fix_systematics=False, param_type="r_Ge77m"):
        def sample(input_data):
                # Ensure all inputs use highest precision
                eps_sel = DEFAULT_FLOAT_TYPE(input_data["eps_sel"])
                unc_sel = DEFAULT_FLOAT_TYPE(input_data["unc_sel"])
                r_b = DEFAULT_FLOAT_TYPE(input_data["r_b"])
                exp = DEFAULT_FLOAT_TYPE(model["exp"])
                ratio = DEFAULT_FLOAT_TYPE(input_data["ratio"])
                
                # Handle scaling factor parameters
                scaling_factor = DEFAULT_FLOAT_TYPE(input_data["scaling_factor"])
                unc_scaling = DEFAULT_FLOAT_TYPE(input_data["unc_scaling"])
                
                if fix_systematics:
                        # Use nominal values without fluctuations
                        true_eps_sel = eps_sel
                        observed_eps_sel = eps_sel
                        true_scaling_factor = scaling_factor
                        observed_scaling_factor = scaling_factor
                else:
                        # Include systematic uncertainties
                        true_eps_sel = DEFAULT_FLOAT_TYPE(np.random.normal(eps_sel, unc_sel))
                        observed_eps_sel = DEFAULT_FLOAT_TYPE(np.random.normal(true_eps_sel, unc_sel))
                        true_scaling_factor = DEFAULT_FLOAT_TYPE(np.random.normal(scaling_factor, unc_scaling))
                        observed_scaling_factor = DEFAULT_FLOAT_TYPE(np.random.normal(true_scaling_factor, unc_scaling))

                if param_type == "r_Ge77m":
                      r_Ge77m = DEFAULT_FLOAT_TYPE(input_data["r_Ge77m"])
                else:  # param_type == "r_Ge77_and_Ge77m"
                      r_Ge77_and_Ge77m = DEFAULT_FLOAT_TYPE(input_data["r_Ge77_and_Ge77m"])
                      r_Ge77m = r_Ge77_and_Ge77m * scaling_factor

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
                        "eps_sel": observed_eps_sel,
                        "unc_sel": unc_sel,
                        "scaling_factor": observed_scaling_factor,
                        "unc_scaling": unc_scaling
                }
                
                return output

        return ak.Array([sample(model) for _ in range(sample_size)])

def generate_model(template_model, param_value, param_type="r_Ge77m"):
        output = template_model.copy()
        output["param_type"] = param_type
        
        if param_type == "r_Ge77m":
                output["r_Ge77m"] = param_value
                # Calculate derived r_Ge77_and_Ge77m
                scaling_factor = output["scaling_factor"]
                output["r_Ge77_and_Ge77m"] = param_value / scaling_factor
        else:  # param_type == "r_Ge77_and_Ge77m"
                output["r_Ge77_and_Ge77m"] = param_value
                # Calculate derived r_Ge77m
                scaling_factor = output["scaling_factor"]
                output["r_Ge77m"] = param_value * scaling_factor
        
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
    model, sample_size, fix_systematics, param_type = args
    return generate_toy_data(model, sample_size=sample_size, fix_systematics=fix_systematics, param_type=param_type)

def compute_nll_for_toy_batch(args):
    """Helper function to compute NLL for a batch of toy data"""
    toy_batch, model, fix_systematics, param_type = args
    return [nllr(data_test, model, fix_systematics=fix_systematics, param_type=param_type) for data_test in toy_batch]


def compute_nll_h0_for_model(args):
    """Helper function to compute NLL H0 for a specific model"""
    toy_h0, model, fix_systematics, param_type = args
    return [nllr(toy_h0_entry, model, fix_systematics=fix_systematics, param_type=param_type) for toy_h0_entry in toy_h0]


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
        range_points=40,
        fix_systematics=False,
        param_type="r_Ge77m",
        scaling_factor=0.5,
        unc_scaling=0.1
        ):

        if n_threads is None:
                n_threads = mp.cpu_count()

        with open(rc_rates_path, 'r') as f:
                rc_rates = yaml.safe_load(f)
        with open(sel_eff_path, 'r') as f:
                sel_eff = yaml.safe_load(f)
        selected_candidates = lh5.read_as("skm",selected_candidates_path,"ak")
        with open(exposure_path, 'r') as f:
                exposure = yaml.safe_load(f)

        n_obs = len(selected_candidates)
        n_b_ext = rc_rates["selection"]["mult_lar_psd"]["cts"]
        ratio = rc_rates["scaling"]
        exp = exposure["total_exposure"]

        selection_effs = [sel_eff["eps_time"], sel_eff["eps_prompt"], sel_eff["eps_delayed"]["energy"], sel_eff["eps_delayed"]["selection"]["mult_lar_psd"]]

        # Calculate with highest precision
        eps_values = [DEFAULT_FLOAT_TYPE(eps["val"]) for eps in selection_effs]
        eps_sel = DEFAULT_FLOAT_TYPE(np.prod(eps_values))
        print("Selection efficiencies:", eps_values)
        
        # Uncertainty propagation with highest precision
        unc_terms = []
        unc_individual = []
        for i in range(len(selection_effs)):
                unc_i = DEFAULT_FLOAT_TYPE(selection_effs[i]["unc"])
                unc_individual.append(unc_i)
                prod_others = DEFAULT_FLOAT_TYPE(np.prod([eps_values[j] for j in range(len(eps_values)) if j != i]))
                unc_terms.append((unc_i * prod_others) ** 2)
        print("Uncertainty individual terms:", unc_individual)
        print("Uncertainty terms:", unc_terms)
        
        unc_sel = DEFAULT_FLOAT_TYPE(np.sqrt(np.sum(unc_terms)))
        print("Selection efficiency:", eps_sel, "Uncertainty:", unc_sel)

        data = {
                "n_sb": DEFAULT_FLOAT_TYPE(n_obs),
                "n_b_ext": DEFAULT_FLOAT_TYPE(n_b_ext),
                "ratio": DEFAULT_FLOAT_TYPE(ratio),
                "exp": DEFAULT_FLOAT_TYPE(exp),
                "eps_sel": eps_sel,
                "unc_sel": unc_sel,
                "scaling_factor": DEFAULT_FLOAT_TYPE(scaling_factor),
                "unc_scaling": DEFAULT_FLOAT_TYPE(unc_scaling),
        }
        print("Data: ", data)
        print(f"Parameter type: {param_type}")

        model_MLE = estimate_MLE(data, fix_systematics=fix_systematics, param_type=param_type)
        print("Estimated MLE parameters:", model_MLE)

        para_range = np.linspace(0, range_max, range_points, dtype=DEFAULT_FLOAT_TYPE)
        bin_width = para_range[1] - para_range[0]
        # Use appropriate parameter for minimum value calculation
        if param_type == "r_Ge77m":
                min_val = np.min([model_MLE["r_Ge77m"], 0.01 * bin_width])
        else:
                min_val = np.min([model_MLE["r_Ge77_and_Ge77m"], 0.01 * bin_width])
        para_range[para_range<min_val] = min_val  # Avoid numerical issues with very small values

        model_range = ak.Array([generate_model(model_MLE, x, param_type) for x in para_range])

        print("Generating toy data...")
        with ProcessPoolExecutor(max_workers=n_threads) as executor:
                # Prepare arguments for parallel processing
                args_list = [(model_range[i], sample_size, fix_systematics, param_type) for i in range(len(model_range))]
                
                # Use list comprehension with executor.map for parallel processing
                toy_range = list(tqdm(executor.map(generate_toy_mcs_batch, args_list), 
                                    total=len(args_list), desc="toy generation"))
        

        print("Computing NLL range using multithreading...")
        with ProcessPoolExecutor(max_workers=n_threads) as executor:
                # Prepare arguments for parallel processing
                args_list = [(toy_range[i], model_range[i], fix_systematics, param_type) for i in range(len(toy_range))]
                
                # Use list comprehension with executor.map for parallel processing
                nll_range = list(tqdm(executor.map(compute_nll_for_toy_batch, args_list), 
                                    total=len(args_list), desc="NLL calculations"))
        
        nll_range = np.array(nll_range)
        data_nll = np.array([nllr(data, model_range[i], verbose=False, fix_systematics=fix_systematics, param_type=param_type) for i in range(len(toy_range))])
        p_value_range = np.array([np.count_nonzero(nll_range[i] >= data_nll[i])/len(nll_range[i]) for i in range(len(toy_range))])

        model_H0 = model_MLE.copy()
        # Set null hypothesis: primary parameter = 0
        if param_type == "r_Ge77m":
                model_H0["r_Ge77m"] = 0
                model_H0["r_Ge77_and_Ge77m"] = 0
        else:  # param_type == "r_Ge77_and_Ge77m"
                model_H0["r_Ge77_and_Ge77m"] = 0
                model_H0["r_Ge77m"] = 0

        toy_H0 = generate_toy_data(model_H0, sample_size=sample_size, fix_systematics=fix_systematics, param_type=param_type)

        print("Computing NLL H0 range using multithreading...")
        with ProcessPoolExecutor(max_workers=n_threads) as executor:
                args_list = [(toy_H0, model_range[i], fix_systematics, param_type) for i in range(len(toy_range))]
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

        import os
        
        # Determine output path - if it ends with .json, replace with .lh5, otherwise assume it's already .lh5
        if output_path.endswith('.json'):
                lh5_path = output_path.replace(".json", ".lh5")
        else:
                lh5_path = output_path
        
        # Remove existing file to avoid shape conflicts
        if os.path.exists(lh5_path):
                os.remove(lh5_path)
        
        # Ensure arrays have consistent dtypes for storage (convert from longdouble to float64)
        nll_range = np.array(nll_range, dtype=STORAGE_FLOAT_TYPE)
        data_nll = np.array(data_nll, dtype=STORAGE_FLOAT_TYPE)
        nll_H0_range = np.array(nll_H0_range, dtype=STORAGE_FLOAT_TYPE)
        para_range = np.array(para_range, dtype=STORAGE_FLOAT_TYPE)
        
        # Convert p-value arrays to consistent dtypes
        p_value_range = np.array(p_value_range, dtype=STORAGE_FLOAT_TYPE)
        median_H0 = np.array(median_H0, dtype=STORAGE_FLOAT_TYPE)
        bound_68_lower = np.array(bound_68_lower, dtype=STORAGE_FLOAT_TYPE)
        bound_68_higher = np.array(bound_68_higher, dtype=STORAGE_FLOAT_TYPE)
        bound_90_lower = np.array(bound_90_lower, dtype=STORAGE_FLOAT_TYPE)
        bound_90_higher = np.array(bound_90_higher, dtype=STORAGE_FLOAT_TYPE)
        
        # Create comprehensive LH5 output with all requested keys
        output_lh5 = {
                # MLE parameters - flatten the nested MLE dictionary
                "MLE_eps_sel": types.Array(np.array([STORAGE_FLOAT_TYPE(model_MLE["eps_sel"])], dtype=STORAGE_FLOAT_TYPE)),
                "MLE_exp": types.Array(np.array([STORAGE_FLOAT_TYPE(model_MLE["exp"])], dtype=STORAGE_FLOAT_TYPE)),
                "MLE_r_Ge77_and_Ge77m": types.Array(np.array([STORAGE_FLOAT_TYPE(model_MLE["r_Ge77_and_Ge77m"])], dtype=STORAGE_FLOAT_TYPE)),
                "MLE_r_Ge77m": types.Array(np.array([STORAGE_FLOAT_TYPE(model_MLE["r_Ge77m"])], dtype=STORAGE_FLOAT_TYPE)),
                "MLE_r_b": types.Array(np.array([STORAGE_FLOAT_TYPE(model_MLE["r_b"])], dtype=STORAGE_FLOAT_TYPE)),
                "MLE_ratio": types.Array(np.array([STORAGE_FLOAT_TYPE(model_MLE["ratio"])], dtype=STORAGE_FLOAT_TYPE)),
                "MLE_total_rate_scaling": types.Array(np.array([STORAGE_FLOAT_TYPE(model_MLE["scaling_factor"])], dtype=STORAGE_FLOAT_TYPE)),
                "MLE_unc_sel": types.Array(np.array([STORAGE_FLOAT_TYPE(model_MLE["unc_sel"])], dtype=STORAGE_FLOAT_TYPE)),
                "MLE_unc_total_rate_scaling": types.Array(np.array([STORAGE_FLOAT_TYPE(model_MLE["unc_scaling"])], dtype=STORAGE_FLOAT_TYPE)),
                
                # NLL data
                "data_nll": types.Array(data_nll),
                "nll_H0_range": types.Array(nll_H0_range),
                "nll_range": types.Array(nll_range),
                
                # P-value data
                "p_value_exp": types.Array(median_H0),
                "p_value_m_exp": types.Array(bound_68_lower),
                "p_value_m_exp_90": types.Array(bound_90_lower),
                "p_value_obs": types.Array(p_value_range),
                "p_value_p_exp": types.Array(bound_68_higher),
                "p_value_p_exp_90": types.Array(bound_90_higher),
                
                # Parameter range
                "para_range": types.Array(para_range)
        }
        
        print("LH5 output shapes:")
        print(f"MLE_eps_sel: {output_lh5['MLE_eps_sel'].nda.shape}")
        print(f"MLE_r_Ge77m: {output_lh5['MLE_r_Ge77m'].nda.shape}")
        print(f"MLE_r_Ge77_and_Ge77m: {output_lh5['MLE_r_Ge77_and_Ge77m'].nda.shape}")
        print(f"nll_range: {nll_range.shape}")
        print(f"data_nll: {data_nll.shape}")  
        print(f"nll_H0_range: {nll_H0_range.shape}")
        print(f"para_range: {para_range.shape}")
        print(f"p_value_obs: {p_value_range.shape}")
        print(f"p_value_exp: {median_H0.shape}")

        lh5.write(types.Struct(output_lh5), name="freq_calc_results", lh5_file=lh5_path)


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
                if isinstance(value, (int, float, np.number)):
                        result[key] = DEFAULT_FLOAT_TYPE(value)
                else:
                        result[key] = value
    elif hasattr(obj, 'fields'):
        # It's an awkward array
        for field in obj.fields:
                if isinstance(obj[field], (int, float, np.number)):
                        result[field] = DEFAULT_FLOAT_TYPE(obj[field])
                else:
                        result[field] = obj[field]

    else:
        # Try to iterate over it as if it's a mapping
        try:
            for key in obj:
                if isinstance(obj[key], (int, float, np.number)):
                        result[key] = DEFAULT_FLOAT_TYPE(obj[key])
                else:
                        result[key] = obj[key]
        except:
            raise ValueError(f"Cannot convert object of type {type(obj)} to dictionary")
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate MLE and confidence intervals for a given dataset.")
    parser.add_argument("output_path", type=str, help="Path to the output fi    le where results will be saved.")
    parser.add_argument("rc_rates_path", type=str, help="Path to the YAML file containing RC rates.")
    parser.add_argument("sel_eff_path", type=str, help="Path to the YAML file containing selection efficiency.")
    parser.add_argument("selected_candidates_path", type=str, help="Path to the LH5 file containing selected candidates.")
    parser.add_argument("exposure_path", type=str, help="Path to the YAML file containing exposure data.")
    parser.add_argument("--n_threads", type=int, default=None, help="Number of threads to use for parallel processing (default: all available cores).")
    parser.add_argument("--sample_size", type=int, default=10000, help="Number of samples for each model (default: 1000).")
    parser.add_argument("--range_max", type=float, default=2.0, help="Maximum range for the parameter (default: 2.0).")
    parser.add_argument("--range_points", type=int, default=40, help="Number of points in the parameter range (default: 40).")
    parser.add_argument("--fix_systematics", action="store_true", help="Fix systematic uncertainties (nuisance parameters) to their nominal values.")
    parser.add_argument("--param_type", type=str, default="r_Ge77m", choices=["r_Ge77m", "r_Ge77_and_Ge77m"], help="Parameter to optimize: 'r_Ge77m' or 'r_Ge77_and_Ge77m' (default: r_Ge77m).")
    parser.add_argument("--scaling_factor", type=float, default=0.5, help="Central value of scaling factor relating r_Ge77m and r_Ge77_and_Ge77m (default: 0.5).")
    parser.add_argument("--unc_scaling", type=float, default=0.1, help="Uncertainty on scaling factor (default: 0.1).")

    args = parser.parse_args()

    print("Starting MLE and CI calculation with the following parameters:")
    print(f"Output Path: {args.output_path}")
    print(f"RC Rates Path: {args.rc_rates_path}")
    print(f"Selection Efficiency Path: {args.sel_eff_path}")
    print(f"Selected Candidates Path: {args.selected_candidates_path}")
    print(f"Exposure Path: {args.exposure_path}")
    print(f"Number of Threads: {args.n_threads}")
    print(f"Sample Size: {args.sample_size}")
    print(f"Fix Systematics: {args.fix_systematics}")
    print(f"Parameter Type: {args.param_type}")
    print(f"Scaling Factor: {args.scaling_factor} ± {args.unc_scaling}")

    calculate_mle_and_ci(
        args.output_path,
        args.rc_rates_path,
        args.sel_eff_path,
        args.selected_candidates_path,
        args.exposure_path,
        n_threads=args.n_threads,
        sample_size=args.sample_size,
        range_max=args.range_max,
        range_points=args.range_points,
        fix_systematics=args.fix_systematics,
        param_type=args.param_type,
        scaling_factor=args.scaling_factor,
        unc_scaling=args.unc_scaling
    )