import lgdo.lh5 as lh5
import numpy as np
import awkward as ak
import yaml
import lgdo.types as types
import os

def calculate_rc_bkg_old(
            skm_file,
            exposure_file,
            output_rates,
            output_distr,
            delayed_e_range,
            on_dataset,
            off_dataset,
            psd_ms_type,
            draw_distribution=False
    ):

    if draw_distribution:
        fig_folder = str(os.path.dirname(exposure_file)).replace("exposure","figs")
    
    # load skm and exposure data
    skm = lh5.read_as("skm/",skm_file,"ak")
    with open(exposure_file,"r") as f:
        exp = yaml.safe_load(f)
    
    # select time in off 
    mask_off_data = (skm["dT"] > off_dataset[0]) & (skm["dT"] < off_dataset[1])
    skm_off = skm[mask_off_data]

    # generate masks for different selections
    mask_mult = skm_off["delayed"]["coinc"]["multiplicity"] == 1
    mask_lar = ~(skm_off["delayed"]["coinc"]["spm_coinc"])
    mask_psd = np.arange(len(skm_off))
    if psd_ms_type != "any":
        mask_psd = skm_off["delayed"]["geds"]["quality"]["psd_is_bb_like"]
    mask_energy = np.sum(np.array([(skm_off["delayed"]["geds"]["energy"] > r[0]) & (skm_off["delayed"]["geds"]["energy"] < r[1]) for r in delayed_e_range]).T,axis=-1) > 0

    # generate combinations of masks
    masks = {
        "all": np.ones(len(skm_off),dtype=bool),
        "mult": mask_mult,
        "mult_lar": mask_mult & mask_lar
    }
    
    if psd_ms_type != "any":
        masks["mult_lar_psd"] = mask_mult & mask_lar & mask_psd

    # calculate relative fraction of events in on dataset compared to off dataset by fitting the dT distribution in the off dataset and extrapolating
    def estimate_relative_fraction_on_vs_off(dataset):
        # get dT distribution of off dataset
        bins = np.linspace(off_dataset[0],off_dataset[1],201)
        bin_width = bins[1] - bins[0]
        hist_off, _ = np.histogram(dataset["dT"], bins)
        area_off = np.sum(hist_off)

        # fit distribution with a linear function
        coeffs = np.polyfit(bins[:-1], hist_off, 1)

        # estimate area (integral) in the on dataset
        area_on = (coeffs[1] * (on_dataset[1] - on_dataset[0]) + coeffs[0] * (on_dataset[1] - on_dataset[0])**2 / 2)/bin_width

        if draw_distribution:
            import matplotlib.pyplot as plt
            plt.hist(dataset["dT"], bins=bins, alpha=0.5, label="Off dataset")
            plt.axvline(on_dataset[1], color='r', linestyle='--')
            x = np.linspace(on_dataset[0], off_dataset[1], 100)
            y = np.polyval(coeffs, x)
            plt.plot(x,y,ls=":",color="black")
            plt.xlabel("dT")
            plt.ylabel("Counts")
            plt.savefig(f"{fig_folder}/dT_distr.pdf")
            plt.clf()

        return area_on / area_off if area_off > 0 else 0

    relative_fraction = estimate_relative_fraction_on_vs_off(skm_off)
    print(f"Relative fraction of on vs off dataset: {relative_fraction:.2f}")

    # calculate background rates
    rates = {}
    unc = {}
    counts = {}
    for key, mask in masks.items():
        print("total selection:", np.sum(mask_energy & mask))
        rates[key] = float(relative_fraction * np.sum(mask_energy & mask) / exp["total_exposure"])
        unc[key] = float(relative_fraction * np.sqrt(np.sum(mask_energy & mask)) / exp["total_exposure"])
        counts[key] = int(np.sum(mask_energy & mask))
        print(f"Rate for {key}: {rates[key]} +/- {unc[key]} cts/(kg yr)")
    


    
    # calculate pdfs
    pdfs = {"bins": np.linspace(0,3000,301)}
    bin_width = pdfs["bins"][1] - pdfs["bins"][0]
    for key, mask in masks.items():
        pdfs[key] = np.histogram(skm_off["delayed"]["geds"]["energy"][mask], bins=pdfs["bins"], weights=np.full(len(skm_off["delayed"]["geds"]["energy"][mask]),1/(exp["total_exposure"]* bin_width) ))[0]

    if draw_distribution:
        import matplotlib.pyplot as plt
        for key, pdf in pdfs.items():
            if key == "bins":
                continue
            plt.step(pdfs["bins"][:-1], pdf, label=key)
        plt.xlabel("Energy (keV)")
        plt.ylabel("Rate (cts/(keV kg yr))")
        plt.legend()
        plt.yscale("log")
        plt.savefig(f"{fig_folder}/energy_distr.pdf")
        plt.clf()


        for key, pdf in pdfs.items():
            if key == "bins":
                continue
            plt.step(pdfs["bins"][:-1], pdf, label=key)
        plt.xlabel("Energy (keV)")
        plt.ylabel("Rate (cts/(keV kg yr))")
        plt.legend()
        plt.yscale("log")
        plt.xlim(1500,2600)
        plt.savefig(f"{fig_folder}/energy_distr_zoom.pdf")
        plt.clf()

    # estimate total rc event rate for difference lower and upper delayed energy ranges
    e_low = np.linspace(1500, 2700, (2700 - 1500) // 100 + 1)
    e_high = np.linspace(1500, 2700, (2700 - 1500) // 100 + 1)
    x, y = np.meshgrid(e_low, e_high)
    scans = {}
    for key, mask in masks.items():
        z = np.zeros_like(x)
        for i in range(len(e_low)):
            for j in range(len(e_high)):
                if e_low[i] < e_high[j]:
                    mask_energy = (skm_off["delayed"]["geds"]["energy"] > e_low[i]) & (skm_off["delayed"]["geds"]["energy"] < e_high[j])
                    z[i, j] = float(relative_fraction * np.sum(mask_energy & mask) / exp["total_exposure"])
        scans[key] = z.tolist()
    scans["x"] = x.tolist()
    scans["y"] = y.tolist()

    output_data = {
        "selection": {key: {"val": rates[key], "unc": unc[key], "cts": counts[key]} for key in rates.keys()},
        "scans": scans,
        "scaling": float(relative_fraction),
    }


    # write out rates as yaml
    with open(output_rates, "w") as f:
        yaml.dump(output_data, f)

    # write out pdfs as lh5 file
    pdfs = {
        key: types.Array(value) for key, value in pdfs.items()
    }
    pdfs = types.Struct(pdfs)

    lh5.write(pdfs, name="pdfs", lh5_file=output_distr) # (output_lh5, name="mgc", lh5_file=str(output))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Calculate background rates from skm and exposure data.")
    parser.add_argument("skm_file", type=str, help="Path to the skm file.")
    parser.add_argument("exposure_file", type=str, help="Path to the exposure file.")
    parser.add_argument("output_rates", type=str, help="Output file for rates.")
    parser.add_argument("output_distr", type=str, help="Output file for distribution.")
    parser.add_argument("--delayed_e_range", default=[[1900.0, 2600.0]], help="Energy range for delayed events.")
    parser.add_argument("--on_dataset", nargs=2, type=float, default=[0.0, 387.0], help="Time range for on dataset.")
    parser.add_argument("--off_dataset", nargs=2, type=float, default=[387.0, 77473.0], help="Time range for off dataset.")
    parser.add_argument("--psd_ms_type", type=str, choices=["any", "low_aoe-coax_rt"], default="low_aoe-coax_rt", help="PSD multiplicity selection type.")
    parser.add_argument("--draw_distribution", action="store_true", help="Draw the dT distribution for the off dataset.")

    args = parser.parse_args()

    print(f"Calculating background rates from {args.skm_file} and {args.exposure_file} with delayed energy range {args.delayed_e_range}, on dataset {args.on_dataset}, off dataset {args.off_dataset}, psd_ms_type {args.psd_ms_type}, draw_distribution={args.draw_distribution}")

    calculate_rc_bkg(
        args.skm_file,
        args.exposure_file,
        args.output_rates,
        args.output_distr,
        args.delayed_e_range,
        args.on_dataset,
        args.off_dataset,
        args.psd_ms_type,
        args.draw_distribution
    )
