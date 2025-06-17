rule draw_mu_diff_vs_energy:
    """
    Once skm_mgc.lh5 is created, this rule draws the mu_diff_vs_energy plot.
    """

    input:
        skm_mgc="gen/mu_hpge_coinc/skm_mgc.lh5"
    output:
        "gen/figs/mu_diff_vs_energy.pdf"
    run:
        
        import lgdo.lh5 as lh5
        import numpy as np
        import matplotlib.pyplot as plt
        import lgdo.types as types
        from ge77m_search_workflow import utils as ut
        # Read the data from the LH5 file
        tmp = lh5.read_as("mgc", str(input.skm_mgc), "ak")

        # narrow plot

        mask = tmp["mu"]["coinc_flags"]["mu"]

        plt.hist2d(np.array(tmp["coinc"]["mu_diff"][mask]), np.array(tmp["geds"]["energy"][mask]), 
           bins=[np.linspace(-5000,10000, 101), 10**np.linspace(1.1,4,101)], 
           cmin=1)
        plt.yscale("log")
        plt.axhline(500, color="red", linestyle="--")
        plt.axvline(-2000, color="red", linestyle="--")
        plt.axvline(5000, color="red", linestyle="--")
        plt.xlabel("HPGe - Muon triggers (ns)")
        plt.ylabel("Energy (keV)")

        plt.savefig(str(output), bbox_inches="tight")
        plt.close()

        # wide plot in mu diff

        plt.hist2d(np.array(tmp["coinc"]["mu_diff"][mask]), np.array(tmp["geds"]["energy"][mask]), 
           bins=[np.linspace(-5000,54000, 101), 10**np.linspace(1.1,4,101)], 
           cmin=1)
        plt.yscale("log")
        plt.axhline(500, color="red", linestyle="--")
        plt.axvline(-2000, color="red", linestyle="--")
        plt.axvline(5000, color="red", linestyle="--")
        plt.xlabel("HPGe - Muon triggers (ns)")
        plt.ylabel("Energy (keV)")
        plt.savefig(str(output).replace(".pdf", "_wide.pdf"), bbox_inches="tight")
        plt.close()