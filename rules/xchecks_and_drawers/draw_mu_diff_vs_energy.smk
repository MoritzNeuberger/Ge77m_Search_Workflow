rule draw_mu_diff_vs_energy:
    """
    Once sum_mgc.lh5 is created, this rule draws the mu_diff_vs_energy plot.
    """

    input:
        sum_mgc="gen/mu_hpge_coinc/sum_mgc.lh5"
    output:
        "gen/figs/mu_diff_vs_energy.pdf"
    run:
        
        import lgdo.lh5 as lh5
        import numpy as np
        import matplotlib.pyplot as plt
        import lgdo.types as types
        from ge77m_search_workflow import utils as ut
        # Read the data from the LH5 file
        tmp = lh5.read_as("mgc", str(input.sum_mgc), "ak")

        # narrow plot

        mask = tmp["mu"]["coinc_flags"]["muon"]

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

        # wide plot 

        plt.hist(np.array(tmp["geds"]["energy"][mask]), 
           bins=np.linspace(0,10000,1001))
        plt.axvline(500, color="red", linestyle="--")
        plt.axvline(6500, color="red", linestyle="--")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Counts")
        plt.yscale("log")
        plt.savefig(str(output).replace(".pdf", "_just_energy.pdf"), bbox_inches="tight")
        plt.close()