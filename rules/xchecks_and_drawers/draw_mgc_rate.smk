rule draw_mgc_rate:
    """
    Once sum_mgc.lh5 is created, this rule draws the mgc rate plot.
    """

    input:
        sum_mgc="gen/mu_hpge_coinc/sum_mgc.lh5"
    output:
        "gen/figs/mgc_rate.pdf"
    run:
        
        import lgdo.lh5 as lh5
        import numpy as np
        import matplotlib.pyplot as plt
        import lgdo.types as types
        from ge77m_search_workflow import utils as ut
        # Read the data from the LH5 file
        tmp = lh5.read_as("mgc", str(input.sum_mgc), "ak")

        # narrow plot

        mask = tmp["geds"]["energy"] > 500

        timestamps = np.array(tmp["coinc"]["id"]["timestamp"][mask])
        plt.hist(timestamps, bins=np.linspace(np.min(timestamps), np.max(timestamps), 101), histtype='step', label='>500 keV', color='black')
        plt.xlabel("Timestamps")
        plt.ylabel("MGC events")
        plt.savefig(str(output), bbox_inches="tight")
        plt.close()

