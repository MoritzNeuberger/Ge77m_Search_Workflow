from setuptools import setup, find_packages

setup(
    name="ge77m_search_workflow",
    version="0.1.0",
    packages=find_packages(where="scripts"),
    package_dir={"": "scripts"},
    install_requires=[
        "numpy",
        "awkward",
        "pandas",
        "snakemake",
        "pygama",
        "legend-pydataobj",
        "pylegendmeta",
        "PyPDF2",
        "h5py",
        "numba_stats",
        "scipy",
        "iminuit",
        "tqdm",
        "matplotlib"
        # (empty, since we list deps in your envs/*.yaml)
    ],
)
