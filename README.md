# Re-evaluating *A Novel Sparse Compositional Technique...*

This is a reanalaysis of Martino et al.'s 2019 manuscript, ["A Novel Sparse Compositional Technique Reveals Microbial Perturbations"](https://journals.asm.org/doi/10.1128/msystems.00016-19) published *mSystems*. The benchmarking code, primarily
written in Python, was obtained from [its GitHub repository](https://github.com/cameronmartino/deicode-benchmarking). I have a copy of the repository as originally posted on January 28, 2019 in the `deicode-benchmarking` directory.


## Running code

Make sure conda and mamba are installed. Installing TinyTex may take some fenagling and on Mac OSX seems to need to be installed in the home directory Library (i.e., `~/Library/TinyTex`) using `install.packages`. See [this issue](https://github.com/rstudio/tinytex/issues/24) for more clues on installing `tinytex`.

```bash
#mamba config --set channel_priority strict
mamba env create -f workflow/envs/XXXXXXX.yml 
conda activate XXXXXXX
```

Use Snakemake (installed in `XXXXXXX` environment) to build project

```bash
snakemake --use-conda --conda-frontend mamba -c1 write_paper
```


## Conda environments

* `workflow/env/martino.yml` - primary environment
* `workflow/env/rpca.yml` - environment for running RPCA-based analysis. Has a number of packages that were required for running different parts of the scripts from the `deicode_benchmarking` repository.


## Code modifications

Python/Numpy was pitching a fit about the use of `np.int` in these files  

  * `~/opt/miniconda3/envs/rclr/lib/python3.11/site-packages/deicode/_optspace.py`  
  * `~/opt/miniconda3/envs/rclr/lib/python3.11/site-packages/gemelli/optspace.py`

I replaced `np.int` with `int` and everything seemed to work fine
