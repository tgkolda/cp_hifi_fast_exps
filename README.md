# GITHUB Repo for CP-HIFI Fast Experiments
This repository contains the examples from Brust and Kolda (2026), "Fast and Accurate CP-HIFI Tensor Decompositions: Exploiting Kronecker Structure". 
It depends on the algorithm implementations available from the [CP-HIFI Code](https://github.com/tgkolda/cp_hifi_code)
and the [Tensor Toolbox for MATLAB](https://www.tensortoolbox.org) (Bader, Kolda et al., Version 3.8, December 7, 2025).

**Reproducibility**: 

* All results in this repository [commit `57cff5a`] were obtained using the [CP-HIFI Code](https://github.com/tgkolda/cp_hifi_code) [commit `a1f74a4`] 
   
## Setup

1. Download version 3.8 of [Tensor Toolbox for MATLAB](https://gitlab.com/tensors/tensor_toolbox/-/releases/v3.8). Unzip it and add its directory to your MATLAB path. Below follows the git instructions to obtain the code, but you can also just download the .zip file from the release page and add it to your MATLAB path.
    ```
    git clone https://gitlab.com/tensors/tensor_toolbox
    cd tensor_toolbox
    git checkout v3.8
    cd ..
    ``` 

2. Obtain the datasets: 

    - The Miranda tensor may be obtained from the [Miranda Turbulent Flow Dataset](https://gitlab.com/tensors/tensor_data_miranda_sim) (Ballard et al., retrieved Fall 2025)

        ```
        git clone https://gitlab.com/tensors/tensor_data_miranda_sim
        cd tensor_data_miranda_sim
        git checkout 5c2d8b7
        cd ..
        ```
  
    - The Vortex tensor may be obtained from [Tensor Data Vortex Shedding](https://gitlab.com/bwlarsen/tensor_data_vortex_shedding) (Larsen, retrieved Fall 2025).
  
        ```
        git clone https://gitlab.com/bwlarsen/tensor_data_vortex_shedding
        cd tensor_data_vortex_shedding
        git checkout a83b235
        cd ..
        ```

3. Clone the [CP-HIFI Code](https://github.com/tgkolda/cp_hifi_code), which contains the implementations of the CP-HIFI algorithms.

    ```
    git clone https://github.com/tgkolda/cp_hifi_code
    cd cp_hifi_code 
    git checkout a1f74a4
    cd ..
    ```

4. Clone this repository, which contains the code for the fast experiments.

    ```
    git clone https://github.com/tgkolda/cp_hifi_fast_exps
    cd cp_hifi_fast_exps
    git checkout 57cff5a
    ```


## Replicate Experiments

* Our experiments were performed on a linux machine with Ubuntu 22.04.5 LTS. 
  The system has 128 GB RAM, 13th Gen Intel Core i9-13900KS × 32 CPUs, and a NVIDIA GeForce RTX 4090 GPU. All computations are in 64-bit on Matlab R2023a.
  Computational times will likely vary due to machine and computational load differences, but relative errors should be unchanged on repeated runs. 

* Be sure that the paths to the Tensor Toolbox, CP-HIFI Code, and both datasets are in your path before running the experiments. Here `/path/to/` should be replaced with the actual path on your machine. 

    ``` matlab
    addpath('/path/to/tensor_toolbox');
    addpath('/path/to/cp_hifi_code');
    addpath('/path/to/tensor_data_miranda_sim');
    addpath('/path/to/tensor_data_vortex_shedding');
    ```

* The folders are organized as follows:
    * **`vortex_experiments`** contains the code for the vortex shedding experiments. The main drivers are `run_aligned_vortex.m` and `run_unaligned_vortex.m`. All the other MATLAB files in this folder are helper functions that are called by the main drivers.
    * **`miranda_experiments`** contains the code for the miranda turbulent flow experiments. The main drivers are `run_aligned_miranda.m` and `run_unaligned_miranda.m`. All the other MATLAB files in this folder are helper functions that are called by the main drivers.
    * **`results`** contains the results of the experiments
    * **`paper_figs`** contains code to generate datafiles for the paper, from the results in `results/`

* To re-run an experiment, navigate into the appropriate folder (`vortex_experiments` or `miranda_experiments`) and use one of two main drivers (`run_aligned_[dataset].m` or `run_unaligned_[dataset].m`). **Note:** By default, the extremely slow direct methods are disabled via  `hasdirect = 0` in the appropriate driver; you can enable them by setting `hasdirect = 1`. The experiments generate three types of data for a combination of `[dataset]_[type]_[solver]`. The value for `[dataset]` is either {vortex,miranda} and for `[type]` it is {aligned, unaligned}.
`[solver]` is among {pcg,direct,direct_decoupled,direct_nonsym}. The generated logs, and run data is stored as:
     * `results/[dataset]_[type]_log_DD-Mmm-YYYY-HH_MM_SS.txt` contains the log of the run
     * `results/[dataset]_[type]_[solver]_all.txt` contains the rank ($r$), runid (run identifier in {1,2,3}), seed (random seed), time (total runtime), rerr (relative error), it (total outer iterations) for all runs using the specified solver
     * `results/[dataset]_[type]_[solver]_best.txt` contains the data for the best run (lowest relative error) for each rank $r$. 
    

    - Vortex aligned experiments:

        ``` matlab
        cd vortex_experiments
        run_aligned_vortex
        cd ..
        ```
    
        This produces the following files:
        * `results/vortex_aligned_log_DD-Mmm-YYYY-HH_MM_SS.txt` 
        * `results/vortex_aligned_[solver]_all.txt` 
        * `results/vortex_aligned_[solver]_best.txt` 
    
    - Vortex unaligned experiments:

        ``` matlab
        cd vortex_experiments
        run_unaligned_vortex
        cd ..
        ```
        
        This produces the following files:
        * `results/vortex_unaligned_log_DD-Mmm-YYYY-HH_MM_SS.txt` 
        * `results/vortex_unaligned_[solver]_all.txt` 
        * `results/vortex_unaligned_[solver]_best.txt`

    - Miranda aligned experiments:

        ``` matlab
        cd miranda_experiments
        run_aligned_miranda
        cd ..
        ```

        This produces the following files:
        * `results/miranda_aligned_log_DD-Mmm-YYYY-HH_MM_SS.txt` 
        * `results/miranda_aligned_[solver]_all.txt` 
        * `results/miranda_aligned_[solver]_best.txt`

    - Miranda unaligned experiments:

        ``` matlab
        cd miranda_experiments
        run_unaligned_miranda
        cd ..
        ```

        This produces the following files:
        * `results/miranda_unaligned_log_DD-Mmm-YYYY-HH_MM_SS.txt` 
        * `results/miranda_unaligned_[solver]_all.txt` 
        * `results/miranda_unaligned_[solver]_best.txt`


## Creating Data for Figures from the Results Files

* The folder `results/` stores outputs from the various runs of `run_[type]_[dataset].m`. 
    * Log files: `[dataset]_[type]_log_[datestamp].txt` 
    * All runs of a particular solver (i.e., three repetitions with different seeds) are in files `[dataset]_[type]_[solver]_all.txt`. These contain the rank ($r$), runid (run identifier in {1,2,3}), seed (random seed), time (total runtime), rerr (relative error), it (total outer iterations) for all runs using the specified solver
    These files contain the rank and seed information to replicate each run.
    * Best runs of a particular solver are in `[dataset]_[type]_[solver]_best.txt`
    * Direct solvers (because they are slow to run) also store a run "trace". These are in `[dataset]_[type]_[solver]_trace.mat`. The `.mat` files are 10x2 cell arrays where the first entry   corresponds to the problem rank. The second entry is a 3x3 cell array, that stores `run_id`, `seed`, and `info` where `info` is a struct that contains the detailed traces for this run.
* The folder `pager_figs/` contains pre-processed data from `results/`. The files are created by running the scripts `genalltabs.m`.  It contains 
    * `[dataset]_errs_[type].txt`
    * `[dataset]_times_[type].txt` 

## Vortex Slice Images

* To generate the vortex slice images:        
    
    ``` matlab
    cd vortex_experiments
    plot_vortex_slices_a
    cd ..
    ```
        
    This produces the following files:
     * `paper_figs/vortex_original_a_exp.pdf` Original frontal slice no. 151 of the vortex tensor
     * `paper_figs/vortex_samp_a_exp.pdf` Original frontal slice no. 151 of the vortex tensor sampled at 50,000 points
     * `paper_figs/vortex_slices_pcg_a_exp.pdf` CP-HIFI frontal slice no. 151 of the vortex tensor computed using pcg 
     * `paper_figs/vortex_direct_nonsym_a_exp.pdf` CP-HIFI frontal slice no. 151 of the vortex tensor computed using direct nonsymmetric 
     