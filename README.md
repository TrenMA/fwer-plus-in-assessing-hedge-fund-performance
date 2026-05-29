# Assessing Hedge Fund Performance with an Information-Based Multiple Test

This repository provides the replication code for the paper:

> Hsu, P.-H., Ma, T., Psaradellis, I., and Sermpinis, G. (2025). *Assessing Hedge Fund Performance with an Information-Based Multiple Test*.

**SSRN Paper:** https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5423519

The `fwer+` (family-wise error rate plus) procedure is an information-based multiple testing method for picking outperformers while controlling the family-wise error rate. This repository contains the code used to generate the main simulation and empirical results reported in the paper. This repository contains the code used to generate the main simulation and empirical results reported in the paper.

## Repository Structure

The repository is organised into two main folders:

* `Simulation part`
* `Empirical part`

### `Simulation part`

This folder contains the code used to generate the simulation results presented in the paper.

Detailed instructions are provided in the corresponding README file within the folder.

### `Empirical part`

This folder contains the code used to generate the empirical results reported in the paper.

Detailed instructions are provided in the corresponding README file within the folder.

## Requirements

The code is written in R and is designed to run in RStudio.

Before running the scripts, please ensure that all required R packages are installed.

## Simulated data

This repository includes the necessary simulated data files in both the `Simulation part` and `Empirical part` folders. These simulated data are provided solely to enable users to run the code smoothly; they are not intended to reproduce the results reported in the paper.

The original data used in the paper were obtained from commercial data providers and publicly available datasets maintained by other researchers. These data are not distributed through this repository. Consequently, the results generated using the simulated data will not match those reported in the paper.

Researchers wishing to reproduce the empirical results reported in the paper must obtain the original data from the relevant commercial data providers (see Data section of the paper) and replace the simulated data files accordingly.

## Citation

If you use the code from this repository, please cite:

> Hsu, P.-H., Ma, T., Psaradellis, I., and Sermpinis, G. (2025). *Assessing Hedge Fund Performance with an Information-Based Multiple Test*. Working Paper.

## Contact

Questions regarding the code should be directed to Tren Ma at [tren.ma@nottingham.ac.uk](mailto:tren.ma@nottingham.ac.uk).
