## Note

To run the simulation, in addition to the provided simulated data, you need to collect:

- The Fung and Hsieh seven factors (see instructions: https://people.duke.edu/~dah7/HFRFData.htm)
- The monthly risk-free rate from the Kenneth R. French Data Library: https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html

To run the code seamlessly, these data must be stored as a dataframe in a file named `Factors.RDS`. The file should contain the following variables:

- `Date`: YYYYMM format  
- `EQ`: Equity market factor  
- `Size`: Size spread factor  
- `DGS10C`: Bond market factor  
- `CS`: Credit spread factor  
- `PTFSBD`: Return of PTFS bond lookback straddle  
- `PTFSFX`: Return of PTFS currency lookback straddle  
- `PTFSCOM`: Return of PTFS commodity lookback straddle  
- `RF`: Monthly risk-free rate (from the Kenneth R. French Data Library)

## Steps to Run the Files

1. Ensure that all files are located in the same folder (e.g., `Simulation part`).
2. Place the `Factors.RDS` file into this folder.
3. Run the file `Simulation Implementation.R`.
