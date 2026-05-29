## Steps to Produce the Empirical Results

1. Ensure that all files are located in the same folder (e.g., `Empirical part`).
2. Open `Main File.R` in RStudio and specify the working directory by replacing the placeholder path in the line marked `"input the path to your 'Empirical part' folder here"`.
3. Run `Main File.R`. This script coordinates all required files and generates the empirical results.

A subfolder named `Results` will be created automatically. It contains two subfolders:

* `pvalue and covariates`
* `tables and figures`

The `pvalue and covariates` folder contains the calculated p-values, estimated alphas, and calculated covariates for all considered factor models. These files serve as inputs for the $fwer^+$ procedure.

The `tables and figures` folder contains all main empirical results reported in the paper, except for the communality heatmap. Tables are saved in HTML format, while figures are saved in PNG format.

After `Main File.R` has finished running, open and run `Communality Plot.R` to generate the communality heatmap.
