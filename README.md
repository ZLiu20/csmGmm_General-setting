# csmGmm_reproduce for the general multiple testing

Code to reproduce simulation results, figures, and real data analysis results for general composite null hypotheses.

## Overview

- Please first install the csmGmm package, for example through the command devtools::install_github("ryanrsun/csmGmm").

- Please navigate to the Data/ folder and either run the download_data.sh script or download the other necessary files manually using the other_data_locations.txt file. These other files are too large to be bundled with the rest of the code, so we keep them separately in Dropbox.

- All code provided is designed to be run on a high-performance computing cluster. We provide batch submission scripts for clusters using LSF.

- In general, one folder is provided for each figure or each table. For example, all the code needed to reproduce Figure 1 is placed inside the Fig1 folder.

- We use the 'here' package to refer to file locations within this project. When running this code, please make sure the working directory is set to somewhere within the csmGmm_General-setting folder to ensure that the here package works correctly.

- Each folder also contains example LSF batch submission scripts that end with the .lsf extension. These submission scripts are used to run the R scripts in parallel on a computing cluster, e.g. to run 200 simulations at once for Figure 3B. A few lines in these scripts will generally need to be slightly modified. For example, the job submission queues will have different names at different institutions. The number of jobs to run will be placed in a comment, e.g. "#Run 1-800" means that an array of jobs with IDs from 1 to 800 should be run.

- After making any necessary adjustments, all .lsf files in a folder should be run (e.g. with the command "bsub <run_Fig3B.lsf"), and then the user should wait until all jobs from the folder are finished running.

- Each folder also contains an R script that starts with the word "plot" such as "plot_Fig3.R." After all jobs in a folder have finished running, this script should be run to produce the final figure or table from the manuscript.


