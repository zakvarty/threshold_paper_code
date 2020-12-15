# threshold_paper_code
R code for plots and statistics used in threshold selection paper (Varty, Tawn, Atkinson and Bierman). 

## Dependencies 
To run this code several R packages are required. They may be installed by running the code below.
    
```
    required_pkgs <- c(
        "dplyr",
        "ggplot2"
        "jsonlite",
        "lubridate",
        "matrixStats",
        "purrr",
        "readr",
        "sp"
        )
        
    install.packages(required_pkgs)
```

## Structure 
The directory structure represents the major sections of the paper, with the core shared code stored in `src/`. Further details of how to run the code for each section are given below.  


### `00_data/`
This directory contains the Groningen field outline and code to scrape earthquake data as a JSON file. This is then reformatted to a CSV earthquake catalogue. 

`00_data/input/` contains the Groningen earthquake outline with and without a 1000 m buffer in CSV forms. 

`00_data/output/` contains full earthquake catalogues in JSON and CSV form. The naming convention is based on the catalogue download time: "YYYY-MM-DD_HH-MM-SS_cat.csv".

`00_data/make_catalogue.R`
This file scrapes the current Groningen catalogue as a JSON file and formats it as a CSV file, saving both forms (lightly adapted from code provided by Stijn) 

To download an updated catalogue, set `download_new_json = FALSE` on line 12 and amend the `JSONpath` in line 15 with the new download date and time.

### 00_src/
This folder contains the source code for R functions that are used within the later directories. 
    
`catalogue_creation/` functions for downloading the Groningen data. 

### 01_introduction/
There is currently no code associated with to the introduction. 

### 02_motivation_and_model/
This folder contains the code required to produce:

Figures of the Groningen earthquake catalogue on natural and index time scales. These are created in ```mativating_data.R``` and stored in ```output/groningen_catalogue.pdf```.

### 03_variable_threshold_benefits/ 
This directory contains the code required for the simulation study that compares parameter and return level estimation using catalogues resulting from a conservative, stepped or extended threshold. 

By running `motivating_example.R` the following will be created in `output/plots/`: 

    - A example catalogue plot, showing which events are include under a conservative and stepped threshold; 
    - The bootstrap MLEs for GPD parameters using each threshold type;
    - The MLE MSE decomposition using each threshold type;
    - The MSE reduction factor relative to the conservative threshold;
    - The estimated conditional return levels (above the conservative threshold) using each threshold.


### 04_threshold_selection_method/
This directory contains the code required to produce the modified QQ and PP plots used when introducing our proposed threshold selection metric. A catalogue is simulated with the same stepped  threshold and hard censoring during the simulation study of `03_variable_selection_benefits/`. 

1000 GPD magnitudes exceeding 1.05ML are simulated. Hard censoring is then applied with the first 500 events have a threshold of 1.65ML which reduces to 1.05ML for the final 500 events. This leaves 582 magnitudes, which are rounded to the neared 0.1ML. 

QQ and PP plots are constructed on Gaussian and exponential margins based on magnitudes exceeding 3 thresholds: one below both threshold levels, one between the threshold levels and one above both. These threshold values are set on line 86 to be `thresholds_vec =  c(0.5, 1.15, 1.85)`. The plots above 1.15ML and 1.85ML are currently used in the paper. 

### 05_threshold_selection simulation/  
This directory contains the code required to run the simulation studies in Section 5 of the paper. These simulations are computationally intensive and  have been written for implementation on a remote computing cluster.

The experiments are initially organised by threshold type: constant, stepped or change-point. Within these, options for hard or phased censoring are investigated for some threshold types. 

For each experiment two directories exist:
 `./directory_name/` is pushed to the computing cluster and selects thresholds for many simulated catalogues of the same structure (structure and selection method varies between experiments). `./directory_name_pulled/` is pulled from the computing cluster. It is in this directory that the plots are constructed by running `plots.R. `

Currently, only the flat and change-point simulations are included in the paper.

_NB: In the git repo these directories will appear to be the same since data files are not version controlled. However, in the local version these are different. The empty (pushable) version is kept to make re-running the experiments easier when revisions are required._ 

To be clear, the workflow when re-running each experiment is: 

    1. Check that `submit.sh` has correct limits on core use and will email on completion; 
    2. Upload the directory without '_pulled' to STORM and 'qsub submit.sh' ;
    3. When code completed, pull the storm directory to 'directory_name_pulled';
    4. Within the pulled directory, run 'plots.R' to create the plots to be included in the paper.


## License 
[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)