# threshold_paper_code
R code for plots and statistics used in threshold selection paper (Varty, Tawn, Atkinson and Bierman). 

## Dependencies 
To run this code several R packages are required. They may be installed by running the code below.
    
```
    required_pkgs <- c(
        "dplyr",
        "jsonlite",
        "lubridate",
        "matrixStats",
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

### `00_src/`
This folder contains the source code for R functions that are used within the later directories. 
    
`catalogue_creation/` functions for downloading the Groningen data. 



### ```01_introduction/```
There is currently no code associated with to the introduction. 

### ```02_motivation_and_model/```
This folder contains the code required to produce:

Figures of the Groningen earthquake catalogue on natural and index time scales. These are created in ```mativating_data.R``` and stored in ```output/groningen_catalogue.pdf```.

### ```03_variable_threshold_benefits``` 
This directory contains the code required for the simulation study that compares parameter and return level estimation using catalogues resulting from a conservative, stepped or extended threshold. 

By running `motivating_example.R` the following will be created in `output/plots/`: 

    - A example catalogue plot, showing which events are include under a conservative and stepped threshold; 
    - The bootstrap MLEs for GPD parameters using each threshold type;
    - The MLE MSE decomposition using each threshold type;
    - The MSE reduction factor relative to the conservative threshold;
    - The estimated conditional return levels (above the conservative threshold) using each threshold.





## License 
[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)