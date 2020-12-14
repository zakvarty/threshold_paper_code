# threshold_paper_code
R code for plots and statistics used in threshold selection paper (Varty, Tawn, Atkinson and Bierman). 

## Structure 
The directory structure represents the major sections of the paper, with the core shared code stored in ```src/```. Further details of how to run the code for each section are given below.  


### ```00_data/```
This directory contains the raw Groningen earthquake data, a script to subset this to the time period (1995-01-01, 2020-01-01] and the final formatted earthquake catalogue. 

```00_data/input/``` contains the Grongingen earthquake outline with and without a 1000m buffer in csv forms. 

```00_data/output/``` contains full earthquake catalogues in .json and .csv form. The naming convention is based on the catalogue download time: "YYYY-MM-DD_HH-MM-SS_cat.csv".

``````00_data/make_catalogue.R```
This file scrapes the current Groningen catalogue as a .json file and formats it as a .csv file. (lightly adapted from code provided by Stijn) 

To download an updated catalogue, set ``download_new_json = FALSE````on line 12 and amend the JSON path in line 15 with the new download date and time.

### ```00_src/```
This folder contains the source code for R functions that are used within the later directories. 

### ```01_introduction/```
There is currently no code associated with to the introduction. 

### ```02_motivation_and_model/```
This folder contains the code required to produce: 

    - Figures of the Groningen earthquake catalogue on natural and index time scales. (Fig 1)

## License 
[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)