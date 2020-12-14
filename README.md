# threshold_paper_code
R code for plots and statistics used in threshold selection paper (Varty, Tawn, Atkinson and Bierman). 

## Structure 
The directory structure represents the major sections of the paper, with the core shared code stored in ```src/```. Further details of how to run the code for each section are given below.  

### ```00_src/```
This folder contains the source code for R functions that are used within the later directories. 

### ```00_data/```
This directory contains the raw Groningen earthquake data, a script to subset this to the time period (1995-01-01, 2020-01-01] and the final formatted earthquake catalogue. 

### ```01_introduction/```
There is currently no code associated with to the introduction. 

### ```02_motivation_and_model/```
This folder contains the code required to produce: 

    - Figures of the Groningen earthquake catalogue on natural and index time scales. (Fig 1)

## License 
[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)