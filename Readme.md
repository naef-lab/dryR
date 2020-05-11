# dryR
`dryR` is an R package that provides the statistical framework to assess differential rhythmicity of temporal RNA-Seq datasets of two and more conditions.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

R uses raw counts 

dryR accepts a raw count dataset . Each row 
To get such a file you might want to use to map. or use https://amp.pharm.mssm.edu/biojupies/

### Installing

To install dryR run the following code in R.
```
install.packages("devtools")
devtools::install_github("benjaminweger/dryR")
```

## example dataset 
`dryR` comes with example data in form of a list called simData. The list contains raw count data simData[["countData"]], a vector with the different conditions/groups simData[["group"]], and a vector indicating Zeitgeber Time simData[["time"]]. The data has been generate using simphony https://github.com/hugheylab/simphony.

## Running an example
```
require("dryR")

# prepare arguments
countData = simData[["countData"]]
group     = simData[["group"]]
time      = simData[["time"]]

# run the analysis
dryList   = dryseq(countData,group,time)

# explore the results
dryList[["results"]]    # data frame summarizing results
dryList[["parameters"]] # coefficients: phase, amplitude and mean for each group
dryList[["ncounts"]]    # normalized counts
dryList[["counts"]]     # raw counts
dryList[["cook"]]       # cook's distance for outlier detection
```

## Help
A documentation using `?dryseq` is available. 
