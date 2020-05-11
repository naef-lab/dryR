# dryR
`dryR` is an R package that provides the statistical framework to assess differential rhythmicity of temporal RNA-Seq datasets of two and more conditions.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

R

Raw count table from fastq files and how to get it. 

### Installing

To install dryR run the following code in R.
```
install.packages("devtools")
devtools::install_github("benjaminweger/dryR")
```

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

