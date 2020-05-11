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

Raw count table from fastq files.

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
install.packages("devtools")
devtools::install_github("benjaminweger/dryR")
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

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
dryList[["counts"]]# raw counts
dryList[["cook"]]  # cook's distance
```

