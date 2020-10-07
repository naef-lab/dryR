# dryR
`dryR` (Differential RhythmicitY analysis in R) is an R package that provides the statistical framework to assess differential rhythmicity of a time series of RNA-Seq data with two and more conditions.

## Getting Started

These instructions will allow you to get `dryR` running on your machine. 

### Prerequisites
You need to install R (see https://www.r-project.org/).

`dryR` accepts count data typically produced by RNA-Seq with a time dimension and several conditions/groups. The input count data should contain only integer values and be organized in a matrix with rows indicating a specific gene and the column refering to a sample. 

The pipelines that can be used to produce these count data from RNA-Seq data (FASTQ files) are described in more detail here:
http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

A user-friendly web application to generate a count table from FASTQ files is provided here: https://amp.pharm.mssm.edu/biojupies/

### Installing

To install `dryR` run the following code in R.
```
install.packages("devtools")
devtools::install_github("naef-lab/dryR")
```
## Quick start
### Example dataset 
`dryR` comes with example data in form of a list called simData. The list contains count data simData[["countData"]], a vector with the different conditions/groups simData[["group"]], and a vector indicating Zeitgeber Time simData[["time"]]. The data was generated using simphony https://github.com/hugheylab/simphony.

### Running an example
```
require("dryR")

# prepare arguments
countData = simData[["countData"]]
group     = simData[["group"]]
time      = simData[["time"]]

# run the analysis for count data (e.g. RNA-Seq data)
dryList   = dryseq(countData,group,time)

# explore the results
dryList[["results"]]    # data frame summarizing results
dryList[["parameters"]] # coefficients: phase, amplitude and mean for each group
dryList[["ncounts"]]    # normalized counts
dryList[["counts"]]     # raw counts
dryList[["cook"]]       # cook's distance for outlier detection

#plot a feature of interest
dry_plot(dryList, "feature_113")

#generate a pdf with a global summary of all models
plot_models_rhythm(dryList, "./")
```


## Normally distributed data
To asses temporal variation of normally distributed measurements, we implemented the function `drylm` that can deal with gaussian noise using linear models. 

```
require("dryR")

# prepare arguments
data      = log(simData[["countData"]]+1)
group     = simData[["group"]]
time      = simData[["time"]]

# run the analysis with normally distributed data
dryList = drylm(data,group,time)

# explore the results
dryList[["results"]]    # data frame summarizing results
dryList[["parameters"]] # coefficients: phase, amplitude and mean for each group

#plot a feature of interest
dry_plot(dryList, "feature_013")

#generate a pdf with a global summary of all models
plot_models_rhythm(dryList, "./")
```

## DryR with a simple vector as input
You can run `drylm` with a simple vector that contains data from a time series of multiple groups. 

```
require("dryR")

# define time and group for each sample
time        = c(1:48,1:48)   # Zeitgeber time or Circadian time in h for each sample
group       = c(rep("KO",48), rep("WT",48))

# run the analysis with normally distributed data
dryList = drylm(vector_example,group,time)

# explore the results
dryList[["results"]]    # data frame summarizing results. Row 2 is a copy of row 1.
dryList[["parameters"]] # coefficients: phase, amplitude and mean for each group. Row 2 is a copy of row 1.

#plot the result of the selected model
dry_plot(dryList, "X1")
```

## Help
A documentation using `?dryseq` pr `?drylm` function is available. 
