% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dryseq_function.R
\name{dryseq}
\alias{dryseq}
\title{Differential rhythmicity analysis for RNA-Seq datasets}
\usage{
dryseq(
  countData,
  group,
  time,
  period = 24,
  sample_name = colnames(countData),
  batch = rep("A", length(sample_name)),
  n.cores = round(parallel::detectCores() * 0.6, 0)
)
}
\arguments{
\item{countData}{matrix containing non-negative integers; each column represents a sample, each row represents a gene/transcript.}

\item{group}{vector containing the name of each sample.}

\item{time}{vector containing numeric values of the time for each sample.}

\item{period}{numeric value to indicate period length of the oscillation. Default: circadian data period of 24 h.}

\item{sample_name}{vector containing sample names. Default: colnames are sample names.}

\item{batch}{vector containing potential batch effects between samples. Default: no batch effect.}

\item{nthreads}{vector numeric value to indicate the threads for parallel computing .Default: 60 \% of detected cores.}
}
\value{
a list that contains the following data.frames: results (summary of results), parameters (rhythmic parameters), ncounts (normalized counts), counts (raw counts), cook (cook's distance)
}
\description{
This function performs a rhythmicity analysis based on generalized linear models with a subsequent models selection. The function accepts raw count data from a temporal RNA-Seq dataset of two or more groups. The function outputs parameters mean, phase and amplitude are for each group.
}
\details{
DryR assesses rhythmicity and mean differences of gene expression in RNA-Seq count data.
    As proposed (Love et al. 2014), a count Y of a gene in a sample s can be modeled as a negative binomial with a fitted mean μ_gs and a gene-specific dispersion parameter θ_g.
    \cr \cr  \emph{Y_gs}~NB(\emph{μ_gs},\emph{θ_g})\cr\cr
    The fitted mean is proportional to the quantity q of fragments that correspond to a gene in a sample scaled by a sample-specific scaling factor s_s (Love et al. 2014).
    This scaling factor depends on the sampling depth of each library and can be estimated using the median-of-ratios method of DESeq2 (Anders and Huber 2010).
    \cr \cr   \emph{μ_gs} = \emph{s_s q_gs}\cr \cr
    DryR estimates gene-specific distribution θg using empirical Bayes shrinkage described by Love et al. (Love et al. 2014).
    and variance is computed from the following relationship to the dispersion parameter θ:
    \cr \cr   Var(\emph{Y_gs}) = E[(\emph{Y_gs} + \emph{θ_g} \emph{μ^2_gs})]\cr \cr
    The fit uses a generalized linear model with a logarithmic link function. Sample specific size factor (λ_s) is defined as an offset. The full GLM is defined as follows:
    \cr \cr log(\emph{μ_gbcs}) = \emph{m_gb} + \emph{m_gc} + \emph{α_gc} cos(\emph{ω t(s)}) + \emph{β_gc} sin(\emph{ω t(s)}) + log(\emph{λ_t(s)})\cr \cr
     μ is the raw count for gene g, condition/group c and Zeitgeber/circadian time t. α and β are coefficients of the cosine and sine functions, respectively. m is a coefficient to describe a mean expression level.
     When necessary, a batch specific mean (m) can be given to the dryseq function to account for technical batch effects.
     A technical batch effect is not allowed to be confounding so the resulting model matrix is fully ranked.
     To select an optimal gene-specific model, dryseq first assesses rhythmicity across the different conditions. To this end, dryR defines different models across all groups.
     Models refined to have either zero (non-rhythmic pattern) or non-zero (rhythmic pattern) α and β coefficients for each analyzed group. Moreover, for some models the values of α and β can be also shared within any combination of all groups
     The coefficients α and β were used to calculate the phase (arctan(α/β)) and amplitude (log2-fold change peak-to-trough; 2sqrt(α^2+β^2) ) of a gene.
     Bayesian information criterion (BIC) based model selection was employed to account for model complexity using the following formula:
     \cr \cr   BIC_j = ln(n)k - 2ln(L̂_ĵ)\cr \cr
     L̂ is defined as the log-likelihood of the model j from the regression, n is the number of data points and k is the number of parameters.
     To assess the confidence of the selected model j we calculated the Schwarz weight (BICW):
     \cr \cr   BICW_j = e^(0.5ΔBIC_j)\ sum(e^0.5 ΔBIC_m), with ΔBIC_j - BIC_j - BIC_m*\cr \cr
     m* is the minimum BIC value in the entire model set. Dryseq consideres the BICW_j as the confidence level for model j. The model with the highest BICW is selected as the optimal model within the set of all defined models.
     In a second iteration step, dryR set the coefficient α and β to the values of the selected model in the first regression.
     dryseq then defined different models for the mean coefficient with differing or shared means between groups. Each model is solved using generalized linear regression and each gene was assigned to a preferred model based on the BICW as described above for the first iteration.
     The model selection is sensitive to outliers: dryseq provide a cook's distance for each gene. A fit for a gene with a maximum cook's distance of higher than 1 should be considered with care.
}
\examples{
countData = simData[["countData"]]
group = simData[["group"]]
time  = simData[["time"]]
dryList = dryseq(countData,group,time)
head(dryList[["results"]])    # data frame summarizing results
head(dryList[["parameters"]]) # coefficients: phase, amplitude and mean for each group
head(dryList[["ncounts"]])    # normalized counts
head(dryList[["counts"]])     # raw counts
head(dryList[["cook"]])       # cook's distance
}
\references{
Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology

Anders, S. and Huber, W. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology
}
