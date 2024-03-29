% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/backgroundMutationRate.R
\name{backgroundMutationRate}
\alias{backgroundMutationRate}
\title{Calculate background mutation rates}
\usage{
backgroundMutationRate(
  maf = NULL,
  coverage = NULL,
  covariate = NULL,
  original_mutation_rate = 1.2e-06,
  max_neighbors = 50,
  ref_genome = NULL,
  category_num = NULL,
  output_file = TRUE,
  quiet = FALSE
)
}
\arguments{
\item{maf}{Mutation Annotation Format, can be whether an R data frame or the path of a txt file.}

\item{coverage}{Coverage file, can be whether an R data frame or the path of a txt file.}

\item{covariate}{Covariate file, can be whether an R data frame or the path of a txt file.}

\item{original_mutation_rate}{The default background mutation rate is 1.2e-6, which is used to
screen a subset of genes through binomial test.}

\item{max_neighbors}{Max number of gene neighbors used to calculate the background mutation rate
of a specified gene.}

\item{ref_genome}{Reference chromosome, either a folder name or a installed BSgenome package.
Default NULL, tries to auto-detect from installed genomes.}

\item{category_num}{Number of mutation categories, default 4. If the number is 0, categories should exist.}

\item{output_file}{Determine whether to export the categories, categorised maf and coverage as txt files.}

\item{quiet}{Whether to show notes during processing.}
}
\value{
The output BMR_out is an input of the DriverGenes function.
}
\description{
This function uses the data output by the preprocessing function in order to get
background mutation rate. It uses the method of bagel. The result will be then
used in varified specific gene detecting methods.
}
\details{
This function first selects a subset of genes whose mutation rates are significantly
greater than original mutation rate according to binomial test, then calculates background mutation
rates for these genes.
}
\examples{
laml_maf <- system.file("extdata", "tcga_laml.maf", package = "DriverGenePathway")
coverage <- system.file("extdata", "coverage.rda", package = "DriverGenePathway")
covariate <- system.file("extdata", "gene.covariates.txt", package = "DriverGenePathway")
load(coverage)
bmr_result <- backgroundMutationRate(laml_maf, coverage, covariate, original_mutation_rate = 1.2e-6,
max_neighbors = 50, ref_genome = NULL, category_num = 1, output_file = TRUE, quiet = FALSE)
}
