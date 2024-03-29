% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mafCategoryAssign.R
\name{mafCategoryAssign}
\alias{mafCategoryAssign}
\title{Assign categories for maf and coverage}
\usage{
mafCategoryAssign(
  maf = NULL,
  coverage = NULL,
  ref_genome = NULL,
  category_num = NULL,
  output_file = TRUE,
  quiet = FALSE
)
}
\arguments{
\item{maf}{Mutation Annotation Format, can be whether an R data frame or the path of a txt file.}

\item{coverage}{Coverage file, can be whether an R data frame or the path of a txt file.}

\item{ref_genome}{Reference chromosome, either a folder name or a installed BSgenome package.
Default NULL, tries to auto-detect from installed genomes.}

\item{category_num}{Number of mutation categories, default 4. If the number is 0, categories should exist.}

\item{output_file}{Determine whether to export the categories, categorised maf and coverage as txt files.}

\item{quiet}{Whether to show notes during processing.}
}
\value{
Categorised maf and coverage.
}
\description{
The mutation categories generated from 192 basic types to designated number of types.
The reference chromosome is necessary if the number of categories is set more than one.
The chromosome can either be name of a folder containing txt files or name of the installed BSgenome package.
Example: BSgenome.Hsapiens.UCSC.hg19.
The preprocess of maf and coverage is integrated in this function.
}
\examples{
laml_maf <- system.file("extdata", "tcga_laml.maf", package = "DriverGenePathway")
coverage <- system.file("extdata", "coverage.rda", package = "DriverGenePathway")
load(coverage)
Categ_result <- mafCategoryAssign(laml_maf, coverage, category_num = 1, output_file = TRUE,
quiet = FALSE)
}
