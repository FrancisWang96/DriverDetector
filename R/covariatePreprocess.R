#' Preprocess covariate file
#'
#' This function guarantees that the covariate file is available for DriverGenePathway. Two kinds of covariate input
#' are feasible. It can be whether a txt file or an R data frame. The function returns the preprocessed
#' covariate. An txt file of the preprocessed covariate is optional to be exported.
#'
#' @param covariate Covariate file, can be whether an R data frame or the path of a txt file.
#' @param maf Mutation Annotation Format, can be whether an R data frame or the path of a txt file.
#' @param output_file Determine whether to export the preprocessed maf as a txt file.
#' @param quiet Whether to show notes during processing.
#' @return The preprocessed covariate file.
#' @examples
#' covariate <- system.file("extdata", "gene.covariates.txt", package = "DriverGenePathway")
#' laml_maf <- system.file("extdata", "tcga_laml.maf", package = "DriverGenePathway")
#' covariate <- covariatePreprocess(covariate, laml_maf)
#' @export

covariatePreprocess <- function(covariate, maf = NULL, output_file = TRUE, quiet = FALSE){

  if(!quiet){cat('Loading covariate ... ')}

  if(is.character(covariate)){
    covariate <- data.table::fread(covariate)
  }
  if(!quiet){cat('success!\n')}

  if(!quiet){cat('Check covariate elements ... ')}

  if(is.na(match('GENE',toupper(names(covariate)))) & is.na(match('HUGO_SYMBOL',toupper(names(covariate))))){
    stop("covariate lacks 'Hugo_Symbol' or 'Gene'")
  }else if(!is.na(match('GENE',toupper(names(covariate)))) & is.na(match('HUGO_SYMBOL',toupper(names(covariate))))){
    names(covariate)[match('GENE',toupper(names(covariate)))] <- 'Hugo_Symbol'
  }else if(!is.na(match('HUGO_SYMBOL',toupper(names(covariate)))) & is.na(match('GENE',toupper(names(covariate))))){
    names(covariate)[match('HUGO_SYMBOL',toupper(names(covariate)))] <- 'Hugo_Symbol'
  }else{
    names(covariate)[match('HUGO_SYMBOL',toupper(names(covariate)))] <- 'Hugo_Symbol'
    covariate <- subset(covariate, select = -match('GENE',toupper(names(covariate))))
  }

  if(!quiet){cat('success!\n')}

  # gene selection

  hgnc <- data.table::fread(system.file('extdata', 'hgnc.txt', package = 'DriverGenePathway'))
  covariate_right <- covariate[covariate$Hugo_Symbol %in% hgnc$symbol,]
  covariate_to_be_corrected <- covariate[!covariate$Hugo_Symbol %in% hgnc$symbol,]

  if(nrow(covariate_to_be_corrected)>0){

    if(!quiet){cat('Hugo symbol correction ... ')}

    covariate_to_be_corrected$old_symbol <- covariate_to_be_corrected$Hugo_Symbol
    covariate_to_be_corrected$Hugo_Symbol <- ''
    covariate_original_rows <- nrow(covariate)
    covariate_removed_rows <- sum(!covariate_to_be_corrected$old_symbol%in%hgnc$old_symbol)

    flag <- match(covariate_to_be_corrected$old_symbol, hgnc$old_symbol)
    covariate_to_be_corrected$Hugo_Symbol <- hgnc$symbol[flag]

    covariate_to_be_corrected <- covariate_to_be_corrected[!covariate_to_be_corrected$Hugo_Symbol=='',]
    covariate_corrected <- subset(covariate_to_be_corrected, select = -ncol(covariate_to_be_corrected))
    covariate <- rbind(covariate_right, covariate_corrected)

    if(!quiet){cat('success!\n')}

    if(covariate_removed_rows>0){
      if(covariate_removed_rows==1){
        if(!quiet){message(sprintf("NOTE: %d/%d Hugo_Symbol is not standard or alias symbol in the HGNC database and has been removed\n",
                        covariate_removed_rows, covariate_original_rows))}
      }else{
        if(!quiet){message(sprintf("NOTE: %d/%d Hugo_Symbols are not standard or alias symbol in the HGNC database and have been removed\n",
                        covariate_removed_rows, covariate_original_rows))}
      }
    }

    if(nrow(covariate) == 0){
      stop("There is no available data left!")
    }
  }else{
    if(!quiet){cat('All Hugo_Symbols look good.\n')}
  }

  # select hugo symbols in maf if maf exists
  if(!is.null(maf)){
    maf <- mafPreprocess(maf,output_file = FALSE, quiet = T)
    covariate <- covariate[covariate$Hugo_Symbol %in% maf$Hugo_Symbol,]
  }


  if(output_file){
    if(!dir.exists('DriverGenePathway_output')){
      dir.create('DriverGenePathway_output')
    }
    setwd('DriverGenePathway_output')
    utils::write.table(covariate,"preprocessed_covariate.txt", sep = "\t", quote = F, row.names = F)
    setwd("..")
  }

  if(!quiet){cat('Preprocess of covariate finishes.')}

  return(covariate)

}

