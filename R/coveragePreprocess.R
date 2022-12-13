#' Preprocess coverage file
#'
#' This function guarantees that the coverage file is available for DriverGenePathway. Two kinds of coverage input
#' are feasible. It can be whether a txt file or an R data frame. The function returns the preprocessed
#' coverage. An txt file of the preprocessed coverage is optional to be exported.
#'
#' @param coverage Coverage file, can be whether an R data frame or the path of a txt file.
#' @param maf Mutation Annotation Format, can be whether an R data frame or the path of a txt file.
#' @param output_file Determine whether to export the preprocessed maf as a txt file.
#' @param quiet Whether to show notes during processing.
#' @return The preprocessed coverage file.
#' @examples
#' coverage <- system.file("extdata", "coverage.rda", package = "DriverGenePathway")
#' load(coverage)
#' laml_maf <- system.file("extdata", "tcga_laml.maf", package = "DriverGenePathway")
#' coverage <- coveragePreprocess(coverage, laml_maf)
#' @export

coveragePreprocess <- function(coverage = NULL, maf = NULL, output_file = TRUE, quiet = FALSE){

  old_symbol = NULL
  num = NULL

  if(!quiet){cat('Loading coverage ... ')}

  if(is.character(coverage)){
    coverage <- data.table::fread(coverage)
  }
  if(!quiet){cat('success!\n')}

  coverage <- as.data.frame(coverage)

  if(!quiet){cat('Check coverage elements ... ')}

  if(is.na(match('GENE',toupper(names(coverage)))) & is.na(match('HUGO_SYMBOL',toupper(names(coverage))))){
    stop("coverage lacks 'Hugo_Symbol' or 'Gene'")
  }else if(!is.na(match('GENE',toupper(names(coverage)))) & is.na(match('HUGO_SYMBOL',toupper(names(coverage))))){
    names(coverage)[match('GENE',toupper(names(coverage)))] <- 'Hugo_Symbol'
  }else if(!is.na(match('HUGO_SYMBOL',toupper(names(coverage)))) & is.na(match('GENE',toupper(names(coverage))))){
    names(coverage)[match('HUGO_SYMBOL',toupper(names(coverage)))] <- 'Hugo_Symbol'
  }else{
    names(coverage)[match('HUGO_SYMBOL',toupper(names(coverage)))] <- 'Hugo_Symbol'
    coverage <- subset(coverage, select = -match('GENE',toupper(names(coverage))))
  }

  if(is.na(match('EFFECT',toupper(names(coverage))))){
    stop("coverage lacks 'effect', which indicates the area where mutation happens")
  }

  if(is.na(match('CATEG',toupper(names(coverage)))) & is.na(match('CATEGORY',toupper(names(coverage))))){
    stop("coverage lacks 'categ' or 'category'")
  }else if(!is.na(match('CATEG',toupper(names(coverage)))) & is.na(match('CATEGORY',toupper(names(coverage))))){
    names(coverage)[match('CATEG',toupper(names(coverage)))] <- 'categ'
  }else if(!is.na(match('CATEGORY',toupper(names(coverage)))) & is.na(match('CATEG',toupper(names(coverage))))){
    names(coverage)[match('CATEGORY',toupper(names(coverage)))] <- 'categ'
  }else{
    names(coverage)[match('CATEG',toupper(names(coverage)))] <- 'categ'
    coverage <- subset(coverage, select = -match('CATEGORY',toupper(names(coverage))))
  }

  if(is.na(match('COVERAGE',toupper(names(coverage))))){
    stop("coverage lacks 'coverage', which indicates the area where mutation happens")
  }

  if(!quiet){cat('success!\n')}

  # gene selection

  hgnc <- data.table::fread(system.file('extdata', 'hgnc.txt', package = 'DriverGenePathway'))
  hgnc_symbol <- hgnc$symbol
  coverage_symbol <- coverage$Hugo_Symbol
  coverage_right <- coverage[coverage_symbol %in% hgnc_symbol,]
  coverage_to_be_corrected <- coverage[!coverage_symbol %in% hgnc_symbol,]

  if(nrow(coverage_to_be_corrected)>0){

    if(!quiet){cat('Hugo symbol correction ... ')}

    coverage_to_be_corrected$old_symbol <- coverage_to_be_corrected$Hugo_Symbol
    coverage_to_be_corrected$Hugo_Symbol <- ''

    coverage_original_rows <- nrow(coverage)
    coverage_removed_rows <- sum(!coverage_to_be_corrected$old_symbol%in%hgnc$old_symbol)
    coverage_can_be_corrected <- coverage_to_be_corrected[which(coverage_to_be_corrected$old_symbol%in%hgnc$old_symbol),]

    flag <- match(coverage_can_be_corrected$old_symbol, hgnc$old_symbol)
    coverage_can_be_corrected$Hugo_Symbol <- hgnc$symbol[flag]

    coverage_corrected <- subset(coverage_can_be_corrected, select = -old_symbol)
    rm(coverage_can_be_corrected,coverage_to_be_corrected)
    coverage <- rbind(coverage_right,coverage_corrected)

    temp <- as.data.frame(table(coverage$Hugo_Symbol))
    gene_more <- as.vector(temp$Var1[which(temp$Freq>576)])
    if(length(gene_more)>0){
      coverage_corrected_right_gene <- coverage[which(!coverage$Hugo_Symbol%in%gene_more),]
      coverage_corrected_not_right_gene <- coverage[which(coverage$Hugo_Symbol%in%gene_more),]
      coverage_corrected_not_right_gene$num <- 0
      for (i in gene_more) {
        coverage_corrected_not_right_gene$num[which(coverage_corrected_not_right_gene$Hugo_Symbol==i)] <- 1:sum(coverage_corrected_not_right_gene$Hugo_Symbol==i)
      }
      coverage_corrected_not_right_gene <- coverage_corrected_not_right_gene[coverage_corrected_not_right_gene$num<=576,]
      coverage_corrected_not_right_gene <- subset(coverage_corrected_not_right_gene,select = -num)
      coverage <- rbind(coverage_corrected_right_gene, coverage_corrected_not_right_gene)
    }

    if(!quiet){cat('success!\n')}

    if(coverage_removed_rows>0){
      if(coverage_removed_rows/192==1){
        if(!quiet){message(sprintf("NOTE: %d/%d Hugo_Symbol is not standard or alias symbol in the HGNC database and has been removed\n",
                                   coverage_removed_rows/192, coverage_original_rows/192))}
      }else{
        if(!quiet){message(sprintf("NOTE: %d/%d Hugo_Symbols are not standard or alias symbol in the HGNC database and have been removed\n",
                                   coverage_removed_rows/192, coverage_original_rows/192))}
      }
    }

    if(nrow(coverage) == 0){
      stop("There is no available data left!")
    }
  }else{
    if(!quiet){cat('All Hugo_Symbols look good.\n')}
  }

  # select hugo symbols in maf if maf exists
  if(!is.null(maf)){
    maf <- mafPreprocess(maf,output_file = FALSE, quiet = T)
    coverage_symbol <- coverage$Hugo_Symbol
    maf_symbol <- maf$Hugo_Symbol
    coverage <- coverage[coverage_symbol %in% maf_symbol,]
  }




  if(output_file){
    if(!dir.exists('DriverGenePathway_output')){
      dir.create('DriverGenePathway_output')
    }
    setwd('DriverGenePathway_output')
    utils::write.table(coverage,"preprocessed_coverage.txt", sep = "\t", quote = F, row.names = F)
    setwd("..")
  }

  if(!quiet){cat('Preprocess of coverage finishes.')}

  return(coverage)

}

