#' Preprocess maf
#'
#' This function guarantees that the maf is available for DriverGenePathway. Two kinds of maf input
#' are feasible. It can be whether a txt file or an R data frame. The function returns the preprocessed
#' maf. An txt file of the preprocessed maf is optional to be exported.
#'
#' @param maf Mutation Annotation Format, can be whether an R data frame or the path of a txt file.
#' @param output_file Determine whether to export the preprocessed maf as a txt file.
#' @param quiet Whether to show notes during processing.
#' @return The preprocessed maf.
#' @examples
#' laml_maf <- system.file("extdata", "tcga_laml.maf", package = "DriverGenePathway")
#' maf <- mafPreprocess(laml_maf)
#' @export

mafPreprocess <- function(maf = NULL, output_file = TRUE, quiet = FALSE){

  old_symbol = NULL

  if(!quiet){cat('Loading maf ... ')}

  if(is.character(maf)){
    maf <- data.table::fread(maf)
  }
  if(!quiet){cat('success!\n')}

  # ensure necessary elements in maf file
  # Hugo_Symbol
  if(!quiet){cat('Check maf elements ... ')}

  if (!"Hugo_Symbol" %in% names(maf)){
    if(!"Gene" %in% names(maf) & !"gene" %in% names(maf)){
        stop("maf lacks 'Hugo_Symbol' or 'Gene'")
    }else if("Gene" %in% names(maf)){
      maf$Hugo_Symbol <- maf$Gene
    }else{
      maf$Hugo_Symbol <- maf$gene
    }
  }

  # Tumor_Sample_Barcode
  if (!"Tumor_Sample_Barcode" %in% names(maf)){
    stop("maf lacks 'Tumor_Sample_Barcode'")
  }

  if (length(unique(maf$Tumor_Sample_Barcode)) < 2){
    stop("DriverGenes is not applicable to single patients.")
  }

  if (!"Variant_Classification" %in% names(maf)){
    stop("maf lacks 'Variant_Classification'")
  }
  if(!quiet){cat('success!\n')}

  if(!quiet){cat('Process variant classification ... ')}

  dict_dir <- system.file("extdata", "mutation_type_dictionary_file.txt", package = "DriverGenePathway")
  dict <- data.table::fread(dict_dir)

  flag_num <- match(toupper(maf$Variant_Classification),
                    toupper(dict$Variant_Classification),nomatch =nrow(dict)+1)
  dict <- rbind(dict,dict[nrow(dict),])
  dict[nrow(dict),] <- "unknown"
  maf$effect <- dict$effect[flag_num]

  if(!quiet){cat('success!\n')}

  bad <- which(maf$effect == "unknown")
  if(length(bad)>0){
    if(length(bad)==1){
      if(!quiet){message(sprintf("NOTE: %d/%d Variant_Classification cannot be matched and has been removed\n",
                length(bad),length(maf$effect)))}
    }else{
      if(!quiet){message(sprintf("NOTE: %d/%d Variant_Classifications cannot be matched and have been removed\n",
                  length(bad),length(maf$effect)))}
    }
    maf <- maf[-bad]
  }
  if(nrow(maf) == 0){
    stop("There is no available data left!")
  }

  # gene selection

  hgnc <- data.table::fread(system.file('extdata', 'hgnc.txt', package = 'DriverGenePathway'))
  maf_right <- maf[maf$Hugo_Symbol %in% hgnc$symbol,]
  maf_to_be_corrected <- maf[!maf$Hugo_Symbol %in% hgnc$symbol,]

  if(nrow(maf_to_be_corrected)>0){

    if(!quiet){cat('Hugo symbol correction ... ')}

    maf_to_be_corrected$old_symbol <- maf_to_be_corrected$Hugo_Symbol
    maf_to_be_corrected$Hugo_Symbol <- ''

    maf_original_rows <- nrow(maf)
    maf_removed_rows <- sum(!maf_to_be_corrected$old_symbol%in%hgnc$old_symbol)
    maf_can_be_corrected <- maf_to_be_corrected[which(maf_to_be_corrected$old_symbol%in%hgnc$old_symbol),]

    flag <- match(maf_can_be_corrected$old_symbol, hgnc$old_symbol)
    maf_can_be_corrected$Hugo_Symbol <- hgnc$symbol[flag]

    maf_corrected <- subset(maf_can_be_corrected, select = -old_symbol)
    rm(maf_can_be_corrected,maf_to_be_corrected)
    maf <- rbind(maf_right, maf_corrected)

    if(!quiet){cat('success!\n')}

    if(maf_removed_rows>0){
      if(maf_removed_rows==1){
        if(!quiet){message(sprintf("NOTE: %d/%d Hugo_Symbol is not standard or alias symbol in the HGNC database and has been removed\n",
                                   maf_removed_rows, maf_original_rows))}
      }else{
        if(!quiet){message(sprintf("NOTE: %d/%d Hugo_Symbols are not standard or alias symbol in the HGNC database and have been removed\n",
                                   maf_removed_rows, maf_original_rows))}
      }
    }

    if(nrow(maf) == 0){
      stop("There is no available data left!")
    }
  }else{
    if(!quiet){cat('All Hugo_Symbols look good.\n')}
  }


  # output preprocessed maf
  if(output_file){
    if(!dir.exists('DriverGenePathway_output')){
      dir.create('DriverGenePathway_output')
      }
    setwd('DriverGenePathway_output')
    utils::write.table(maf,"preprocessed_maf.txt", sep = "\t", quote = F, row.names = F)
    setwd("..")
  }

  if(!quiet){cat('Preprocess of maf finishes.')}
  return(maf)
}

