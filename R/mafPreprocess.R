#' Preprocess maf
#'
#' This function guarantees that the maf is available for DriverGenePathway. Two kinds of maf input
#' are feasible. It can be whether a txt file or an R data frame. The function returns the preprocessed
#' maf. An txt file of the preprocessed maf is optional to be exported.
#'
#' @param maf Mutation Annotation Format, can be whether an R data frame or the path of a txt file.
#' @param output_file Determine whether to export the preprocessed maf as a txt file.
#' @return The preprocessed maf.
#' @examples
#' laml_maf <- system.file("extdata", "tcga_laml.maf", package = "DriverGenePathway")
#' maf <- mafPreprocess(laml_maf)
#' @export

mafPreprocess <- function(maf, output_file = TRUE){

  cat('Loading maf ... ')

  if(is.character(maf)){
    maf <- data.table::fread(maf)
  }
  cat('success!\n')

  # ensure necessary elements in maf file
  # Hugo_Symbol
  cat('Check maf elements ... ')

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
  cat('success!\n')

  cat('Process variant classification ... ')

  dict_dir <- system.file("extdata", "mutation_type_dictionary_file.txt", package = "DriverGenePathway")
  dict <- data.table::fread(dict_dir)

  flag_num <- match(toupper(maf$Variant_Classification),
                    toupper(dict$Variant_Classification),nomatch =nrow(dict)+1)
  dict <- rbind(dict,dict[nrow(dict),])
  dict[nrow(dict),] <- "unknown"
  maf$effect <- dict$effect[flag_num]
  bad <- which(maf$effect == "unknown")
  if(length(bad)>0){
    if(length(bad)==1){
      cat(sprintf("WARNING: %d/%d Variant_Classification cannot be matched and has been removed\n",
                length(bad),length(maf$effect)))
    }else{
      cat(sprintf("WARNING: %d/%d Variant_Classifications cannot be matched and have been removed\n",
                  length(bad),length(maf$effect)))
    }
    maf <- maf[-bad]
  }
  if(nrow(maf) == 0){
    stop("There is no available data left!")
  }

  cat('success!\n')

  # output preprocessed maf
  if(output_file){
    if(!dir.exists('DriverGenePathway_output')){
      dir.create('DriverGenePathway_output')
      }
    setwd('DriverGenePathway_output')
    utils::write.table(maf,"preprocessed_maf.txt", sep = "\t", quote = F, row.names = F)
    setwd("..")
  }

  cat('Preprocess of maf finishes.')
  return(maf)
}

