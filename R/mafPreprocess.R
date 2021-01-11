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

  cat('Loading maf...')

  if(is.character(maf)){
    maf <- data.table::fread(maf)
  }
  cat('success!\n')

  # ensure necessary elements in maf file
  # Hugo_Symbol
  cat('Check maf elements...')

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

  cat('Process variant classification...')
  dict <- dplyr::tibble(Variant_Classification = c("Silent","Synonymous","Missense",
          "Missense_Mutation","Nonsense","Nonsense_Mutation","Nonstop_Mutation",
          "Read-through","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
          "Splice","Splice_Region","Splice_Site","Splice_Site_Del","Splice_Site_DNP",
          "Splice_Site_Ins","Splice_Site_ONP","Splice_Site_SNP","Start_Codon_Del",
          "Start_Codon_DNP","Start_Codon_Ins","Start_Codon_ONP","Stop_Codon_Del",
          "Stop_Codon_DNP","Stop_Codon_Ins","Translation_Start_Site","De_novo_Start",
          "De_novo_Start_InFrame","De_novo_Start_OutOfFrame","IGR","Intron","3'Flank",
          "3'Promoter","3'UTR","3'-UTR","5'Flank","5'-Flank","5'Promoter","5'UTR","5'-UTR",
          "downstream","miRNA","NCSD","Non-coding_Transcript","Promoter","RNA","upstream",
          "upstream;downstream"), effect = c("silent","silent","nonsilent","nonsilent",
          "null","null","null","null","null","null","null","null","null","null","null",
          "null","null","null","null","null","null","null","null","null","null","null",
          "null","null","null","null","null","noncoding","noncoding","noncoding","noncoding",
          "noncoding","noncoding","noncoding","noncoding","noncoding","noncoding","noncoding",
          "noncoding","noncoding","noncoding","noncoding","noncoding","noncoding","noncoding",
          "noncoding"))

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
    write.table(maf,"preprocessed_maf.txt", sep = "\t", quote = F, row.names = F)
    setwd("..")
  }

  cat('Preprocess of maf finishes.')
  return(maf)
}

