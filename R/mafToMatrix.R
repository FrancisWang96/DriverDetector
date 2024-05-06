#'  Transform maf to mutation matrix
#'
#' This function transforms maf files to mutation matrix in order to search driver pathways.
#'
#' @param maf Mutation Annotation Format, can be whether an R data frame or the path of a txt file.
#' @param threshold Minimum mutation count for genes to be analysed, default to be 5.
#' @param analyse_genes Genes to be analysed.
#' @param exclude_genes Genes excluded from analysis.
#' @param analyse_patients Patrients to be analysed.
#' @param exclude_patients Patrients excluded from analysis.
#' @param output_file Determine whether to export the preprocessed maf as a txt file.
#' @param quiet Whether to show notes during processing.
#' @examples
#' laml_maf <- system.file("extdata", "tcga_laml.maf", package = "DriverGenePathway")
#' mafToMatrix(maf = laml_maf)
#' @export


mafToMatrix <- function(maf = NULL, threshold = 5, analyse_genes = NULL, exclude_genes = NULL,
                        analyse_patients = NULL, exclude_patients = NULL, output_file = TRUE, quiet = FALSE){
  if(!quiet){cat('Loading maf ... ')}

  if(is.character(maf)){
    maf <- data.table::fread(maf)
  }
  if(!quiet){cat('success!\n')}

  temp <- dplyr::select(maf, Hugo_Symbol, Tumor_Sample_Barcode)
  # 删掉重复
  temp$merged <- stringr::str_c(temp$Hugo_Symbol, temp$Tumor_Sample_Barcode, sep = "_")
  temp1 <- unique(temp$merged)
  temp2 <- stringr::str_split_fixed(temp1, "_", n=2)
  temp2 <- dplyr::as_tibble(temp2)

  # 删掉threshold次以下
  temp2_gene_count <- dplyr::as_tibble(table(temp2$V1), .name_repair = unique)
  names(temp2_gene_count) <- c('c1','c2')
  flag <- which(temp2_gene_count$c2>=threshold)
  selected_genes <- temp2_gene_count$c1[flag]

  # 筛选基因
  if(!is.null(analyse_genes)){
    selected_genes <- selected_genes[selected_genes %in% analyse_genes]
  }
  # 删掉无效基因
  if(is.null(exclude_genes)){
    exclude_genes <- c('KRTAP5-5', 'HEATR7B2', 'EML2', 'CC2D1A', 'DUS3L',
                     'GABRA6', 'FLG', 'HYDIN', 'SUSD3', 'CROCCL1')
  }
  if(is.vector(exclude_genes)){
    selected_genes <- selected_genes[which(!selected_genes %in% exclude_genes)]
  }

  selected_patients = NULL

  # 筛选病人
  patients <- unique(temp2$V2)
  if(!is.null(analyse_patients)){
    selected_patients <- patients[patients %in% analyse_patients]
  }
  # 删掉无效病人
  if(!is.null(exclude_patients)){
    selected_patients <- patients[!patients %in% analyse_patients]

  }

  if(!is.null(selected_genes)){
    pre_mutation_matrix <- temp[temp$Hugo_Symbol %in% selected_genes, c(1,2)]
  }

  if(!is.null(selected_patients)){
    pre_mutation_matrix <- pre_mutation_matrix[Tumor_Sample_Barcode %in% selected_patients, ]
  }

  Hugo_Symbol <- unique(pre_mutation_matrix$Hugo_Symbol)
  Hugo_Symbol_table <- dplyr::tibble(Hugo_Symbol = Hugo_Symbol, flag = 1:length(Hugo_Symbol))

  Tumor_Sample_Barcode <- unique(pre_mutation_matrix$Tumor_Sample_Barcode)
  Tumor_Sample_Barcode_table <- dplyr::tibble(Tumor_Sample_Barcode = Tumor_Sample_Barcode, flag = 1:length(Tumor_Sample_Barcode))

  mutation_matrix <- matrix(0,length(Tumor_Sample_Barcode) , length(Hugo_Symbol))
  row.names(mutation_matrix) <- Tumor_Sample_Barcode
  colnames(mutation_matrix) <- Hugo_Symbol

  pre_mutation_matrix$r <- Tumor_Sample_Barcode_table$flag[match(pre_mutation_matrix$Tumor_Sample_Barcode,
                                                                 Tumor_Sample_Barcode_table$Tumor_Sample_Barcode)]
  pre_mutation_matrix$c <- Hugo_Symbol_table$flag[match(pre_mutation_matrix$Hugo_Symbol,
                                                        Hugo_Symbol_table$Hugo_Symbol)]
  if(!quiet){cat('Process: ')}
  for (i in 1:nrow(pre_mutation_matrix)) {
    if(i %% floor(nrow(pre_mutation_matrix)/10) == 0){
      if(!quiet){cat(stringr::str_c(i %/% floor(nrow(pre_mutation_matrix)/10),'0% '))}
    }
    mutation_matrix[pre_mutation_matrix$r[i],pre_mutation_matrix$c[i]] <- 1
  }
  if(!quiet){cat('\n')}
  utils::write.table(mutation_matrix, 'mutation_matrix.txt', quote = F)
  return(mutation_matrix)
}
