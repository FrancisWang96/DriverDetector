#' Calculate background mutation rates
#'
#' This function uses the data output by the preprocessing function in order to get
#' background mutation rate. It uses the method of bagel. The result will be then
#' used in varified specific gene detecting methods.
#'
#' @param maf Mutation Annotation Format, can be whether an R data frame or the path of a txt file.
#' @param coverage Coverage file, can be whether an R data frame or the path of a txt file.
#' @param covariate Covariate file, can be whether an R data frame or the path of a txt file.
#' @param original_mutation_rate The default background mutation rate is 1.2e-6, which is used to
#' screen a subset of genes through binomial test.
#' @param max_neighbors Max number of gene neighbors used to calculate the background mutation rate
#' of a specified gene.
#' @param ref_genome Reference chromosome, either a folder name or a installed BSgenome package.
#' Default NULL, tries to auto-detect from installed genomes.
#' @param category_num Number of mutation categories, default 4. If the number is 0, categories should exist.
#' @param output_file Determine whether to export the categories, categorised maf and coverage as txt files.
#' @param quiet Whether to show notes during processing.
#' @details This function first selects a subset of genes whose mutation rates are significantly
#' greater than original mutation rate according to binomial test, then calculates background mutation
#' rates for these genes.
#' @return The output BMR_out is an input of the DriverGenes function.
#' @examples
#' laml_maf <- system.file("extdata", "tcga_laml.maf", package = "DriverGenePathway")
#' coverage <- system.file("extdata", "coverage.rda", package = "DriverGenePathway")
#' covariate <- system.file("extdata", "gene.covariates.txt", package = "DriverGenePathway")
#' load(coverage)
#' bmr_result <- backgroundMutationRate(laml_maf, coverage, covariate, original_mutation_rate = 1.2e-6,
#' max_neighbors = 50, ref_genome = NULL, category_num = 1, output_file = TRUE, quiet = FALSE)
#' @export

backgroundMutationRate <- function(maf = NULL, coverage = NULL, covariate = NULL,
                                   original_mutation_rate = 1.2e-6, max_neighbors = 50, ref_genome = NULL,
                                   category_num = NULL, output_file = TRUE, quiet = FALSE){
  options(digits=15)

  gene_idx = NULL
  background_mutaion_rate = NULL

  # 判断maf是否经过了categoryAssgin
  if(!"categ" %in% colnames(maf)){
    if(!quiet){cat('Maf lacks categories, now assigning...')}
    Categ_result <- mafCategoryAssign(maf, coverage, ref_genome = ref_genome,
                      category_num = category_num, output_file = output_file, quiet = quiet)
    maf <- Categ_result$maf
    coverage <- Categ_result$coverage
  }

  # 判断coverage是否经过了categoryAssgin
  if(!"categ" %in% colnames(coverage)||
     length(intersect(coverage$categ, maf$categ))!=length(unique(coverage$categ))){
    if(!quiet){cat('coverage lacks categories or categories mismatch with maf, now reassigning...')}
    Categ_result <- mafCategoryAssign(maf, coverage, ref_genome = ref_genome,
                                      category_num = category_num, output_file = output_file,
                                      quiet = quiet)
    maf <- Categ_result$maf
    coverage <- Categ_result$coverage
  }

  # 预处理coverage
  covariate <- covariatePreprocess(covariate, maf = maf, output_file = output_file, quiet = quiet)

  G <- data.frame(Hugo_Symbol=as.character(unique(coverage$Hugo_Symbol)))
  f <- colnames(covariate)
  cvnames <- f[2:ncol(covariate)]
  nv <- length(cvnames)
  gidx <- match(G$Hugo_Symbol,covariate$Hugo_Symbol)
  # 把covariate赋给G
  for (i in 1:nv){
    V <- covariate[[cvnames[i]]]
    G[[cvnames[i]]] <- V[gidx]
  }
  f <- colnames(coverage)
  coverage_patient_names <- f[4:ncol(coverage)]

  #remove any genes that we don't have coverage for
  bad_gene  <- setdiff(maf$Hugo_Symbol,coverage$Hugo_Symbol)
  if (length(bad_gene)!=0){
    cat(sprintf("NOTE: %d/%d Hugo_Symbol could not be mapped to coverage. Excluding them.\n",
                length(bad_gene), length(unique(maf$Hugo_Symbol))))
    flag_remove <- which(!(maf$Hugo_Symbol %in% bad_gene))
    maf <- maf[flag_remove,]
  }

  #category编号
  K <- data.frame(name = stringr::str_sort(unique(coverage$categ),locale = "C"))
  coverage$categ_idx <- match(coverage$categ,K$name)
  ncat <- length(K$name)
  maf$categ_idx <- match(maf$categ,K$name)

  #make sure there is a null+indel category
  if (is.numeric(unique(coverage$categ))){
    null_categ <- length(K$name)
  }else{
    null_categ <- grep('null|indel',K$name)
    if (length(null_categ)==0){
      stop("ERROR: no null/indel category. \n")
    }
    if (length(null_categ)>1){
      stop("ERROR: multiple null/indel categories. \n")
    }
  }

  # make sure C is sorted by the same Hugo_Symbol order as in G
  coverage$gene_idx <- match(coverage$Hugo_Symbol,G$Hugo_Symbol)
  coverage <- dplyr::arrange(coverage,gene_idx)

  #map genes
  maf$gene_idx <- match(maf$Hugo_Symbol,G$Hugo_Symbol)

  #regularize the sample name in the mutation table
  name_before <- maf$Tumor_Sample_Barcode
  maf$Tumor_Sample_Barcode <- gsub("-Tumor$","",maf$Tumor_Sample_Barcode)
  if (any(name_before != maf$Tumor_Sample_Barcode)){
    cat("NOTE: Trimming '-Tumor' from Tumor_Sample_Barcode. \n")
  }
  name_before <- maf$Tumor_Sample_Barcode
  maf$Tumor_Sample_Barcode <- gsub("-","_",maf$Tumor_Sample_Barcode)
  if (any(name_before != maf$Tumor_Sample_Barcode)){
    cat("NOTE: Converting '-' to '_' in Tumor_Sample_Barcode. \n")
  }

  Tumor_Sample_Barcode <- stringr::str_sort(unique(maf$Tumor_Sample_Barcode),locale = "C")
  pat <- data.frame(name = Tumor_Sample_Barcode)
  pat$cov_idx <- match(pat$name,coverage_patient_names)
  maf$patient_idx <- match(maf$Tumor_Sample_Barcode,Tumor_Sample_Barcode)
  np <- nrow(pat)
  if (np < 2){
    stop("This function is not applicable to single patients. \n")
  }

  #is generic coverage data given?
  generic_column_name <- "coverage"
  if ((length(coverage_patient_names)>1) || (length(coverage_patient_names)==1 &
                                             (coverage_patient_names[1]!=generic_column_name))){
    if (any(coverage_patient_names %in% generic_column_name)){
      stop(gettextf("reserved name '%s' cannot appear in list of patient names",
                    generic_column_name))
    }
    if (length(coverage_patient_names)!=length(unique(coverage_patient_names))){
      stop("patient names in coverage_file must be unique")
    }
    #make sure all patients are accounted for in coverage file
    if (any(is.na(pat$cov_idx))){
      stop("some patients in mutation_file are not accounted for in
                 coverage_file")
    }
    generic_coverage_flag <- FALSE
  }else{
    pat$cov_idx <- 1
    generic_coverage_flag <- TRUE
  }

  #BUILD n and N tables
  # cat("Building n and N tables...\n")

  Hugo_Symbol <- sort(as.character(unique(coverage$Hugo_Symbol)))
  categ_name <- as.character(K$name)
  categ_name[length(categ_name)+1] <- "total"
  ng <- length(unique(coverage$Hugo_Symbol))
  n_silent <- array(0,c(ng,ncat+1,np),dimnames=list(Hugo_Symbol,categ_name,Tumor_Sample_Barcode))
  n_nonsilent <- array(0,c(ng,ncat+1,np),dimnames=list(Hugo_Symbol,categ_name,Tumor_Sample_Barcode))
  n_noncoding <- array(0,c(ng,ncat+1,np),dimnames=list(Hugo_Symbol,categ_name,Tumor_Sample_Barcode))

  for (i in 1:ncat){
    for (j in 1:np){
      silent_hist <- graphics::hist(maf$gene_idx[(maf$effect %in% "silent") &
                                                 (maf$categ_idx %in% i) &
                                                 (maf$Tumor_Sample_Barcode %in% Tumor_Sample_Barcode[j])],
                                    c(0:length(Hugo_Symbol)), plot =FALSE)
      n_silent[,i,j] <-  silent_hist$counts
    }
  }

  for (i in 1:ncat){
    for (j in 1:np){
      nonsilent_hist <- graphics::hist(maf$gene_idx[((maf$effect %in% "nonsilent")|
                                                     (maf$effect %in% "null")) & (maf$categ_idx %in% i) &
                                                    (maf$Tumor_Sample_Barcode %in%
                                                       Tumor_Sample_Barcode[j])],c(0:length(Hugo_Symbol)),
                                       plot =FALSE)
      n_nonsilent[,i,j] <-  nonsilent_hist$counts
    }
  }


  for (i in 1:ncat){
    for (j in 1:np){
      noncoding_hist <- graphics::hist(maf$gene_idx[(maf$effect %in% "noncoding") &
                                                    (maf$categ_idx %in% i) &
                                                    (maf$Tumor_Sample_Barcode %in% Tumor_Sample_Barcode[j])],
                                       c(0:length(Hugo_Symbol)), plot =FALSE)
      n_noncoding[,i,j] <-  noncoding_hist$counts
    }
  }


  N_silent <- array(0,c(ng,ncat+1,np),dimnames=list(Hugo_Symbol,categ_name,Tumor_Sample_Barcode))
  N_nonsilent <- array(0,c(ng,ncat+1,np),dimnames=list(Hugo_Symbol,categ_name,Tumor_Sample_Barcode))
  N_noncoding <- array(0,c(ng,ncat+1,np),dimnames=list(Hugo_Symbol,categ_name,Tumor_Sample_Barcode))

  for (i in 1:ncat){
    for (j in 1:np)
    {
      N_silent[,i,j] <-  coverage$coverage[(coverage$categ %in% categ_name[i]) & (coverage$effect %in% "silent")]
      N_nonsilent[,i,j] <-  coverage$coverage[(coverage$categ %in% categ_name[i]) & (coverage$effect %in% "nonsilent")]
      N_noncoding[,i,j] <-  coverage$coverage[(coverage$categ %in% categ_name[i]) & (coverage$effect %in% "noncoding")]
    }
  }

  #MAKE SURE ALL NUMBERS ARE INTEGERS
  n_silent <- round(n_silent)
  n_nonsilent  <- round(n_nonsilent)
  n_noncoding  <- round(n_noncoding)
  N_silent  <- round(N_silent)
  N_nonsilent  <- round(N_nonsilent)
  N_noncoding  <- round(N_noncoding)

  #REMOVE MUTATIONS IN BINS WITH EXTREMELY LOW COVERAGE
  n_silent[n_silent>N_silent] <- 0
  n_nonsilent[n_nonsilent>N_nonsilent] <- 0
  n_noncoding[n_noncoding>N_noncoding] <- 0

  #SANITY CHECKS ON TOTALS
  tot_n_nonsilent <- sum(n_nonsilent)
  tot_N_nonsilent <- sum(N_nonsilent)
  tot_n_silent <- sum(n_silent)
  tot_N_silent <- sum(N_silent)
  tot_n_noncoding <- sum(n_noncoding)
  tot_N_noncoding <- sum(N_noncoding)
  tot_rate_nonsilent <- tot_n_nonsilent/tot_N_nonsilent
  tot_rate_silent <- tot_n_silent/tot_N_silent
  tot_rate_noncoding <- tot_n_noncoding/tot_N_noncoding
  tot_rate_coding <- (tot_n_nonsilent+tot_n_silent)/
    (tot_N_nonsilent+tot_N_silent)

  min_tot_n_nonsilent <- 50
  min_tot_n_silent <- 50
  min_tot_n_noncoding <- 50
  min_rate_nonsilent <- 1e-9
  max_rate_nonsilent <- 1e-3
  min_rate_silent <- 1e-9
  max_rate_silent <- 1e-3
  min_rate_noncoding <- 1e-9
  max_rate_noncoding <- 1e-3
  max_abs_log2_difference_nonsilent_silent <- 1.0
  max_abs_log2_difference_noncoding_coding <- 1.0

  #see if silent and nonsilent are OK: if not, give warning
  # if (tot_n_nonsilent<min_tot_n_nonsilent || tot_n_silent<min_tot_n_silent)
  #   stop("not enough mutations to analyze")
  if (tot_rate_nonsilent<min_rate_nonsilent || tot_rate_nonsilent>max_rate_nonsilent){
    stop("nonsilent mutation rate out of range")
  }
  if (tot_rate_silent<min_rate_silent || tot_rate_silent>max_rate_silent){
    stop("silent mutation rate out of range")
  }
  abs_log2_difference_nonsilent_silent = abs(log2(tot_rate_nonsilent/tot_rate_silent))
  if (abs_log2_difference_nonsilent_silent>max_abs_log2_difference_nonsilent_silent){
    cat('Warning: silent and nonsilent rates are too different.')
  }

  ## see if noncoding is OK: if not, give warning and zero it all out
  ok = FALSE
  if (tot_n_noncoding==0){
    cat('NOTE:  no noncoding mutations.')
  }else{
    if (tot_n_noncoding<min_tot_n_noncoding){
      cat('WARNING:  not enough noncoding mutations to analyze')
    }else{
        if (tot_rate_noncoding<min_rate_noncoding ||tot_rate_noncoding>max_rate_noncoding){
          cat('WARNING:  noncoding mutation rate out of range')
        }
        else{
          abs_log2_difference_noncoding_coding <- abs(log2(tot_rate_noncoding/tot_rate_coding))
        if (abs_log2_difference_noncoding_coding > max_abs_log2_difference_noncoding_coding){
          cat("WARNING:  coding and noncoding rates are too different")
        }else{
          ok = TRUE
        }
      }
    }
  }

  if (!ok){
    cat('Zeroing out all noncoding mutations and coverage for the rest
              of the calculation.\n')
    n_noncoding[,,]= 0
    N_noncoding[,,] = 0
  }

  #add total columns
  N_silent[,ncat+1,] <- N_silent[,null_categ,]
  N_nonsilent[,ncat+1,] <- N_nonsilent[,null_categ,]
  N_noncoding[,ncat+1,] <- N_noncoding[,null_categ,]
  n_silent[,ncat+1,] <- apply(n_silent,c(1,3),sum)
  n_nonsilent[,ncat+1,] <- apply(n_nonsilent,c(1,3),sum)
  n_noncoding[,ncat+1,] <- apply(n_noncoding,c(1,3),sum)

  #total across patients, save in G
  G$N_nonsilent <- apply(N_nonsilent[,ncat+1,],1,sum)
  G$N_silent <- apply(N_silent[,ncat+1,],1,sum)
  G$N_noncoding <- apply(N_noncoding[,ncat+1,],1,sum)
  G$n_nonsilent <- apply(n_nonsilent[,ncat+1,],1,sum)
  G$n_silent <- apply(n_silent[,ncat+1,],1,sum)
  G$n_noncoding <- apply(n_noncoding[,ncat+1,],1,sum)

  # PROCESS COVARIATES
  cat("Processing covariates... \n")

  V <- matrix(data = NaN,nrow = ng,ncol = nv)
  for (i in 1:nv) {
    V[,i] <- G[[cvnames[i]]]
  }
  colnames(V) <- cvnames

  # Find Bagels
  cat(sprintf("Finding bagels...  \n"))
  qual_min <- 0.05

  Z <- matrix(data = NA,nrow = nrow(V),ncol = ncol(V))

  for (vi in 1:nv) {
    missingv <- (is.na(V[,vi]) | is.infinite(V[,vi]))
    mn <- mean(V[!missingv,vi])
    sd <- stats::sd(V[!missingv,vi])
    Z[!missingv,vi] <- (V[!missingv,vi]-mn)/sd
  }
  #Z <- scale(V)

  G$n_neighbor <- NaN
  G$x <- NaN
  G$X <- NaN

  G1 <- G
  G1$pb <- NaN
  # 需要突变率大于背景突变率的基因，所以被择是greater
  for (i in 1:nrow(G)) {
    G1$pb[i] <- stats::binom.test(G1$n_nonsilent[i],G1$N_nonsilent[i],original_mutation_rate,
                                  alternative = "greater")$p.value
  }


  ord_pb <- order(G1$pb)
  G1 <- G1[ord_pb,]
  analysis_gene <- G1$Hugo_Symbol[which(G1$pb <= 0.05)]

  for (g in 1:length(analysis_gene)){
    if (g %% 1000 == 0) print(g)
    gi <- which(G$Hugo_Symbol %in% analysis_gene[g])

    df2= (Z-matrix(rep(Z[gi,],ng),nrow=ng,ncol=nv,byrow = TRUE))^2
    dfz  = df2
    dfz[is.na(dfz)] <-0
    dist2 <- apply(dfz,1,sum)/apply(!is.na(df2),1,sum)
    dist2 <- round(dist2,15)
    ord <- order(dist2,decreasing = FALSE,na.last = TRUE)
    ord <- c(gi,ord[ord!=gi])
    nfit <-0
    Nfit <-0
    for (ni in 0:max_neighbors)
    {
      gidx = ord[ni+1]
      ngene = G$n_silent[gidx]+ G$n_noncoding[gidx]
      Ngene = G$N_silent[gidx] + G$N_noncoding[gidx]
      if (ni==0)
      {
        ngene0=ngene
        Ngene0=Ngene
      }

      hc <- fun_hc(ngene,Ngene,ngene0,Ngene0)
      qual_left = min(hc,1-hc)
      qual <- 2*qual_left

      nfit=nfit+ngene
      Nfit=Nfit+Ngene
      if (ni>0 && qual<qual_min)
        break
      G$n_neighbor[gi] <- ni
      G$x[gi] <- nfit
      G$X[gi] <- Nfit
    }

  } # of for (gi in 1:ng)

  cat("Expanding to (x,X).gcp...\n")
  n_gcp = n_nonsilent + n_silent + n_noncoding
  N_gcp = N_nonsilent + N_silent + N_noncoding

  n_cp = apply(n_gcp,c(2,3),sum)
  N_cp = apply(N_gcp,c(2,3),sum)

  n_c = apply(n_cp,1,sum)
  N_c = apply(N_cp,1,sum)
  mu_c = n_c/N_c

  n_tot = n_c[length(n_c)]
  N_tot = N_c[length(N_c)]
  mu_tot = n_tot/N_tot
  f_c = mu_c/mu_tot
  f_Nc = N_c/N_tot

  n_p = n_cp[nrow(n_cp),]
  N_p = N_cp[nrow(N_cp),]
  mu_p = n_p/N_p
  f_p = mu_p/mu_tot
  f_Np = N_p/mean(N_p)

  x_gcp = array(rep(G$x,(ncat+1)*np),c(nrow(G),ncat+1,np))
  X_gcp = array(rep(G$X,(ncat+1)*np),c(nrow(G),ncat+1,np))

  for (ci in 1:(ncat+1)) {
    x_gcp[,ci,] <- x_gcp[,ci,]*(f_c[ci]*f_Nc[ci])
    X_gcp[,ci,] <- X_gcp[,ci,]*f_Nc[ci]
  }
  for (pi in 1:np) {
    x_gcp[,,pi] <- x_gcp[,,pi]*(f_p[pi]*f_Np[pi])
    X_gcp[,,pi] <- X_gcp[,,pi]*f_Np[pi]
  }

  G1 <- G[ord_pb,]
  G1$background_mutaion_rate <- G1$x/G1$X
  if(length(analysis_gene)>20){
    G1 <- G1[1:20,]
  }else{
    G1 <- G1[1:length(analysis_gene),]
  }
  G1 <- plyr::arrange(G1,plyr::desc(background_mutaion_rate))
  if(!quiet){cat('Background Mutation Rate result:\n')}
  if(!quiet){print(G1)}
  if(!quiet){cat('\n')}

  if(output_file){
    if(!dir.exists('DriverGenePathway_output')){
      dir.create('DriverGenePathway_output')
    }
    setwd('DriverGenePathway_output')
    utils::write.table(G1,"background_mutation_rate.txt", sep = "\t", quote = F, row.names = F)
    utils::write.table(analysis_gene, 'backgroundGenes.txt', row.names = F, col.names = F, quote = F)

    bmr_pic <- ggplot2::ggplot(G1,ggplot2::aes(x=Hugo_Symbol,y=background_mutaion_rate))+
      ggplot2::geom_point(size = 3, ggplot2::aes(color = background_mutaion_rate))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::geom_hline(yintercept=original_mutation_rate, linetype="dashed", color = "red")+
      ggplot2::geom_text(x = 10, y = original_mutation_rate, vjust = -1.5, label = "Original mutation rate", color="red", fontface = "italic", size = 5)+
      ggplot2::theme_minimal()+
      ggplot2::theme(axis.text.x = element_text(angle=25))+
      ggplot2::labs(title="Plot of background mutation rates", x ="Genes")

    ggplot2::ggsave(bmr_pic, file="background_mutation_rate.pdf", width = 9,height = 4)

    setwd("..")
  }

  return(list(G=G,maf=maf,N_silent=N_silent,n_silent=n_silent,
              N_nonsilent=N_nonsilent,n_nonsilent=n_nonsilent,
              N_noncoding=N_noncoding,n_noncoding=n_noncoding,
              X_gcp=X_gcp,x_gcp=x_gcp,null_categ=null_categ,
              analysis_gene=analysis_gene))
}
