#' Identify significante genes with hypothetical test methods.
#'
#' This function Applies different hypothetical test methods to identify significant genes. It outputs
#' the genes and their q-values.
#'
#' @param bmr A list output by the backgroundMutationRate function.
#' @param p_class Hypothetical test methods. "betaBinomial" represents beta binomial test;
#' "fisherBinomial" represents Fisher combined P-value test; "likelihoodRatio" represents likelihood
#' ratio test; "convolution" represents convolution test; "projection" represents 2D projection method;
#' "all" represents applying all above methods.
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
#' @param output_file Whether to export the results as csv files and pngs.
#' @param filter Whether to filter genes by q-values.
#' @param sigThreshold The threshhold of q-value to judge if the gene is significant
#' @param quiet Whether to show notes during processing.
#' @details This function searches the significant genes using different hypothetical test methods,
#' including beta binomial test, Fisher combined P-value test, likelihood ratio test, convolution test and
#' 2D projection method. There are two ways to run driverGenes. Since driverGenes can be treated as
#' the next step of the backgroundMutationRate function, the output of backgroundMutationRate can
#' one form of input. The other way is by inputting the maf, coverage and covariate files, and driverGenes
#' automaticaly runs the backgroundMutationRate procedure.
#' The output of backgroundMutationRate corresponds to the argument "bmr".
#' @return Signifiant genes
#' @examples
#' laml_maf <- system.file("extdata", "tcga_laml.maf", package = "DriverGenePathway")
#' coverage <- system.file("extdata", "coverage.rda", package = "DriverGenePathway")
#' covariate <- system.file("extdata", "gene.covariates.txt", package = "DriverGenePathway")
#' load(coverage)
#' bmr_result <- backgroundMutationRate(laml_maf, coverage, covariate, original_mutation_rate = 1.2e-6,
#'                                      max_neighbors = 50, ref_genome = NULL, category_num = 1,
#'                                      output_file = TRUE, quiet = FALSE)
#' driverGenes(bmr_result, p_class = "betaBinomial", output_file = TRUE, filter = TRUE,
#'             sigThreshold = 0.1, quiet = FALSE)
#'
#' driverGenes(p_class = "betaBinomial", maf = laml_maf, coverage = coverage, covariate = covariate,
#' original_mutation_rate = 1.2e-6, max_neighbors = 50, ref_genome = NULL, category_num = 1,
#' output_file = TRUE, filter = TRUE, sigThreshold = 0.1, quiet = FALSE)
#' @export

driverGenes <- function(bmr = NULL, p_class = "all", maf = NULL, coverage = NULL, covariate = NULL,
                        original_mutation_rate = 1.2e-6, max_neighbors = 50, ref_genome = NULL,
                        category_num = NULL, output_file = TRUE, filter = TRUE, sigThreshold = 0.1, quiet = FALSE)
{

  Hugo_Symbol = NULL
  lgq.CT = NULL
  lgq.LRT = NULL
  lgq.PJ = NULL
  lgq.btBinom = NULL
  lgq.fisher = NULL
  q.btBinom = NULL
  q.ct = NULL
  q.fisher = NULL
  q.lrt = NULL
  q.projection = NULL
  p.btBinom = NULL
  q.betaBinomial = NULL
  q.convolution = NULL
  q.fisherBinomial = NULL
  q.likelihoodRatio = NULL
  q.stouffer = NULL
  lgq.stouffer = NULL
  q.logit = NULL
  lgq.logit = NULL
  q.brown = NULL
  lgq.brown = NULL
  q.kost = NULL
  lgq.kost = NULL
  q.harmonic = NULL
  lgq.harmonic = NULL
  method = NULL
  repel = NULL

  if(is.null(bmr)){
    if(!quiet){cat('Calculating Background Mutation Rates...\n')}
    bmr <- backgroundMutationRate(maf, coverage, covariate, original_mutation_rate,
                                  max_neighbors, ref_genome, category_num, output_file, quiet)
  }

  if(output_file){
    if(!dir.exists('DriverGenePathway_output')){
      dir.create('DriverGenePathway_output')
    }
    setwd('DriverGenePathway_output')
  }
  options(digits = 15)
  G <- bmr$G
  N_silent_gc <- apply(bmr$N_silent,c(1,2),sum)
  N_nonsilent_gc <- apply(bmr$N_nonsilent,c(1,2),sum)
  N_noncoding_gc <- apply(bmr$N_noncoding,c(1,2),sum)
  n_silent_gc <- apply(bmr$n_silent,c(1,2),sum)
  n_nonsilent_gc <- apply(bmr$n_nonsilent,c(1,2),sum)
  n_noncoding_gc <- apply(bmr$n_noncoding,c(1,2),sum)

  alln_gc <- n_silent_gc + n_nonsilent_gc + n_noncoding_gc
  allN_gc <- N_silent_gc + N_nonsilent_gc + N_noncoding_gc

  X_gc <- apply(bmr$X_gcp,c(1,2),sum)
  x_gc <- apply(bmr$x_gcp,c(1,2),sum)

  # use different BMR for different genes
  cat_num <- ncol(x_gc)-1

  BMR_N <- X_gc
  BMR_n <- x_gc
  BMR <- x_gc/X_gc

  N_nonsilent_g <- G$N_nonsilent
  n_nonsilent_g <- G$n_nonsilent
  N_silent_g <- G$N_silent
  n_silent_g <- G$n_silent
  N_noncoding_g <- G$N_noncoding
  n_noncoding_g <- G$n_noncoding
  X_g <- G$X
  x_g <- G$x
  analysis_gene <- bmr$analysis_gene
  bmr_maf <- bmr$maf

  if(!p_class %in% c("betaBinomial","fisherBinomial","likelihoodRatio","convolution","projection","stouffer","logit","brown","kost","harmonic","denovo","all")){
    stop("p_class must be one of 'betaBinomial','fisherBinomial','likelihoodRatio','convolution','projection', 'stouffer','logit','brown','kost','harmonic','denovo', or 'all'")
  }

  if(p_class == "betaBinomial"){

    sig_betabinomial <- as.data.frame(matrix(numeric(0),nrow = length(x_g),ncol = 3))
    colnames(sig_betabinomial) <- c("Hugo_Symbol","p.btBinom","q.btBinom")
    sig_betabinomial$Hugo_Symbol <- G$Hugo_Symbol
    sig_betabinomial$p.btBinom <- NA
    for (g in analysis_gene) {
      i <- which(G$Hugo_Symbol %in% g)
      sig_betabinomial$p.btBinom[i] <- 1-fun_hc(n_nonsilent_g[i],
                                                N_nonsilent_g[i],
                                                x_g[i],X_g[i])

    }
    sig_betabinomial <- sig_betabinomial[which(!is.na(sig_betabinomial$p.btBinom)),]
    sig_betabinomial$q.btBinom <- stats::p.adjust(sig_betabinomial$p.btBinom,method="BH",nrow(G))
    if(filter == TRUE){
      sig_betabinomial <- sig_betabinomial[which(sig_betabinomial$q.btBinom <= sigThreshold),]

    }

    sig_betabinomial <- sig_betabinomial[order(sig_betabinomial$q.btBinom,decreasing = F),]
    if(output_file){
      utils::write.csv(sig_betabinomial,file = "beta_binomial_significant_genes.csv",row.names = F)
    }
    result_BB <- subset(sig_betabinomial,q.btBinom<0.1)
    if(nrow(result_BB)>30){
      result_BB1 <- result_BB[1:30,]
    }else{
      result_BB1 <- result_BB
    }
    result_BB1$lgq.btBinom <- -log10(result_BB1$q.btBinom+0.00000001)
    BB <- ggplot2::ggplot(result_BB1,ggplot2::aes(x=Hugo_Symbol,y=lgq.btBinom))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.btBinom))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by beta-binomial", x ="Genes", y = "-log(q.beta-binomial)")

    if(output_file){
      ggplot2::ggsave(BB, file="beta_binomial_significant_genes.pdf", width = 20)
    }
  }

  if(p_class == "fisherBinomial"){

    # Fisher combined binomial p-value
    sig_fisher <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = 3))
    colnames(sig_fisher) <- c("Hugo_Symbol","p.fisher","q.fisher")
    sig_fisher$Hugo_Symbol <- G$Hugo_Symbol
    sig_fisher$p.fisher <- NA
    FCPT_binom <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))

    for(g in analysis_gene){
      g <- which(G$Hugo_Symbol %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_fisher[g,]$p.fisher <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & (N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num]) & BMR[g,1:cat_num]>0)
        for(i in flag){
          FCPT_binom[g,i] <- stats::binom.test(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],alternative = "greater")$p.value
        }
        q = ( -2 ) * sum( log( FCPT_binom[g,flag] ))
        df = 2 * cat_num
        sig_fisher[g,]$p.fisher <- 1 - stats::pchisq( q, df )
      }
    }
    sig_fisher <- sig_fisher[which(!is.na(sig_fisher$p.fisher)),]
    stats::p.adjust( sig_fisher$p.fisher, method="BH" ,nrow(G)) -> sig_fisher$q.fisher

    if(filter == TRUE){
      sig_fisher <- sig_fisher[which(sig_fisher$q.fisher <= sigThreshold),]
    }

    sig_fisher <- sig_fisher[order(sig_fisher$q.fisher,decreasing = F),]
    if(output_file){
      utils::write.csv(sig_fisher,file = "fisher_combined_binomial_significant_genes.csv",row.names = F)
    }

    result_FCPT <- subset(sig_fisher,q.fisher<0.1)
    if(nrow(result_FCPT)>30){
      result_FCPT1 <- result_FCPT[1:30,]
    }else{
      result_FCPT1 <- result_FCPT
    }
    result_FCPT1$lgq.fisher <- -log10(result_FCPT1$q.fisher+0.00000001)
    FCPT <- ggplot2::ggplot(result_FCPT1,ggplot2::aes(x=Hugo_Symbol,y=lgq.fisher))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.fisher))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by fisher combined binomial", x ="Genes", y = "-log(q.fisher)")

    if(output_file){
      ggplot2::ggsave(FCPT, file="fisher_combined_binomial_significant_genes.pdf", width = 20)
    }

  } # end of if(p_class == "fisher")


  if(p_class == "stouffer"){

    # Fisher combined binomial p-value
    sig_stouffer <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = 3))
    colnames(sig_stouffer) <- c("Hugo_Symbol","p.stouffer","q.stouffer")
    sig_stouffer$Hugo_Symbol <- G$Hugo_Symbol
    sig_stouffer$p.stouffer <- NA
    stouffer_binom <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))

    for(g in analysis_gene){
      g <- which(G$Hugo_Symbol %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_stouffer[g,]$p.stouffer <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & (N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num]) & BMR[g,1:cat_num]>0)
        for(i in flag){
          stouffer_binom[g,i] <- stats::binom.test(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],alternative = "greater")$p.value
        }
        weight <- rep(1, length(stouffer_binom[g,flag]))
        z <- stats::qnorm(as.numeric(stouffer_binom[g,flag]))
        sig_stouffer[g,]$p.stouffer <- 1 - stats::pnorm(sum(weight * z) / sqrt(sum(weight^2)), lower.tail=FALSE)
      }
    }
    sig_stouffer <- sig_stouffer[which(!is.na(sig_stouffer$p.stouffer)),]
    stats::p.adjust( sig_stouffer$p.stouffer, method="BH" ,nrow(G)) -> sig_stouffer$q.stouffer

    if(filter == TRUE){
      sig_stouffer <- sig_stouffer[which(sig_stouffer$q.stouffer <= sigThreshold),]
    }

    sig_stouffer <- sig_stouffer[order(sig_stouffer$q.stouffer,decreasing = F),]
    if(output_file){
      utils::write.csv(sig_stouffer,file = "stouffer_significant_genes.csv",row.names = F)
    }

    result_stouffer <- subset(sig_stouffer,q.stouffer<0.1)
    if(nrow(result_stouffer)>30){
      result_stouffer1 <- result_stouffer[1:30,]
    }else{
      result_stouffer1 <- result_stouffer
    }
    result_stouffer1$lgq.stouffer <- -log10(result_stouffer1$q.stouffer+0.00000001)
    stouffer <- ggplot2::ggplot(result_stouffer1,ggplot2::aes(x=Hugo_Symbol,y=lgq.stouffer))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.stouffer))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by Stouffer's method", x ="Genes", y = "-log(q.stouffer)")
    if(output_file){
      ggplot2::ggsave(stouffer, file="stouffer_significant_genes.pdf", width = 20)
    }

  } # end of if(p_class == "stouffer")

  if(p_class == "logit"){

    # Fisher combined binomial p-value
    sig_logit <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = 3))
    colnames(sig_logit) <- c("Hugo_Symbol","p.logit","q.logit")
    sig_logit$Hugo_Symbol <- G$Hugo_Symbol
    sig_logit$p.logit <- NA
    logit_binom <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))

    for(g in analysis_gene){
      g <- which(G$Hugo_Symbol %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_logit[g,]$p.logit <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & (N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num]) & BMR[g,1:cat_num]>0)
        for(i in flag){
          logit_binom[g,i] <- stats::binom.test(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],alternative = "greater")$p.value
        }
        tt <- (- sum(log(logit_binom[g,flag] / (1 - logit_binom[g,flag])))) / sqrt(cat_num * pi^2 * (5 * cat_num + 2) / (3 * (5 * cat_num + 4)))
        sig_logit[g,]$p.logit <- stats::pt(tt,df=5*cat_num+4, lower.tail=FALSE)
      }
    }
    sig_logit <- sig_logit[which(!is.na(sig_logit$p.logit)),]
    stats::p.adjust( sig_logit$p.logit, method="BH" ,nrow(G)) -> sig_logit$q.logit

    if(filter == TRUE){
      sig_logit <- sig_logit[which(sig_logit$q.logit <= sigThreshold),]
    }

    sig_logit <- sig_logit[order(sig_logit$q.logit,decreasing = F),]
    if(output_file){
      utils::write.csv(sig_logit,file = "logit_significant_genes.csv",row.names = F)
    }

    result_logit <- subset(sig_logit,q.logit<0.1)
    if(nrow(result_logit)>30){
      result_logit1 <- result_logit[1:30,]
    }else{
      result_logit1 <- result_logit
    }
    result_logit1$lgq.logit <- -log10(result_logit1$q.logit+0.00000001)
    logit <- ggplot2::ggplot(result_logit1,ggplot2::aes(x=Hugo_Symbol,y=lgq.logit))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.logit))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by logit method", x ="Genes", y = "-log(q.logit)")
    if(output_file){
      ggplot2::ggsave(logit, file="logit_significant_genes.pdf", width = 20)
    }

  } # end of if(p_class == "fisher")

  if(p_class == "brown"){

    # Fisher combined binomial p-value
    sig_brown <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = 3))
    colnames(sig_brown) <- c("Hugo_Symbol","p.brown","q.brown")
    sig_brown$Hugo_Symbol <- G$Hugo_Symbol
    sig_brown$p.brown <- NA
    brown_binom <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))
    covariate <- system.file("extdata", "preprocessed_covariate.txt", package = "DriverGenePathway")
    covariate <- data.table::fread(covariate)

    for(gene in analysis_gene){
      if(!is.na(match(gene,covariate$Hugo_Symbol))){
        covariate_g <- covariate[match(gene,covariate$Hugo_Symbol),-1]
        g <- which(G$Hugo_Symbol %in% gene)
        if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
          sig_brown[g,]$p.brown <- 1
        }else{
          flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & (N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num]) & BMR[g,1:cat_num]>0)
          for(i in flag){
            brown_binom[g,i] <- stats::binom.test(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],alternative = "greater")$p.value
          }

          sig_brown[g,]$p.brown <- EmpiricalBrownsMethod::empiricalBrownsMethod(covariate_g, p_values = as.numeric(brown_binom[g,flag]), extra_info = FALSE)
        }
      }
    }
    sig_brown <- sig_brown[which(!is.na(sig_brown$p.brown)),]
    stats::p.adjust( sig_brown$p.brown, method="BH" ,nrow(G)) -> sig_brown$q.brown

    if(filter == TRUE){
      sig_brown <- sig_brown[which(sig_brown$q.brown <= sigThreshold),]
    }

    sig_brown <- sig_brown[order(sig_brown$q.brown,decreasing = F),]
    if(output_file){
      utils::write.csv(sig_brown,file = "brown_significant_genes.csv",row.names = F)
    }

    result_brown <- subset(sig_brown,q.brown<0.1)
    if(nrow(result_brown)>30){
      result_brown1 <- result_brown[1:30,]
    }else{
      result_brown1 <- result_brown
    }
    result_brown1$lgq.brown <- -log10(result_brown1$q.brown+0.00000001)
    brown <- ggplot2::ggplot(result_brown1,ggplot2::aes(x=Hugo_Symbol,y=lgq.brown))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.brown))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by brown method", x ="Genes", y = "-log(q.brown)")
    if(output_file){
      ggplot2::ggsave(brown, file="brown_significant_genes.pdf", width = 20)
    }

  } # end of if(p_class == "fisher")

  if(p_class == "kost"){

    # Fisher combined binomial p-value
    sig_kost <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = 3))
    colnames(sig_kost) <- c("Hugo_Symbol","p.kost","q.kost")
    sig_kost$Hugo_Symbol <- G$Hugo_Symbol
    sig_kost$p.kost <- NA
    kost_binom <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))
    covariate <- system.file("extdata", "preprocessed_covariate.txt", package = "DriverGenePathway")
    covariate <- data.table::fread(covariate)

    for(gene in analysis_gene){
      if(!is.na(match(gene,covariate$Hugo_Symbol))){
        covariate_g <- covariate[match(gene,covariate$Hugo_Symbol),-1]
        g <- which(G$Hugo_Symbol %in% gene)
        if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
          sig_kost[g,]$p.kost <- 1
        }else{
          flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & (N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num]) & BMR[g,1:cat_num]>0)
          for(i in flag){
            kost_binom[g,i] <- stats::binom.test(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],alternative = "greater")$p.value
          }
          covariate_g <- rbind(covariate_g,covariate_g)
          # for (i in 1:(flag-1)) {
          #   covariate_g1 <- rbind(covariate_g1,covariate_g)
          # }
          sig_kost[g,]$p.kost <- EmpiricalBrownsMethod::kostsMethod(as.matrix(covariate_g), p_values = as.numeric(kost_binom[g,flag]), extra_info = FALSE)
        }
      }
    }
    sig_kost <- sig_kost[which(!is.na(sig_kost$p.kost)),]
    stats::p.adjust( sig_kost$p.kost, method="BH" ,nrow(G)) -> sig_kost$q.kost

    if(filter == TRUE){
      sig_kost <- sig_kost[which(sig_kost$q.kost <= sigThreshold),]
    }

    sig_kost <- sig_kost[order(sig_kost$q.kost,decreasing = F),]
    if(output_file){
      utils::write.csv(sig_kost,file = "kost_significant_genes.csv",row.names = F)
    }

    result_kost <- subset(sig_kost,q.kost<0.1)
    if(nrow(result_kost)>30){
      result_kost1 <- result_kost[1:30,]
    }else{
      result_kost1 <- result_kost
    }
    result_kost1$lgq.kost <- -log10(result_kost1$q.kost+0.00000001)
    kost <- ggplot2::ggplot(result_kost1,ggplot2::aes(x=Hugo_Symbol,y=lgq.kost))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.kost))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by kost method", x ="Genes", y = "-log(q.kost)")
    if(output_file){
      ggplot2::ggsave(kost, file="kost_significant_genes.pdf", width = 20)
    }

  } # end of if(p_class == "fisher")

  if(p_class == "harmonic"){

    # Fisher combined binomial p-value
    sig_harmonic <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = 3))
    colnames(sig_harmonic) <- c("Hugo_Symbol","p.harmonic","q.harmonic")
    sig_harmonic$Hugo_Symbol <- G$Hugo_Symbol
    sig_harmonic$p.harmonic <- NA
    harmonic_binom <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))
    covariate <- system.file("extdata", "preprocessed_covariate.txt", package = "DriverGenePathway")
    covariate <- data.table::fread(covariate)

    for(gene in analysis_gene){
      if(!is.na(match(gene,covariate$Hugo_Symbol))){
        covariate_g <- covariate[match(gene,covariate$Hugo_Symbol),-1]
        g <- which(G$Hugo_Symbol %in% gene)
        if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
          sig_harmonic[g,]$p.harmonic <- 1
        }else{
          flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & (N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num]) & BMR[g,1:cat_num]>0)
          for(i in flag){
            harmonic_binom[g,i] <- stats::binom.test(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],alternative = "greater")$p.value
          }
          covariate_g <- rbind(covariate_g,covariate_g)
          # for (i in 1:(flag-1)) {
          #   covariate_g1 <- rbind(covariate_g1,covariate_g)
          # }
          sig_harmonic[g,]$p.harmonic <- harmonicmeanp::p.hmp(as.numeric(harmonic_binom[g,flag]), L=length(flag))

        }
      }
    }
    sig_harmonic <- sig_harmonic[which(!is.na(sig_harmonic$p.harmonic)),]
    stats::p.adjust( sig_harmonic$p.harmonic, method="BH" ,nrow(G)) -> sig_harmonic$q.harmonic

    if(filter == TRUE){
      sig_harmonic <- sig_harmonic[which(sig_harmonic$q.harmonic <= sigThreshold),]
    }

    sig_harmonic <- sig_harmonic[order(sig_harmonic$q.harmonic,decreasing = F),]
    if(output_file){
      utils::write.csv(sig_harmonic,file = "harmonic_significant_genes.csv",row.names = F)
    }

    result_harmonic <- subset(sig_harmonic,q.harmonic<0.1)
    if(nrow(result_harmonic)>30){
      result_harmonic1 <- result_harmonic[1:30,]
    }else{
      result_harmonic1 <- result_harmonic
    }
    result_harmonic1$lgq.harmonic <- -log10(result_harmonic1$q.harmonic+0.00000001)
    harmonic <- ggplot2::ggplot(result_harmonic1,ggplot2::aes(x=Hugo_Symbol,y=lgq.harmonic))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.harmonic))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by harmonic method", x ="Genes", y = "-log(q.harmonic)")
    if(output_file){
      ggplot2::ggsave(harmonic, file="harmonic_significant_genes.pdf", width = 20)
    }

  } # end of if(p_class == "fisher")

  if(p_class == "denovo"){

    driverPathway(maf = bmr_maf, threshold = 5, analyse_genes = analysis_gene,
                              driver_size=3,
                              pop_size=200,
                              iters=200,
                              permut_time=500,
                              output_file = output_file,
                              quiet = quiet)

  } # end of if(p_class == "fisher")

  if(p_class == "likelihoodRatio"){

    sig_lrt <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = 3))
    colnames(sig_lrt) <- c("Hugo_Symbol","p.lrt","q.lrt")
    sig_lrt$Hugo_Symbol <- G$Hugo_Symbol
    sig_lrt$p.lrt <- NA
    # Likelihood ratio test
    LRT_lh0 <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))
    LRT_lh1 <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))

    for(g in analysis_gene){
      g <- which(G$Hugo_Symbol %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_lrt[g,]$p.lrt <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num] & BMR[g,1:cat_num]>0)
        for(i in flag){
          LRT_lh0[g,i] <- stats::dbinom(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],log = T)
          BMR_raw_n <- n_nonsilent_gc[g,i] + n_silent_gc[g,i] + n_noncoding_gc[g,i]
          BMR_raw_N <- N_nonsilent_gc[g,i] + N_silent_gc[g,i] + N_noncoding_gc[g,i]
          LRT_lh1[g,i] <- stats::dbinom(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR_raw_n/BMR_raw_N,log = T)
        }
        q = 2 * ( sum( LRT_lh1[g,flag] ) - sum( LRT_lh0[g,flag]))
        df = sum( LRT_lh1[g,flag] != 0 );
        if( df > 0 ) p_lrt = 1 - stats::pchisq( q, df )
        if( df == 0 ) p_lrt = 1
        sig_lrt$p.lrt[g] <- p_lrt
      }
    }
    sig_lrt <- sig_lrt[which(!is.na(sig_lrt$p.lrt)),]
    stats::p.adjust(sig_lrt$p.lrt, method="BH" , nrow(G)) -> sig_lrt$q.lrt
    if(filter == TRUE){
      sig_lrt <- sig_lrt[which(sig_lrt$q.lrt <= sigThreshold),]
    }
    sig_lrt <- sig_lrt[order(sig_lrt$q.lrt,decreasing = F),]
    if(output_file){
      utils::write.csv(sig_lrt,file = "likelihood_ratio_significant_genes.csv",row.names = F)
    }

    result_LRT <- subset(sig_lrt,q.lrt<0.1)
    if(nrow(result_LRT)>30){
      result_LRT1 <- result_LRT[1:30,]
    }else{
      result_LRT1 <- result_LRT
    }
    result_LRT1$lgq.LRT <- -log10(result_LRT1$q.lrt+0.00000001)
    LRT <- ggplot2::ggplot(result_LRT1,ggplot2::aes(x=Hugo_Symbol,y=lgq.LRT))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.LRT))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by likelihood ratio", x ="Genes", y = "-log(q.lrt)")
    if(output_file){
      ggplot2::ggsave(LRT, file="likelihood_ratio_significant_genes.pdf", width = 20)
    }
  }

  if(p_class == "convolution"){

    # Convolution test
    sig_ct <- as.data.frame(matrix(data = NA,nrow = nrow(G),ncol = 3))
    colnames(sig_ct) <- c("Hugo_Symbol","p.ct","q.ct")
    sig_ct$Hugo_Symbol <- G$Hugo_Symbol
    sig_ct$p.ct <- NA

    LRT_lh0 <- as.data.frame(matrix(data = NA,nrow = nrow(G),ncol = cat_num))

    for(g in analysis_gene){
      g <- which(G$Hugo_Symbol %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_ct[g,]$p.ct <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num] & BMR[g,1:cat_num]>0)
        for(i in flag){
          LRT_lh0[g,i] <- stats::dbinom(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],log = T)
        }
      }
    }

    for(g in analysis_gene){
      g <- which(G$Hugo_Symbol %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_ct[g,]$p.ct <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num] & BMR[g,1:cat_num]>0)
        for(i in flag){
          ni <- N_nonsilent_gc[g,i]
          ei <- BMR[g,i]
          xmax = 100
          hmax = 25
          bin = 0.001
          gethist(xmax,ni,ei,ptype = "positive_log") -> bi
          binit(bi,hmax, bin ) -> bi
          if( i == flag[1] ) { hist0 = bi; }
          if( i > flag[1] & i < flag[length(flag)] ) { hist0 = convolute_b( hist0, bi ); binit( hist0, hmax, bin ) -> hist0 }
          if( i ==  flag[length(flag)]) { hist0 = convolute_b( hist0, bi ) }
        }
        # Convolution test
        bx = -sum( LRT_lh0[g,flag] )
        p.ct = sum( exp( -hist0[hist0>=bx] ))
        qc = sum( exp( -hist0 ))
        sig_ct$p.ct[g] <- p.ct
      }
    }
    sig_ct <- sig_ct[which(!is.na(sig_ct$p.ct)),]
    stats::p.adjust( sig_ct$p.ct, method="BH" ,nrow(G)) -> sig_ct$q.ct

    if(filter == TRUE){
      sig_ct <- sig_ct[which(sig_ct$q.ct <= sigThreshold),]
    }

    sig_ct <- sig_ct[order(sig_ct$q.ct,decreasing = F),]
    if(output_file){
      utils::write.csv(sig_lrt,file = "convolution_significant_genes.csv",row.names = F)
    }

    result_CT <- subset(sig_ct,q.ct<0.1)
    if(nrow(result_CT)>30){
      result_CT1 <- result_CT[1:30,]
    }else{
      result_CT1 <- result_CT
    }
    result_CT1$lgq.CT <- -log10(result_CT1$q.ct+0.00000001)
    CT <- ggplot2::ggplot(result_CT1,ggplot2::aes(x=Hugo_Symbol,y=lgq.CT))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.CT))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by convolution", x ="Genes", y = "-log(q.ct)")
    if(output_file){
      ggplot2::ggsave(CT, file="convolution_significant_genes.pdf", width = 20)
    }
  }

  if(p_class == "projection"){

    N_nonsilent=bmr$N_nonsilent
    n_nonsilent=bmr$n_nonsilent
    X_gcp=bmr$X_gcp
    x_gcp=bmr$x_gcp
    null_categ=bmr$null_categ

    ######### 2D projection to calculate p value and q ##########

    #PROJECTION
    if(!quiet){cat("Calculating p.value using 2D Projection method...\n")}
    null_score_boost = 3
    min_effect_size =1.25
    convolution_numbins =1000

    ncat <- dim(X_gcp)[2]-1
    np <- dim(X_gcp)[3]

    G$p <-NaN
    for (g in 1:length(analysis_gene)){
      if (g %% 1000 == 0) print(g)

      g <- which(G$Hugo_Symbol %in% analysis_gene[g])
      P0 <- matrix(0,nrow=np,ncol = ncat)
      P1 <- matrix(0,nrow=np,ncol = ncat)
      priority <- matrix(0,nrow=np,ncol = ncat)
      N <- t(N_nonsilent[g,1:ncat,])
      n <- t(n_nonsilent[g,1:ncat,])
      x <- t(x_gcp[g,1:ncat,])
      X <- t(X_gcp[g,1:ncat,])

      P0 <- fun_hp(0,N,x,X)
      P1 <- fun_hp(1,N,x,X)
      priority <- t(apply(P1,1,order,decreasing=T))

      shft <- priority - matrix(rep(1:ncat,np),nrow=np,ncol=ncat,byrow = TRUE)
      map <- matrix(1:(np*ncat),nrow = np,ncol = ncat,byrow = FALSE)
      newmap <- map +shft*np
      P0 <- matrix(P0[newmap[1:length(newmap)]],nrow = np)
      P1 <- matrix(P1[newmap[1:length(newmap)]],nrow = np)
      P2 <- 1-(P0+P1)
      P2[P2<0] <- 0

      Pdeg <- array(0,c(np,ncat+1,ncat+1))

      for (d1 in 0:ncat){
        for (d2 in 0:d1)
        {
          if (d1 == ncat)
            p <- matrix(1,nrow = np,ncol = 1)
          else
            p <- apply(as.matrix(P0[,(d1+1):ncol(P0)]),1,prod)
          if (d1>0)
          { if (d1==d2)
            p = p*P2[,d1]
          else
          {
            p = p * P1[,d1]
            if (d2<(d1-1))
              p = p * apply(as.matrix(P0[,(d2+1):(d1-1)]),1,prod)
            if (d2>0)
              p = p * (P1[,d2]+P2[,d2])
          }
          }
          Pdeg[,d1+1,d2+1] <- p
        }
      }


      Sdeg = array(0,c(np,ncat+1,ncat+1))
      for (d1 in 1:ncat)
      {
        for (d2 in 0:d1)
        {
          if (d1==d2){
            p = P2[,d1]
          }else
          {
            if (d2>0)
              p = P1[,d1]*P1[,d2]
            else
              p = P1[,d1]
          }

          Sdeg[,d1+1,d2+1] = -log10(p)
        }
      }

      priority2 <- cbind(matrix(0,nrow = np,ncol = 1),priority)

      #Sdeg[priority2==null_categ] = Sdeg[priority2==null_categ] + null_score_boost
      s1 <- Sdeg[,,1]
      Sdeg[,,1][priority2==null_categ] <- s1[priority2==null_categ]+null_score_boost
      degree <-  matrix(0,nrow=np,ncol=2)
      score_obs = 0

      for(p in  1:np)
      {
        i = 1
        for (d in ncat:1)
        {
          c = priority[p,d]
          if (i==1)
          {
            if (n[p,c]>=2)
            {
              degree[p,] = c(d,d)
              i = 3
            }
            else if (n[p,c]==1)
            {
              degree[p,i] = d
              i=i+1
            }
          }
          else if (i==2)
          {
            if (n[p,c]>=1)
            {
              degree[p,i] = d
              i=i+1
            }
            else
              break
          }
        }
        score_sample = Sdeg[p,degree[p,1]+1,degree[p,2]+1]
         if(!is.finite(score_sample)){
        #   cat(sprintf("g=%d",g))
        #   cat(sprintf("p=%d",p))
           flag_finite <- is.finite(Sdeg[,degree[p,1]+1,degree[p,2]+1])
           score_sample <- max(Sdeg[,degree[p,1]+1,degree[p,2]+1][flag_finite])
        #   cat(sprintf("maxscore=%f",score_sample))
         }
        score_obs = score_obs + score_sample
      }

      score_obs = score_obs / min_effect_size

      if (score_obs<=0)
      {
        G$p[g]=1
        next
      }
      numbins = convolution_numbins
      binsize = score_obs / numbins
      H = matrix(0,nrow=numbins,ncol=1)
      H[1] <- 1

      offset = pmin(array(numbins,c(np,ncat+1,ncat+1)), round(Sdeg/binsize))
      ncols = (ncat+1)*(ncat+2)/2
      newH = matrix(0,nrow=numbins,ncol=ncols)
      for (p in 1:np)
      {
        newH[,] <- 0
        col=1
        for (d1 in 0:ncat)
          for (d2 in 0:d1)
          {
            o = offset[p,d1+1,d2+1]
            if (o<length(H))
              newH[(o+1):nrow(newH),col] = Pdeg[p,d1+1,d2+1] * H[1:(length(H)-o)]
            col=col+1
          }
        H = apply(newH,1,sum,drop=FALSE)
      }

      G$p[g] = max(0,1-sum(H))

    }  # of for (g in 1:ng)
    nG <- nrow(G)
    G <- G[which(!is.na(G$p)),]
    G$q <- stats::p.adjust(G$p,method="BH",nG)
    G <- G[order(G$q),]
    sig_projection <- subset(G,select=c("Hugo_Symbol","p","q"))

    if(filter == TRUE){
      sig_projection <- sig_projection[which(sig_projection$q <= sigThreshold),]
    }

    colnames(sig_projection) <- c("Hugo_Symbol","p.projection","q.projection")
    sig_projection <- sig_projection[which(!is.na(sig_projection$p.projection)),]
    if(output_file){
      utils::write.csv(sig_projection, file = "2d_projection_significant_genes.csv",row.names = F)
    }

    result_PJ <- subset(sig_projection,q.projection<0.1)
    if(nrow(result_PJ)>30){
      result_PJ1 <- result_PJ[1:30,]
    }else{
      result_PJ1 <- result_PJ
    }
    gene_projection <- result_PJ[,1]
    result_PJ1$lgq.PJ <- -log10(result_PJ1$q.projection+0.00000001)
    PJ <- ggplot2::ggplot(result_PJ1,ggplot2::aes(x=Hugo_Symbol,y=lgq.PJ))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.PJ))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by 2D projection", x ="Genes", y = "-log(q.projection)")
    if(output_file){
      ggplot2::ggsave(PJ, file="2d_projection_significant_genes.pdf", width = 20)
    }
  }

  if(p_class == "all"){

    # betaBinomial
    sig_betabinomial <- as.data.frame(matrix(numeric(0),nrow = length(x_g),ncol = 3))
    colnames(sig_betabinomial) <- c("Hugo_Symbol","p.btBinom","q.btBinom")
    sig_betabinomial$Hugo_Symbol <- G$Hugo_Symbol
    sig_betabinomial$p.btBinom <- NA
    for (g in analysis_gene) {
      i <- which(G$Hugo_Symbol %in% g)
      sig_betabinomial$p.btBinom[i] <- 1-fun_hc(n_nonsilent_g[i],
                                                N_nonsilent_g[i],
                                                x_g[i],X_g[i])

    }
    sig_betabinomial <- sig_betabinomial[which(!is.na(sig_betabinomial$p.btBinom)),]
    sig_betabinomial$q.btBinom <- stats::p.adjust(sig_betabinomial$p.btBinom,method="BH",nrow(G))
    if(filter == TRUE){
      sig_betabinomial <- sig_betabinomial[which(sig_betabinomial$q.btBinom <= sigThreshold),]

    }

    sig_betabinomial <- sig_betabinomial[order(sig_betabinomial$q.btBinom,decreasing = F),]
    if(output_file){
      utils::write.csv(sig_betabinomial,file = "beta_binomial_significant_genes.csv",row.names = F)
    }
    result_BB <- subset(sig_betabinomial,q.btBinom<0.1)
    if(nrow(result_BB)>30){
      result_BB1 <- result_BB[1:30,]
    }else{
      result_BB1 <- result_BB
    }
    result_BB1$lgq.btBinom <- -log10(result_BB1$q.btBinom+0.00000001)
    BB <- ggplot2::ggplot(result_BB1,ggplot2::aes(x=Hugo_Symbol,y=lgq.btBinom))+
      ggplot2::geom_point(size = 3, ggplot2::aes(color = lgq.btBinom))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::theme(axis.text.x = element_text(angle=35))+
      ggplot2::labs(title="Plot of result genes tested by beta-binomial", x ="Genes", y = "-log(q.beta-binomial)")

    if(output_file){
      ggplot2::ggsave(BB, file="beta_binomial_significant_genes.pdf", width = 9, height = 4)
    }


    # Fisher combined binomial p-value
    sig_fisher <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = 3))
    colnames(sig_fisher) <- c("Hugo_Symbol","p.fisher","q.fisher")
    sig_fisher$Hugo_Symbol <- G$Hugo_Symbol
    sig_fisher$p.fisher <- NA
    FCPT_binom <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))

    for(g in analysis_gene){
      g <- which(G$Hugo_Symbol %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_fisher[g,]$p.fisher <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & (N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num]) & BMR[g,1:cat_num]>0)
        for(i in flag){
          FCPT_binom[g,i] <- stats::binom.test(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],alternative = "greater")$p.value
        }
        q = ( -2 ) * sum( log( FCPT_binom[g,flag] ))
        df = 2 * cat_num
        sig_fisher[g,]$p.fisher <- 1 - stats::pchisq( q, df )
      }
    }
    sig_fisher <- sig_fisher[which(!is.na(sig_fisher$p.fisher)),]
    stats::p.adjust( sig_fisher$p.fisher, method="BH" ,nrow(G)) -> sig_fisher$q.fisher

    if(filter == TRUE){
      sig_fisher <- sig_fisher[which(sig_fisher$q.fisher <= sigThreshold),]
    }

    sig_fisher <- sig_fisher[order(sig_fisher$q.fisher,decreasing = F),]
    if(output_file){
      utils::write.csv(sig_fisher,file = "fisher_combined_binomial_significant_genes.csv",row.names = F)
    }

    result_FCPT <- subset(sig_fisher,q.fisher<0.1)
    if(nrow(result_FCPT)>30){
      result_FCPT1 <- result_FCPT[1:30,]
    }else{
      result_FCPT1 <- result_FCPT
    }
    result_FCPT1$lgq.fisher <- -log10(result_FCPT1$q.fisher+0.00000001)
    FCPT <- ggplot2::ggplot(result_FCPT1,ggplot2::aes(x=Hugo_Symbol,y=lgq.fisher))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.fisher))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by fisher combined binomial", x ="Genes", y = "-log(q.fisher)")

    if(output_file){
      ggplot2::ggsave(FCPT, file="fisher_combined_binomial_significant_genes.pdf", width = 20)
    }

    # Stouffer's Z-score
    sig_stouffer <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = 3))
    colnames(sig_stouffer) <- c("Hugo_Symbol","p.stouffer","q.stouffer")
    sig_stouffer$Hugo_Symbol <- G$Hugo_Symbol
    sig_stouffer$p.stouffer <- NA
    stouffer_binom <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))

    for(g in analysis_gene){
      g <- which(G$Hugo_Symbol %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_stouffer[g,]$p.stouffer <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & (N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num]) & BMR[g,1:cat_num]>0)
        for(i in flag){
          stouffer_binom[g,i] <- stats::binom.test(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],alternative = "greater")$p.value
        }
        weight <- rep(1, length(stouffer_binom[g,flag]))
        z <- stats::qnorm(as.numeric(stouffer_binom[g,flag]))
        sig_stouffer[g,]$p.stouffer <- 1 - stats::pnorm(sum(weight * z) / sqrt(sum(weight^2)), lower.tail=FALSE)
      }
    }
    sig_stouffer <- sig_stouffer[which(!is.na(sig_stouffer$p.stouffer)),]
    stats::p.adjust( sig_stouffer$p.stouffer, method="BH" ,nrow(G)) -> sig_stouffer$q.stouffer

    if(filter == TRUE){
      sig_stouffer <- sig_stouffer[which(sig_stouffer$q.stouffer <= sigThreshold),]
    }

    sig_stouffer <- sig_stouffer[order(sig_stouffer$q.stouffer,decreasing = F),]
    if(output_file){
      utils::write.csv(sig_stouffer,file = "stouffer_significant_genes.csv",row.names = F)
    }

    result_stouffer <- subset(sig_stouffer,q.stouffer<0.1)
    if(nrow(result_stouffer)>30){
      result_stouffer1 <- result_stouffer[1:30,]
    }else{
      result_stouffer1 <- result_stouffer
    }
    result_stouffer1$lgq.stouffer <- -log10(result_stouffer1$q.stouffer+0.00000001)
    stouffer <- ggplot2::ggplot(result_stouffer1,ggplot2::aes(x=Hugo_Symbol,y=lgq.stouffer))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.stouffer))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by Stouffer's method", x ="Genes", y = "-log(q.stouffer)")
    if(output_file){
      ggplot2::ggsave(stouffer, file="stouffer_significant_genes.pdf", width = 20)
    }

    # logit
    sig_logit <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = 3))
    colnames(sig_logit) <- c("Hugo_Symbol","p.logit","q.logit")
    sig_logit$Hugo_Symbol <- G$Hugo_Symbol
    sig_logit$p.logit <- NA
    logit_binom <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))

    for(g in analysis_gene){
      g <- which(G$Hugo_Symbol %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_logit[g,]$p.logit <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & (N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num]) & BMR[g,1:cat_num]>0)
        for(i in flag){
          logit_binom[g,i] <- stats::binom.test(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],alternative = "greater")$p.value
        }
        tt <- (- sum(log(logit_binom[g,flag] / (1 - logit_binom[g,flag])))) / sqrt(cat_num * pi^2 * (5 * cat_num + 2) / (3 * (5 * cat_num + 4)))
        sig_logit[g,]$p.logit <- stats::pt(tt,df=5*cat_num+4, lower.tail=FALSE)
      }
    }
    sig_logit <- sig_logit[which(!is.na(sig_logit$p.logit)),]
    stats::p.adjust( sig_logit$p.logit, method="BH" ,nrow(G)) -> sig_logit$q.logit

    if(filter == TRUE){
      sig_logit <- sig_logit[which(sig_logit$q.logit <= sigThreshold),]
    }

    sig_logit <- sig_logit[order(sig_logit$q.logit,decreasing = F),]
    if(output_file){
      utils::write.csv(sig_logit,file = "logit_significant_genes.csv",row.names = F)
    }

    result_logit <- subset(sig_logit,q.logit<0.1)
    if(nrow(result_logit)>30){
      result_logit1 <- result_logit[1:30,]
    }else{
      result_logit1 <- result_logit
    }
    result_logit1$lgq.logit <- -log10(result_logit1$q.logit+0.00000001)
    logit <- ggplot2::ggplot(result_logit1,ggplot2::aes(x=Hugo_Symbol,y=lgq.logit))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.logit))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by logit method", x ="Genes", y = "-log(q.logit)")
    if(output_file){
      ggplot2::ggsave(logit, file="logit_significant_genes.pdf", width = 20)
    }

    # Brown's method
    sig_brown <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = 3))
    colnames(sig_brown) <- c("Hugo_Symbol","p.brown","q.brown")
    sig_brown$Hugo_Symbol <- G$Hugo_Symbol
    sig_brown$p.brown <- NA
    brown_binom <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))
    covariate <- system.file("extdata", "preprocessed_covariate.txt", package = "DriverGenePathway")
    covariate <- data.table::fread(covariate)

    for(gene in analysis_gene){
      if(!is.na(match(gene,covariate$Hugo_Symbol))){
        covariate_g <- covariate[match(gene,covariate$Hugo_Symbol),-1]
        g <- which(G$Hugo_Symbol %in% gene)
        if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
          sig_brown[g,]$p.brown <- 1
        }else{
          flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & (N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num]) & BMR[g,1:cat_num]>0)
          for(i in flag){
            brown_binom[g,i] <- stats::binom.test(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],alternative = "greater")$p.value
          }

          sig_brown[g,]$p.brown <- EmpiricalBrownsMethod::empiricalBrownsMethod(covariate_g, p_values = as.numeric(brown_binom[g,flag]), extra_info = FALSE)
        }
      }
    }
    sig_brown <- sig_brown[which(!is.na(sig_brown$p.brown)),]
    stats::p.adjust( sig_brown$p.brown, method="BH" ,nrow(G)) -> sig_brown$q.brown

    if(filter == TRUE){
      sig_brown <- sig_brown[which(sig_brown$q.brown <= sigThreshold),]
    }

    sig_brown <- sig_brown[order(sig_brown$q.brown,decreasing = F),]
    if(output_file){
      utils::write.csv(sig_brown,file = "brown_significant_genes.csv",row.names = F)
    }

    result_brown <- subset(sig_brown,q.brown<0.1)
    if(nrow(result_brown)>30){
      result_brown1 <- result_brown[1:30,]
    }else{
      result_brown1 <- result_brown
    }
    result_brown1$lgq.brown <- -log10(result_brown1$q.brown+0.00000001)
    brown <- ggplot2::ggplot(result_brown1,ggplot2::aes(x=Hugo_Symbol,y=lgq.brown))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.brown))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by brown method", x ="Genes", y = "-log(q.brown)")
    if(output_file){
      ggplot2::ggsave(brown, file="brown_significant_genes.pdf", width = 20)
    }

    # Kost's method
    sig_kost <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = 3))
    colnames(sig_kost) <- c("Hugo_Symbol","p.kost","q.kost")
    sig_kost$Hugo_Symbol <- G$Hugo_Symbol
    sig_kost$p.kost <- NA
    kost_binom <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))
    covariate <- system.file("extdata", "preprocessed_covariate.txt", package = "DriverGenePathway")
    covariate <- data.table::fread(covariate)

    for(gene in analysis_gene){
      if(!is.na(match(gene,covariate$Hugo_Symbol))){
        covariate_g <- covariate[match(gene,covariate$Hugo_Symbol),-1]
        g <- which(G$Hugo_Symbol %in% gene)
        if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
          sig_kost[g,]$p.kost <- 1
        }else{
          flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & (N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num]) & BMR[g,1:cat_num]>0)
          for(i in flag){
            kost_binom[g,i] <- stats::binom.test(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],alternative = "greater")$p.value
          }
          covariate_g <- rbind(covariate_g,covariate_g)
          # for (i in 1:(flag-1)) {
          #   covariate_g1 <- rbind(covariate_g1,covariate_g)
          # }
          sig_kost[g,]$p.kost <- EmpiricalBrownsMethod::kostsMethod(as.matrix(covariate_g), p_values = as.numeric(kost_binom[g,flag]), extra_info = FALSE)
        }
      }
    }
    sig_kost <- sig_kost[which(!is.na(sig_kost$p.kost)),]
    stats::p.adjust( sig_kost$p.kost, method="BH" ,nrow(G)) -> sig_kost$q.kost

    if(filter == TRUE){
      sig_kost <- sig_kost[which(sig_kost$q.kost <= sigThreshold),]
    }

    sig_kost <- sig_kost[order(sig_kost$q.kost,decreasing = F),]
    if(output_file){
      utils::write.csv(sig_kost,file = "kost_significant_genes.csv",row.names = F)
    }

    result_kost <- subset(sig_kost,q.kost<0.1)
    if(nrow(result_kost)>30){
      result_kost1 <- result_kost[1:30,]
    }else{
      result_kost1 <- result_kost
    }
    result_kost1$lgq.kost <- -log10(result_kost1$q.kost+0.00000001)
    kost <- ggplot2::ggplot(result_kost1,ggplot2::aes(x=Hugo_Symbol,y=lgq.kost))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.kost))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by kost method", x ="Genes", y = "-log(q.kost)")
    if(output_file){
      ggplot2::ggsave(kost, file="kost_significant_genes.pdf", width = 20)
    }

    # Harmonic mean p-value
    sig_harmonic <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = 3))
    colnames(sig_harmonic) <- c("Hugo_Symbol","p.harmonic","q.harmonic")
    sig_harmonic$Hugo_Symbol <- G$Hugo_Symbol
    sig_harmonic$p.harmonic <- NA
    harmonic_binom <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))
    covariate <- system.file("extdata", "preprocessed_covariate.txt", package = "DriverGenePathway")
    covariate <- data.table::fread(covariate)

    for(gene in analysis_gene){
      if(!is.na(match(gene,covariate$Hugo_Symbol))){
        covariate_g <- covariate[match(gene,covariate$Hugo_Symbol),-1]
        g <- which(G$Hugo_Symbol %in% gene)
        if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
          sig_harmonic[g,]$p.harmonic <- 1
        }else{
          flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & (N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num]) & BMR[g,1:cat_num]>0)
          for(i in flag){
            harmonic_binom[g,i] <- stats::binom.test(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],alternative = "greater")$p.value
          }
          covariate_g <- rbind(covariate_g,covariate_g)
          # for (i in 1:(flag-1)) {
          #   covariate_g1 <- rbind(covariate_g1,covariate_g)
          # }
          sig_harmonic[g,]$p.harmonic <- harmonicmeanp::p.hmp(as.numeric(harmonic_binom[g,flag]), L=length(flag))

        }
      }
    }
    sig_harmonic <- sig_harmonic[which(!is.na(sig_harmonic$p.harmonic)),]
    stats::p.adjust( sig_harmonic$p.harmonic, method="BH" ,nrow(G)) -> sig_harmonic$q.harmonic

    if(filter == TRUE){
      sig_harmonic <- sig_harmonic[which(sig_harmonic$q.harmonic <= sigThreshold),]
    }

    sig_harmonic <- sig_harmonic[order(sig_harmonic$q.harmonic,decreasing = F),]
    if(output_file){
      utils::write.csv(sig_harmonic,file = "harmonic_significant_genes.csv",row.names = F)
    }

    result_harmonic <- subset(sig_harmonic,q.harmonic<0.1)
    if(nrow(result_harmonic)>30){
      result_harmonic1 <- result_harmonic[1:30,]
    }else{
      result_harmonic1 <- result_harmonic
    }
    result_harmonic1$lgq.harmonic <- -log10(result_harmonic1$q.harmonic+0.00000001)
    harmonic <- ggplot2::ggplot(result_harmonic1,ggplot2::aes(x=Hugo_Symbol,y=lgq.harmonic))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.harmonic))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by harmonic method", x ="Genes", y = "-log(q.harmonic)")
    if(output_file){
      ggplot2::ggsave(harmonic, file="harmonic_significant_genes.pdf", width = 20)
    }


    # Likelihood ratio test
    sig_lrt <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = 3))
    colnames(sig_lrt) <- c("Hugo_Symbol","p.lrt","q.lrt")
    sig_lrt$Hugo_Symbol <- G$Hugo_Symbol
    sig_lrt$p.lrt <- NA
    # Likelihood ratio test
    LRT_lh0 <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))
    LRT_lh1 <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))

    for(g in analysis_gene){
      g <- which(G$Hugo_Symbol %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_lrt[g,]$p.lrt <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num] & BMR[g,1:cat_num]>0)
        for(i in flag){
          LRT_lh0[g,i] <- stats::dbinom(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],log = T)
          BMR_raw_n <- n_nonsilent_gc[g,i] + n_silent_gc[g,i] + n_noncoding_gc[g,i]
          BMR_raw_N <- N_nonsilent_gc[g,i] + N_silent_gc[g,i] + N_noncoding_gc[g,i]
          LRT_lh1[g,i] <- stats::dbinom(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR_raw_n/BMR_raw_N,log = T)
        }
        q = 2 * ( sum( LRT_lh1[g,flag] ) - sum( LRT_lh0[g,flag]))
        df = sum( LRT_lh1[g,flag] != 0 );
        if( df > 0 ) p_lrt = 1 - stats::pchisq( q, df )
        if( df == 0 ) p_lrt = 1
        sig_lrt$p.lrt[g] <- p_lrt
      }
    }
    sig_lrt <- sig_lrt[which(!is.na(sig_lrt$p.lrt)),]
    stats::p.adjust(sig_lrt$p.lrt, method="BH" , nrow(G)) -> sig_lrt$q.lrt
    if(filter == TRUE){
      sig_lrt <- sig_lrt[which(sig_lrt$q.lrt <= sigThreshold),]
    }
    sig_lrt <- sig_lrt[order(sig_lrt$q.lrt,decreasing = F),]
    if(output_file){
      utils::write.csv(sig_lrt,file = "likelihood_ratio_significant_genes.csv",row.names = F)
    }

    result_LRT <- subset(sig_lrt,q.lrt<0.1)
    if(nrow(result_LRT)>30){
      result_LRT1 <- result_LRT[1:30,]
    }else{
      result_LRT1 <- result_LRT
    }
    result_LRT1$lgq.LRT <- -log10(result_LRT1$q.lrt+0.00000001)
    LRT <- ggplot2::ggplot(result_LRT1,ggplot2::aes(x=Hugo_Symbol,y=lgq.LRT))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.LRT))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by likelihood ratio", x ="Genes", y = "-log(q.lrt)")
    if(output_file){
      ggplot2::ggsave(LRT, file="likelihood_ratio_significant_genes.pdf", width = 20)
    }


    # Convolution test
    sig_ct <- as.data.frame(matrix(data = NA,nrow = nrow(G),ncol = 3))
    colnames(sig_ct) <- c("Hugo_Symbol","p.ct","q.ct")
    sig_ct$Hugo_Symbol <- G$Hugo_Symbol
    sig_ct$p.ct <- NA

    LRT_lh0 <- as.data.frame(matrix(data = NA,nrow = nrow(G),ncol = cat_num))

    for(g in analysis_gene){
      g <- which(G$Hugo_Symbol %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_ct[g,]$p.ct <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num] & BMR[g,1:cat_num]>0)
        for(i in flag){
          LRT_lh0[g,i] <- stats::dbinom(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],log = T)
        }
      }
    }

    for(g in analysis_gene){
      g <- which(G$Hugo_Symbol %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_ct[g,]$p.ct <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num] & BMR[g,1:cat_num]>0)
        for(i in flag){
          ni <- N_nonsilent_gc[g,i]
          ei <- BMR[g,i]
          xmax = 100
          hmax = 25
          bin = 0.001
          gethist(xmax,ni,ei,ptype = "positive_log") -> bi
          binit(bi,hmax, bin ) -> bi
          if( i == flag[1] ) { hist0 = bi; }
          if( i > flag[1] & i < flag[length(flag)] ) { hist0 = convolute_b( hist0, bi ); binit( hist0, hmax, bin ) -> hist0 }
          if( i ==  flag[length(flag)]) { hist0 = convolute_b( hist0, bi ) }
        }
        # Convolution test
        bx = -sum( LRT_lh0[g,flag] )
        p.ct = sum( exp( -hist0[hist0>=bx] ))
        qc = sum( exp( -hist0 ))
        sig_ct$p.ct[g] <- p.ct
      }
    }
    sig_ct <- sig_ct[which(!is.na(sig_ct$p.ct)),]
    stats::p.adjust( sig_ct$p.ct, method="BH" ,nrow(G)) -> sig_ct$q.ct

    if(filter == TRUE){
      sig_ct <- sig_ct[which(sig_ct$q.ct <= sigThreshold),]
    }

    sig_ct <- sig_ct[order(sig_ct$q.ct,decreasing = F),]
    if(output_file){
      utils::write.csv(sig_lrt,file = "convolution_significant_genes.csv",row.names = F)
    }

    result_CT <- subset(sig_ct,q.ct<0.1)
    if(nrow(result_CT)>30){
      result_CT1 <- result_CT[1:30,]
    }else{
      result_CT1 <- result_CT
    }
    result_CT1$lgq.CT <- -log10(result_CT1$q.ct+0.00000001)
    CT <- ggplot2::ggplot(result_CT1,ggplot2::aes(x=Hugo_Symbol,y=lgq.CT))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.CT))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by convolution", x ="Genes", y = "-log(q.ct)")
    if(output_file){
      ggplot2::ggsave(CT, file="convolution_significant_genes.pdf", width = 20)
    }

    #PROJECTION
    N_nonsilent=bmr$N_nonsilent
    n_nonsilent=bmr$n_nonsilent
    X_gcp=bmr$X_gcp
    x_gcp=bmr$x_gcp
    null_categ=bmr$null_categ
    if(!quiet){cat("Calculating p.value using 2D Projection method...\n")}
    null_score_boost = 3
    min_effect_size =1.25
    convolution_numbins =1000

    ncat <- dim(X_gcp)[2]-1
    np <- dim(X_gcp)[3]

    G$p <-NaN
    for (g in 1:length(analysis_gene)){
      if (g %% 1000 == 0) print(g)

      g <- which(G$Hugo_Symbol %in% analysis_gene[g])
      P0 <- matrix(0,nrow=np,ncol = ncat)
      P1 <- matrix(0,nrow=np,ncol = ncat)
      priority <- matrix(0,nrow=np,ncol = ncat)
      N <- t(N_nonsilent[g,1:ncat,])
      n <- t(n_nonsilent[g,1:ncat,])
      x <- t(x_gcp[g,1:ncat,])
      X <- t(X_gcp[g,1:ncat,])

      P0 <- fun_hp(0,N,x,X)
      P1 <- fun_hp(1,N,x,X)
      priority <- t(apply(P1,1,order,decreasing=T))

      shft <- priority - matrix(rep(1:ncat,np),nrow=np,ncol=ncat,byrow = TRUE)
      map <- matrix(1:(np*ncat),nrow = np,ncol = ncat,byrow = FALSE)
      newmap <- map +shft*np
      P0 <- matrix(P0[newmap[1:length(newmap)]],nrow = np)
      P1 <- matrix(P1[newmap[1:length(newmap)]],nrow = np)
      P2 <- 1-(P0+P1)
      P2[P2<0] <- 0

      Pdeg <- array(0,c(np,ncat+1,ncat+1))

      for (d1 in 0:ncat){
        for (d2 in 0:d1)
        {
          if (d1 == ncat)
            p <- matrix(1,nrow = np,ncol = 1)
          else
            p <- apply(as.matrix(P0[,(d1+1):ncol(P0)]),1,prod)
          if (d1>0)
          { if (d1==d2)
            p = p*P2[,d1]
          else
          {
            p = p * P1[,d1]
            if (d2<(d1-1))
              p = p * apply(as.matrix(P0[,(d2+1):(d1-1)]),1,prod)
            if (d2>0)
              p = p * (P1[,d2]+P2[,d2])
          }
          }
          Pdeg[,d1+1,d2+1] <- p
        }
      }


      Sdeg = array(0,c(np,ncat+1,ncat+1))
      for (d1 in 1:ncat)
      {
        for (d2 in 0:d1)
        {
          if (d1==d2){
            p = P2[,d1]
          }else
          {
            if (d2>0)
              p = P1[,d1]*P1[,d2]
            else
              p = P1[,d1]
          }

          Sdeg[,d1+1,d2+1] = -log10(p)
        }
      }

      priority2 <- cbind(matrix(0,nrow = np,ncol = 1),priority)

      #Sdeg[priority2==null_categ] = Sdeg[priority2==null_categ] + null_score_boost
      s1 <- Sdeg[,,1]
      Sdeg[,,1][priority2==null_categ] <- s1[priority2==null_categ]+null_score_boost
      degree <-  matrix(0,nrow=np,ncol=2)
      score_obs = 0

      for(p in  1:np)
      {
        i = 1
        for (d in ncat:1)
        {
          c = priority[p,d]
          if (i==1)
          {
            if (n[p,c]>=2)
            {
              degree[p,] = c(d,d)
              i = 3
            }
            else if (n[p,c]==1)
            {
              degree[p,i] = d
              i=i+1
            }
          }
          else if (i==2)
          {
            if (n[p,c]>=1)
            {
              degree[p,i] = d
              i=i+1
            }
            else
              break
          }
        }
        score_sample = Sdeg[p,degree[p,1]+1,degree[p,2]+1]
        if(!is.finite(score_sample)){
          #   cat(sprintf("g=%d",g))
          #   cat(sprintf("p=%d",p))
          flag_finite <- is.finite(Sdeg[,degree[p,1]+1,degree[p,2]+1])
          score_sample <- max(Sdeg[,degree[p,1]+1,degree[p,2]+1][flag_finite])
          #   cat(sprintf("maxscore=%f",score_sample))
        }
        score_obs = score_obs + score_sample
      }

      score_obs = score_obs / min_effect_size

      if (score_obs<=0)
      {
        G$p[g]=1
        next
      }
      numbins = convolution_numbins
      binsize = score_obs / numbins
      H = matrix(0,nrow=numbins,ncol=1)
      H[1] <- 1

      offset = pmin(array(numbins,c(np,ncat+1,ncat+1)), round(Sdeg/binsize))
      ncols = (ncat+1)*(ncat+2)/2
      newH = matrix(0,nrow=numbins,ncol=ncols)
      for (p in 1:np)
      {
        newH[,] <- 0
        col=1
        for (d1 in 0:ncat)
          for (d2 in 0:d1)
          {
            o = offset[p,d1+1,d2+1]
            if (o<length(H))
              newH[(o+1):nrow(newH),col] = Pdeg[p,d1+1,d2+1] * H[1:(length(H)-o)]
            col=col+1
          }
        H = apply(newH,1,sum,drop=FALSE)
      }

      G$p[g] = max(0,1-sum(H))

    }  # of for (g in 1:ng)
    nG <- nrow(G)
    G <- G[which(!is.na(G$p)),]
    G$q <- stats::p.adjust(G$p,method="BH",nG)
    G <- G[order(G$q),]
    sig_projection <- subset(G,select=c("Hugo_Symbol","p","q"))

    if(filter == TRUE){
      sig_projection <- sig_projection[which(sig_projection$q <= sigThreshold),]
    }

    colnames(sig_projection) <- c("Hugo_Symbol","p.projection","q.projection")
    sig_projection <- sig_projection[which(!is.na(sig_projection$p.projection)),]
    if(output_file){
      utils::write.csv(sig_projection, file = "2d_projection_significant_genes.csv",row.names = F)
    }

    result_PJ <- subset(sig_projection,q.projection<0.1)
    if(nrow(result_PJ)>30){
      result_PJ1 <- result_PJ[1:30,]
    }else{
      result_PJ1 <- result_PJ
    }
    gene_projection <- result_PJ[,1]
    result_PJ1$lgq.PJ <- -log10(result_PJ1$q.projection+0.00000001)
    PJ <- ggplot2::ggplot(result_PJ1,ggplot2::aes(x=Hugo_Symbol,y=lgq.PJ))+
      ggplot2::geom_point(size = 5, ggplot2::aes(color = lgq.PJ))+
      ggplot2::scale_color_gradient(low="blue", high="red")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.background = ggplot2::element_rect(fill="white"))+
      ggplot2::labs(title="Plot of result genes tested by 2D projection", x ="Genes", y = "-log(q.projection)")
    if(output_file){
      ggplot2::ggsave(PJ, file="2d_projection_significant_genes.pdf", width = 20)
    }

    denovo_result <- driverPathway(maf = bmr_maf, threshold = 5, analyse_genes = analysis_gene,
                                   driver_size=3,
                                   pop_size=200,
                                   iters=200,
                                   permut_time=500,
                                   output_file = output_file,
                                   quiet = quiet)
    if(!is.matrix(denovo_result$driver_geneset)){
      denovo_result_genes <- denovo_result$driver_geneset
    }else{
      denovo_result_genes <- c()
      for (i in nrow(denovo_result$driver_geneset)) {
        denovo_result_genes <- unique(append(denovo_result_genes, denovo_result$driver_geneset[i,]))
      }
    }

    sigGenes_combined <- c(
      as.character(sig_betabinomial$Hugo_Symbol[sig_betabinomial$q.btBinom <= sigThreshold]),
      as.character(sig_fisher$Hugo_Symbol[sig_fisher$q.fisher <= sigThreshold]),
      as.character(sig_stouffer$Hugo_Symbol[sig_stouffer$q.stouffer <= sigThreshold]),
      as.character(sig_logit$Hugo_Symbol[sig_logit$q.logit <= sigThreshold]),
      as.character(sig_brown$Hugo_Symbol[sig_brown$q.brown <= sigThreshold]),
      as.character(sig_kost$Hugo_Symbol[sig_kost$q.kost <= sigThreshold]),
      as.character(sig_harmonic$Hugo_Symbol[sig_harmonic$q.harmonic <= sigThreshold]),
      as.character(sig_ct$Hugo_Symbol[sig_ct$q.ct <= sigThreshold]),
      as.character(sig_lrt$Hugo_Symbol[sig_lrt$q.lrt <= sigThreshold]),
      as.character(sig_projection$Hugo_Symbol[sig_projection$q.projection <= sigThreshold]),
      as.character(denovo_result))
    sigGenes_table <- as.data.frame(table(sigGenes_combined))
    sigGenes_table <- sigGenes_table[order(sigGenes_table$Freq,decreasing = T),]

    sigGenes_final <- sigGenes_table[which(sigGenes_table$Freq == max(sigGenes_table$Freq)),1]

    q_result <- sig_betabinomial
    gene_betabinomial <- q_result[,1]

    q_result <- subset(q_result,q.btBinom<=0.1,select=-p.btBinom)

    q_new <- sig_fisher
    gene_fisher <- q_new[,1]

    q_new <- subset(q_new,q.fisher<=0.1)
    q_union <- data.frame(Hugo_Symbol=union(q_new$Hugo_Symbol,q_result$Hugo_Symbol))
    flag_old <- match(q_union$Hugo_Symbol,q_result$Hugo_Symbol)
    flag_new <- match(q_union$Hugo_Symbol,q_new$Hugo_Symbol)
    q_union[,2] <- q_result[flag_old,2]
    q_union$q.fisher <- q_new[flag_new,3]
    q_result <- q_union

    q_new <- sig_stouffer
    gene_fisher <- q_new[,1]

    q_new <- subset(q_new,q.stouffer<=0.1)
    q_union <- data.frame(Hugo_Symbol=union(q_new$Hugo_Symbol,q_result$Hugo_Symbol))
    flag_old <- match(q_union$Hugo_Symbol,q_result$Hugo_Symbol)
    flag_new <- match(q_union$Hugo_Symbol,q_new$Hugo_Symbol)
    q_union[,2:3] <- q_result[flag_old,2:3]
    q_union$q.stouffer <- q_new[flag_new,3]
    q_result <- q_union

    q_new <- sig_logit
    gene_fisher <- q_new[,1]

    q_new <- subset(q_new,q.logit<=0.1)
    q_union <- data.frame(Hugo_Symbol=union(q_new$Hugo_Symbol,q_result$Hugo_Symbol))
    flag_old <- match(q_union$Hugo_Symbol,q_result$Hugo_Symbol)
    flag_new <- match(q_union$Hugo_Symbol,q_new$Hugo_Symbol)
    q_union[,2:4] <- q_result[flag_old,2:4]
    q_union$q.logit <- q_new[flag_new,3]
    q_result <- q_union

    q_new <- sig_brown
    gene_fisher <- q_new[,1]

    q_new <- subset(q_new,q.brown<=0.1)
    q_union <- data.frame(Hugo_Symbol=union(q_new$Hugo_Symbol,q_result$Hugo_Symbol))
    flag_old <- match(q_union$Hugo_Symbol,q_result$Hugo_Symbol)
    flag_new <- match(q_union$Hugo_Symbol,q_new$Hugo_Symbol)
    q_union[,2:5] <- q_result[flag_old,2:5]
    q_union$q.brown <- q_new[flag_new,3]
    q_result <- q_union

    q_new <- sig_kost
    gene_fisher <- q_new[,1]

    q_new <- subset(q_new,q.kost<=0.1)
    q_union <- data.frame(Hugo_Symbol=union(q_new$Hugo_Symbol,q_result$Hugo_Symbol))
    flag_old <- match(q_union$Hugo_Symbol,q_result$Hugo_Symbol)
    flag_new <- match(q_union$Hugo_Symbol,q_new$Hugo_Symbol)
    q_union[,2:6] <- q_result[flag_old,2:6]
    q_union$q.kost <- q_new[flag_new,3]
    q_result <- q_union

    q_new <- sig_harmonic
    gene_fisher <- q_new[,1]

    q_new <- subset(q_new,q.harmonic<=0.1)
    q_union <- data.frame(Hugo_Symbol=union(q_new$Hugo_Symbol,q_result$Hugo_Symbol))
    flag_old <- match(q_union$Hugo_Symbol,q_result$Hugo_Symbol)
    flag_new <- match(q_union$Hugo_Symbol,q_new$Hugo_Symbol)
    q_union[,2:7] <- q_result[flag_old,2:7]
    q_union$q.harmonic <- q_new[flag_new,3]
    q_result <- q_union

    q_new <- sig_lrt
    gene_lrt <- q_new[,1]

    q_new <- subset(q_new,q.lrt<=0.1)
    q_union <- data.frame(Hugo_Symbol=union(q_new$Hugo_Symbol,q_result$Hugo_Symbol))
    flag_old <- match(q_union$Hugo_Symbol,q_result$Hugo_Symbol)
    flag_new <- match(q_union$Hugo_Symbol,q_new$Hugo_Symbol)
    q_union[,2:8] <- q_result[flag_old,2:8]
    q_union$q.lrt <- q_new[flag_new,3]
    q_result <- q_union

    q_new <- sig_ct
    gene_ct <- q_new[,1]

    q_new <- subset(q_new,q.ct<=0.1)
    q_union <- data.frame(Hugo_Symbol=union(q_new$Hugo_Symbol,q_result$Hugo_Symbol))
    flag_old <- match(q_union$Hugo_Symbol,q_result$Hugo_Symbol)
    flag_new <- match(q_union$Hugo_Symbol,q_new$Hugo_Symbol)
    q_union[,2:9] <- q_result[flag_old,2:9]
    q_union$q.ct <- q_new[flag_new,3]
    q_result <- q_union

    q_new <- sig_projection
    gene_projection <- q_new[,1]

    q_new <- subset(q_new,q.projection<=0.1)
    q_union <- data.frame(Hugo_Symbol=union(q_new$Hugo_Symbol,q_result$Hugo_Symbol))
    flag_old <- match(q_union$Hugo_Symbol,q_result$Hugo_Symbol)
    flag_new <- match(q_union$Hugo_Symbol,q_new$Hugo_Symbol)
    q_union[,2:10] <- q_result[flag_old,2:10]
    q_union$q.projection <- q_new[flag_new,3]
    q_result <- q_union

    names(q_result) <- c("Hugo_Symbol","q.betaBinomial","q.fisherBinomial","q.stouffer","q.logit",
                         "q.brown","q.kost","q.harmonic","q.likelihoodRatio",
                         "q.convolution","q.projection")

    flag <- c()
    for (i in 1:nrow(q_result)) {
      if(sum(is.na(q_result[i,]))<3){
        flag <- append(flag,i)
      }
    }
    q_result_sub <- q_result[flag,]

    q_result1 <- q_result_sub[,c(1,2)]
    q_result_names <- colnames(q_result_sub)
    for (i in 3:ncol(q_result_sub)) {
      colnames(q_result1) <- q_result_names[c(1,i)]
      q_result1 <- rbind(q_result1, q_result_sub[,c(1,i)])
    }
    colnames(q_result1) <- c("Hugo_Symbol","q")
    flag <- !is.na(q_result1$q)
    q_result1 <- q_result1[flag,]

    q_box <- ggplot2::ggplot(q_result1, ggplot2::aes(x = Hugo_Symbol, y = q))+
      ggplot2::geom_boxplot(ggplot2::aes(color = Hugo_Symbol))+
      ggplot2::theme_minimal()+
      ggplot2::labs(title="Box plot of q values tested by all methods", x = "Genes", y = "q_values")


    q_result2 <- q_result_sub[,c(1,2)]
    q_result_names <- colnames(q_result_sub)
    for (i in 3:ncol(q_result_sub)) {
      colnames(q_result2) <- q_result_names[c(1,i)]
      q_result2 <- rbind(q_result2, q_result_sub[,c(1,i)])
    }
    colnames(q_result2) <- c("Hugo_Symbol","q")
    q_result2$method <- c(rep("q.betaBinomial",nrow(q_result_sub)),rep("q.fisherBinomial",nrow(q_result_sub)),rep("q.stouffer",nrow(q_result_sub)),
                          rep("q.logit",nrow(q_result_sub)),rep("q.brown",nrow(q_result_sub)),rep("q.kost",nrow(q_result_sub)),
                          rep("q.harmonic",nrow(q_result_sub)),rep("q.likelihoodRatio",nrow(q_result_sub)),rep("q.convolution",nrow(q_result_sub)),
                          rep("q.projection",nrow(q_result_sub)))
    flag <- !is.na(q_result2$q)
    q_result2 <- q_result2[flag,]
    q_result2$repel <- NA
    for (i in unique(q_result2$Hugo_Symbol)) {
      flag <- q_result2$Hugo_Symbol==i
      j <- sample(1:sum(flag),1)
      q_result2$repel[which(flag)[j]] <- q_result2$Hugo_Symbol[which(flag)[j]]
    }

    all_results <- ggplot2::ggplot(q_result2, ggplot2::aes(x=method, y=-log(q), group=Hugo_Symbol))+
      ggplot2::geom_line(ggplot2::aes(color = Hugo_Symbol),linetype = "dotted") +
      ggplot2::geom_point(ggplot2::aes(color = Hugo_Symbol))+
      ggrepel::geom_text_repel(ggplot2::aes(label = repel), max.overlaps = getOption("ggrepel.max.overlaps", default = 10), size = 3.5)+
      ggplot2::theme_minimal()+
      ggplot2::theme(axis.text.x = element_text(angle=25))+
      ggplot2::labs(title="Significant genes tested by all methods", x = "Genes", y = "-log(q)")

    if(output_file){
      utils::write.csv(q_result,file = "significant_gene_q_values.csv",quote = F,row.names = F)
      ggplot2::ggsave(q_box, file="box_q_values.pdf", width = 9, height = 4)
      ggplot2::ggsave(all_results, file="genes_by_all.pdf", width = 9, height = 4)
    }


    return(as.character(sigGenes_final))

  }
  if(output_file){
    setwd("..")
  }

} # of MutSig_runCV function



