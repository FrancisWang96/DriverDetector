#' Assign categories for maf and coverage
#'
#' The mutation categories generated from 192 basic types to designated number of types.
#' The reference chromosome is necessary if the number of categories is set more than one.
#' The chromosome can either be name of a folder containing txt files or name of the installed BSgenome package.
#' Example: BSgenome.Hsapiens.UCSC.hg19.
#' The preprocess of maf and coverage is integrated in this function.
#'
#' @param maf Mutation Annotation Format, can be whether an R data frame or the path of a txt file.
#' @param coverage Coverage file, can be whether an R data frame or the path of a txt file.
#' @param ref_genome Reference chromosome, either a folder name or a installed BSgenome package.
#' Default NULL, tries to auto-detect from installed genomes.
#' @param category_num Number of mutation categories, default 4. If the number is 0, categories should exist.
#' @param output_file Determine whether to export the categories, categorised maf and coverage as txt files.
#' @param quiet Whether to show notes during processing.
#' @return Categorised maf and coverage.
#' @examples
#' laml_maf <- system.file("extdata", "tcga_laml.maf", package = "DriverGenePathway")
#' coverage <- system.file("extdata", "coverage.rda", package = "DriverGenePathway")
#' load(coverage)
#' Categ_result <- mafCategoryAssign(laml_maf, coverage, category_num = 1, output_file = TRUE,
#' quiet = FALSE)
#' @export

mafCategoryAssign <- function(maf = NULL, coverage = NULL, ref_genome = NULL, category_num = NULL,
                              output_file = TRUE, quiet = FALSE){

  Hugo_Symbol = NULL
  categ = NULL
  categ_idx = NULL
  effect = NULL
  pkgname = NULL
  # if not preprocessed, preprocess
  # check if maf is to preprocess
  if(is.character(maf)){
    if(!quiet){cat('Loading maf ... ')}
    maf <- data.table::fread(maf)
    if(!quiet){cat('success!\n')}
  }

  if(!"effect" %in% names(maf)){
    if(!quiet){cat('Preprocessing maf ... ')}
    maf <- mafPreprocess(maf, output_file = FALSE, quiet = TRUE)
    if(!quiet){cat('success!\n')}
  }

  # check if coverage is to preprocess

  if(is.character(coverage)){
    if(!quiet){cat('Loading coverage ... ')}
    coverage <- data.table::fread(coverage)
    if(!quiet){cat('success!\n')}
  }

  if(!quiet){cat('Preprocessing coverage ... ')}
  coverage <- coveragePreprocess(coverage, maf, output_file = FALSE, quiet = TRUE)
  if(!quiet){cat('success!\n')}

  # categ match between maf and coverage, if so categs_already_present is set TRUE.
  ucc <- unique(coverage$categ)
  bases <- c("A","C","G","T")
  names192 <- data.frame(names=0)
  i=1
  for (from in 1:4) {
    for (left in 1:4) {
      for (right in 1:4) {
        for (to in 1:4) {
          if(from==to){
            next
          }
          names192[i,1] <- paste(bases[left],"(",bases[from],"->",bases[to],")",
                                 bases[right],sep = "")
          i <- i+1
        }
      }
    }
  }

  categs_already_present <- FALSE

  if("categ" %in% colnames(maf)){
    mcc <- unique(maf$categ)
    ucc <- unique(coverage$categ)
    if(length(ucc) != length(mcc) || sum(!(ucc %in% mcc)) > 0 ||
       sum(!(mcc %in% ucc)) > 0){
      if(!quiet){cat("The current category in maf is not available, a new one will be generated.\n")}
      maf <- subset(maf,select = -categ)
    }else{
      categs_already_present <- TRUE
    }
  }

  # find out if we can do category discovery
  can_do_category_discovery <- TRUE
  coverage_is_on_full192 <- (length(ucc) == 192 &
                             length(intersect(ucc,as.matrix(names192))) == 192)
  if(!coverage_is_on_full192){
    can_do_category_discovery <- FALSE
    reason <- "coverage does not cover all 192 triplets"
  }
  # 如果都没有呢？
  if(!("chr" %in% tolower(colnames(maf))) & "chromosome" %in%
     tolower(colnames(maf))){
    colnum <- which(tolower(colnames(maf)) %in% "chromosome")
    maf$chr <- maf[[colnum]]
  }else if("chr" %in% tolower(colnames(maf))){
    colnum <- which(tolower(colnames(maf)) %in% "chr")
    maf$chr <- maf[[colnum]]
  }else{
    can_do_category_discovery <- FALSE
    reason <- "Chromosome/chr column missing in maf"
  }

  if(stringr::str_detect(tolower(maf$chr[1]), 'chr')){
    maf$chr <- stringr::str_sub(maf$chr,4,-1)
  }

  if(!("start" %in% tolower(colnames(maf))) & "start_position" %in%
     tolower(colnames(maf))){
    colnum <- which(tolower(colnames(maf)) %in% "start_position")
    maf$start <- maf[[colnum]]
  }else if("start" %in% tolower(colnames(maf))){
    colnum <- which(tolower(colnames(maf)) %in% "start")
    maf$start <- maf[[colnum]]
  }else{
    can_do_category_discovery <- FALSE
    reason <- "Start_position/start column missing in maf"
  }

  # 优先级tumor_seq_allele2 > tumor_seq_allele1 > reference_allele
  if(!("ref_allele" %in% colnames(maf)) & "reference_allele" %in%
     tolower(colnames(maf))){
    colnum <- which(tolower(colnames(maf)) %in% "reference_allele")
    maf$ref_allele <- maf[[colnum]]
  }else if("ref_allele" %in% tolower(colnames(maf))){
    colnum <- which(tolower(colnames(maf)) %in% "ref_allele")
    maf$ref_allele <- maf[[colnum]]
  }else{
    can_do_category_discovery <- FALSE
    reason <- "Reference_allele/ref_allele column missing in maf"
  }

  if(!("newbase" %in% tolower(colnames(maf)))){
    if("tumor_seq_allele1" %in% tolower(colnames(maf))){
      colnum <- which(tolower(colnames(maf)) %in% "tumor_seq_allele1")
      maf$newbase <- maf[[colnum]]
      if("tumor_seq_allele2" %in% tolower(colnames(maf))){
        colnum_allele1 <- which(tolower(colnames(maf)) %in% "tumor_seq_allele1")
        idx <- which(maf$ref_allele == maf[[colnum_allele1]])
        colnum_allele2 <- which(tolower(colnames(maf)) %in% "tumor_seq_allele2")
        maf$newbase[idx] <- maf[[colnum_allele2]][idx]
      }
    }
  }else if("newbase" %in% tolower(colnames(maf))){
    colnum <- which(tolower(colnames(maf)) %in% "newbase")
    maf$newbase <- maf[[colnum]]
  }else{
    can_do_category_discovery <- FALSE
    reason <- "Reference_allele&tumor_seq_allele1 & tumor_seq_allele2 column missing in maf"
  }

  # decide method
  if(is.null(category_num)){
    if(can_do_category_discovery){
      method <- 3
      ncategs <- 4
    }else if(categs_already_present){
      method <- 1
    }else{
      cat(sprintf("NOTE: unable to perform category discovery, because %s. ",
                  reason))
      method<- 2
    }
  }else if(category_num == 0){
    if(!categs_already_present){
      stop("When setting category_num==0, categ column must be already present in
           maf.")
    }
    method <- 1
  }else if(category_num == 1){
    method <-2
  }else if(category_num > 1){
    if(category_num>6){
      cat("NOTE: maximum categories that can be discovered is 6. \n")
      category_num <- 6
    }
    if(!can_do_category_discovery){
      stop(sprintf("unable to perform category discovery, because %s", reason))
    }else{
      method <- 3
      ncategs <- category_num
    }
  }

  f <- colnames(coverage)
  coverage_patient_names <- f[4:ncol(coverage)]

  # assign category
  if(method==1){
    # Using the categories already presented.
    K <- data.frame(names = rep(1,each=length(unique(maf$categ))))
    K$names <- sort(unique(maf$categ))
  }else if(method==2){
    # Will use two categories: missense and null+indel.

    K <- dplyr::tibble(left = c('ACGT','ACGT'),
                    from = c('AC','AC'),
                    change = c('in','in'),
                    right = c('ACGT','ACGT'),
                    autoname = c('missense','null+indel'),
                    name = c('missense','null+indel'),
                    type = c('point','non-point'))

    # assign categories
    maf$categ <- "---"
    flag_null <- which(maf$effect == "null")
    maf$categ[flag_null] <- K$name[2]
    flag_nonull <- which(maf$effect != "null")
    maf$categ[flag_nonull] <- K$name[1]

    # collapse coverage
    temp <- unique(coverage$categ)
    coverage$categ_idx <- match(coverage$categ,temp)

    coverage <- dplyr::arrange(coverage,Hugo_Symbol,effect,categ_idx)
    #order as Hugo_Symbol, effect, categ_idx, if decrease, then desc(Hugo_Symbol)

    ug <- unique(coverage$Hugo_Symbol)
    ng <- length(ug)
    ue <- unique(coverage$effect)
    ne <- length(ue)
    nk <- nrow(K)
    # 为什么只取前两列？因为coverage会去掉！哈哈哈
    idx <- which(coverage$categ_idx <= nk)
    C2 <- coverage[idx,]
    C2$categ[which(C2$categ_idx == 1)] <- K$name[1]
    C2$categ[which(C2$categ_idx == 2)] <- K$name[2]
    C2 <- subset(C2,select = c(Hugo_Symbol,effect,categ))
    np <- length(coverage_patient_names)

    for (p in 1:np) {
      oldcov <- array(as.integer(unlist(coverage[[coverage_patient_names[p]]])),
                      c(192,ne,ng))
      newcov <- rep(colSums(oldcov),each=2)
      C2[[coverage_patient_names[p]]] <- newcov
    }
    coverage <- C2
    remove(C2)
  }else if(method == 3){

    bad <- which(is.nan(maf$start))
    if(length(bad) > 0){
      if(!quiet){
        cat(sprintf("WARNING: %d/%d mutations have non-numeric Start_positions and are
                  excluded from analysis. \n",
                    length(bad),length(maf$start)))
      }
      maf <- maf[-bad]
    }
    if(nrow(maf) == 0){
      stop("There is no available data left!")
    }

    if(length(list.files(ref_genome))){
      situation <- 'folder'
    }else{
      hsgs.installed = BSgenome::installed.genomes(splitNameParts = TRUE)
      data.table::setDT(x = hsgs.installed)

      if(is.null(ref_genome)){
        if(nrow(hsgs.installed) == 0){
          stop("Could not find any installed BSgenomes.\n Use BSgenome::available.genomes() for options.")
        }else{
          if(!quiet){
            cat("-Found following BSgenome installtions. Using first entry\n")
          }
          print(hsgs.installed)
          situation <- 'bsg'
          ref_genome = hsgs.installed[,pkgname][1]
        }
      }else{
        if(sum(hsgs.installed[,pkgname] %in% ref_genome) == 0){
          if(!quiet){
            cat(paste0("-Could not find BSgenome "), ref_genome, "\n")
          }
          if(nrow(hsgs.installed) == 0){
            stop("Could not find any installed BSgenomes either.\nUse BSgenome::available.genomes() for options.")
          }else{
            if(!quiet){
              cat("-Found following BSgenome installtions. Correct ref_genome argument if necessary\n")
            }
            print(hsgs.installed)
            stop()
          }
        }else{
          situation <- 'bsg'
        }
      }
      ref_genome = BSgenome::getBSgenome(genome = ref_genome)
    }

    if(situation == 'bsg'){
      f1_flag <- grep("^chr",names(ref_genome))
      f1 <- names(ref_genome)[f1_flag]
      uchr <- unique(maf$chr)
      uchr <- uchr[order(uchr)]
      match_table <- data.frame(uchr = uchr, chr = 1:length(uchr))
      maf$chr_idx <- match(maf$chr,uchr)
      # 把"chr"去掉，并转为字符
      uchr <- gsub("chr","",uchr)
      f2 <- paste("chr",uchr, sep = "")
      chr_file_available <- f2 %in% f1
      if(sum(chr_file_available) == 0){
        stop('ref_genome mismatch with maf')
      }else{
        bad <- which(!chr_file_available[maf$chr_idx])
        if(length(bad) > 0){
          if(!quiet){
            cat(sprintf("WARNING: %d/%d mutations cannot match with ref_genome and are excluded from analysis. \n",
                        length(bad),nrow(maf)))
          }
          maf <- maf[-bad]
        }
        if(nrow(maf) == 0){
          stop("There is no available data left!")
        }
      }
    }else{
      f1_flag <- grep("^chr.*txt$",list.files(ref_genome))
      f1 <- list.files(ref_genome)[f1_flag]
      uchr <- unique(maf$chr)
      uchr <- uchr[order(uchr)]
      maf$chr_idx <- match(maf$chr,uchr)
      # 把"chr"去掉，并转为字符
      uchr <- gsub("chr","",uchr)
      f2 <- paste("chr",uchr,".txt",sep = "")
      chr_file_available <- f2 %in% f1
      if(sum(chr_file_available) == 0){
        stop('no chr files available')
      }else{
        bad <- which(!chr_file_available[maf$chr_idx])
        if(length(bad) > 0){
          if(!quiet){
            cat(sprintf("WARNING: %d/%d mutations cannot match with ref_genome and are excluded from analysis. \n",
                        length(bad),nrow(maf)))
          }
          maf <- maf[-bad]
        }
        if(nrow(maf) == 0){
          stop("There is no available data left!")
        }
      }
    }



    # cat(sprintf("Generating mutation categories...\n"))
    C_coverage <- dplyr::as_tibble(coverage[,4:ncol(coverage)])
    if(length(C_coverage) == 1){
      coverage$totcov <- C_coverage$value
    }else{
      coverage$totcov <- apply(C_coverage,1,sum)
    }

    # coverage*病人数
    npm <- length(unique(maf$Tumor_Sample_Barcode))
    if((length(coverage_patient_names)==1) & (npm > 1)){
      coverage$totcov <- coverage$totcov * npm
    }

    # will use only the coding mutations+coverage to do this
    coverage$is_coding[which((coverage$effect == "nonsilent")|(coverage$effect == "silent"))] <- 1
    coverage$is_coding[which(!((coverage$effect == "nonsilent")|(coverage$effect == "silent")))] <- 0

    # collapse coverage to 192
    X <- data.frame(categ=rep(1,each=length(unique(coverage$categ))))
    X$categ <- unique(coverage$categ)
    coverage$categ_idx <- match(coverage$categ,unique(coverage$categ),nomatch = 0)
    X$left <- X$categ
    X$split <- strsplit(X$categ,split = "")
    X$left <- substr(X$categ,start=1,stop = 1)
    X$from <- substr(X$categ,start = 3,stop = 3)
    X$to <- substr(X$categ,start = 6,stop = 6)
    X$right <- substr(X$categ,start = 8,stop = 8)

    # 对coding区的192个categ加起来
    X$yname <- paste(X$from,"in",X$left,sep = " ")
    X$yname <- paste(X$yname,"_",X$right,sep = "")
    C_iscoding <- coverage[which(coverage$is_coding == 1),]
    X$N <- tapply(C_iscoding$totcov,C_iscoding$categ_idx,sum)
    X$newbase_idx <-  match(X$to,c("A","C","G","T"))

    # 生成64个"A in A_A"还有一个"any N"
    x <- c("A","C","G","T")
    x_in <- rep(x,each=16)
    x_left <- rep(rep(x,each=4),times=4)
    x_right <- rep(x,times=16)
    y <- paste(x_in,"in",x_left,sep = " ")
    y <- paste(y,"_",x_right,sep = "")
    y[65] <- "any N"
    Y_name <- y
    rm(x,x_in,x_left,x_right,y)
    Y <- data.frame(num=1:65,name=Y_name)
    X$context65 <- match(X$yname,Y_name,nomatch = 65)
    # 为了Y的N和X的N一样，但是没必要啊
    Y$N <- round(rbind(as.matrix(tapply(X$N,X$context65,sum)),0)/3)

    # STEPE 2
    #MUTATIONS：get context65 by looking up from reference genome
    if(!quiet){cat('Searching trinucleotide using chromosome information... \n')}

    #f2 <- paste("chr",uchr,".txt",sep = "")
    if(situation == 'folder'){
      f2 <- paste(ref_genome,f2,sep = "/")
      triplet <- data.frame(name = rep(1,each=nrow(maf)))
      for (ci in 1:sum(chr_file_available)){
        if(!quiet){
          cat(sprintf("%d/%d ", ci,sum(chr_file_available)))
        }
        midx <- which(maf$chr_idx == ci)
        chrfile <- f2[ci]
        d <- file.info(chrfile)
        if(is.na(d$size)){
          next
        }
        # filesize <- d$size
        ff <- data.table::fread(f2[ci],header = F)
        #find mutation base and its left and right neighbours
        triplet$name[midx] <- stringr::str_sub(ff[1],maf$start[midx]-1,maf$start[midx]+1)
      }
    }else{
      triplet <- data.frame(name = rep(1, each=nrow(maf)))
      for (ci in 1:sum(chr_file_available)){
        if(!quiet){
          cat(sprintf("%d/%d ", ci,sum(chr_file_available)))
        }
        midx <- which(maf$chr_idx == ci)
        genome_match <- which(stringr::str_c('chr', match_table$uchr[ci])==ref_genome@user_seqnames)
        ff <- ref_genome[[genome_match]]
        if(is.na(ff@length)){
          next
        }
        #find mutation base and its left and right neighbours
        triplet$name[midx] <- stringr::str_sub(ff,maf$start[midx]-1,maf$start[midx]+1)
      }
    }

    flag_non_triplet <- which(triplet == 1)
    if(length(flag_non_triplet)){
      triplet$name[flag_non_triplet] <- "---"
    }
    maf$triplet <- toupper(triplet$name)
    maf$triplet_middle <- stringr::str_sub(maf$triplet,2, 2)
    # 信息全的
    midx <- which((maf$ref_allele != "-") & (maf$newbase != "-") &
                    (!is.na(maf$ref_allele)) & (!is.na(maf$newbase)) & (!is.na(maf$triplet_middle)))
    # maf中的参考位点与chr中的位点
    matchfrac <- sum(maf$ref_allele[midx] == maf$triplet_middle[midx])/length(midx)
    #matchfrac
    if(matchfrac < 0.9){
      adj <- "possible"
      if(matchfrac < 0.7){
        adj <- "probable"
      }
      if(!quiet){
        cat(sprintf("\n NOTE: %s build mismatch between mutation_file and chr_files", adj))
      }
    }
    maf$yname <- paste(substr(maf$triplet,2,2),"in",substr(maf$triplet,1,1),sep = " ")
    maf$yname <- paste(maf$yname,substr(maf$triplet,3,3),sep = "_")
    maf$context65 <- match(maf$yname,Y$name,nomatch = 65)
    maf$newbase_idx <- match(substr(maf$newbase,1,1),c("A","C","G","T"),
                           nomatch = NaN)

    midx <- which((maf$ref_allele != "-") & (maf$newbase != "-") &
                    (maf$context65 >= 1) & (maf$context65 <= 65) &
                    (maf$newbase_idx >= 1) & (maf$newbase_idx <= 4))
    M_Na <- subset(maf,select = c("context65","newbase_idx"))[midx,]
    M_N <- plyr::count(M_Na, names(M_Na))
    # 没用到
    bases <- c("A","C","G","T")
    Y$A <- 0
    Y$C <- 0
    Y$G <- 0
    Y$T <- 0

    for (i in 1:4) {
      M_N_i <- M_N[which(M_N$newbase_idx == i),]
      Y[[bases[i]]][M_N_i$context65] <- M_N_i$freq
    }

    #STEP3 Category Discovery
    #Category Discovery

    # 把G in A_C to T 与C in G_T to A归位一类
    Nn <- preprocessCollapseNn65to32(Y)
    P <- data.frame(max_k = ncategs, mutcategs_report_filename =
                      "mafCategoryAssign result.txt")

    Ks <- preprocessFindMutCateg(Nn,P)
    # 找完了
    K <- Ks[[ncategs]]
    c <- preprocessAssignCateg(K)
    X$kidx <- matrix(data = NaN,nrow = nrow(X),ncol = 1)
    for (i in 1:nrow(X)) {
      X$kidx[i] <- which(c[X$context65[i],,X$newbase[i]]==1)
    }

    #STEP4
    #assign mutation categories
    if(!quiet){
      cat('\n')
      cat(sprintf("Assigning mutation categories..."))

    }
    maf$categ <- rep("---",times=nrow(maf))
    for (i in 1:nrow(X)) {
      idx <- which((maf$context65==X$context65[i]) &
                     (maf$newbase_idx==X$newbase_idx[i]))
      maf$categ[idx] <- rep(as.character(K$name[X$kidx[i]]),each=length(idx))
    }

    # add null+indel category
    K <- rbind(K,dplyr::tibble(left = c('ACGT'),
                               from = c('AC'),
                               change = c('in'),
                               right = c('ACGT'),
                               N=sum(K$N[1:(nrow(K))]),
                               n=length(which(maf$effect == "null")),
                               rate=length(which(maf$effect == "null"))/sum(K$N[1:(nrow(K))]),
                               relrate=length(which(maf$effect == "null"))/sum(K$N[1:(nrow(K))])/K$rate[1]*K$relrate[1],
                               autoname = c('null+indel'),
                               name = c('null+indel'),
                               type = c('non-point')))

    midx <- which(maf$effect == "null")
    maf$categ[midx] <- 'null+indel'
    if(!quiet){cat('success!\n')}

    # STEP5
    # collapse coverage
    if(!quiet){cat(sprintf("Collapsing coverages..."))}
    coverage <- dplyr::arrange(coverage,Hugo_Symbol,effect,categ_idx)
    #order as Hugo_Symbol, effect, categ_idx, if decrease, then desc(Hugo_Symbol)

    ug <- unique(coverage$Hugo_Symbol)
    ng <- length(ug)
    ue <- unique(coverage$effect)
    ne <- length(ue)
    nk <- nrow(K)

    idx <- which(coverage$categ_idx <= nk)

    K_name <- unlist(unique(K$name))
    C2 <- coverage[idx,]
    # find categ of Coverage
    #C2_categ <- rep(NaN,times=nrow(C2))
    idx1 <- which((!is.nan(C2$categ_idx))&(C2$categ_idx>=1) &
                    (C2$categ_idx<=length(K_name)))
    C2$categ <- K_name[C2$categ_idx[idx1]]
    C2 <- subset(C2,select = c(Hugo_Symbol,effect,categ))

    np <- length(coverage_patient_names)

    for (p in 1:np) {
      oldcov <- array(as.integer(unlist(coverage[[coverage_patient_names[p]]])),
                      c(192,ne,ng))
      newcov <- array(data = NaN,c(nk,ne,ng))
      for (ki in 1:nk) {
        if(ki==nk){
          cidx <- c(1:192)
        }else{
          cidx <- which(X$kidx==ki)
        }
        newcov[ki,,] <- apply(oldcov[cidx,,], c(2,3), sum)
      }

      C2[[coverage_patient_names[p]]] <- as.vector(newcov)
    }

    coverage <- C2
    remove(C2)
    if(!quiet){cat('success!\n')}
  }

  if(!quiet){cat('Mutation categories:\n')}
  if(method==1 | method==2){
    K1 <- K
  }else{
    rownames(K) <- 1:nrow(K)
    K1 <- dplyr::select(K,c('name','autoname','n','N','rate','type'))
    names(K1) <- c('categ_name','category','n','N','rate','type')
  }
  if(!quiet){print(K1)}
  if(!quiet){cat('\n')}

  if(output_file){
    if(!dir.exists('DriverGenePathway_output')){
      dir.create('DriverGenePathway_output')
    }
    setwd('DriverGenePathway_output')
    utils::write.table(maf,"categorised_maf.txt", sep = "\t", quote = F, row.names = F)
    utils::write.table(coverage,"categorised_coverage.txt", sep = "\t", quote = F, row.names = F)
    utils::write.table(as.matrix(K1),"categories.txt", sep = "\t", quote = F, row.names = F)
    setwd("..")
  }

  if(!quiet){cat('Preprocess of category assign finishes.')}

  return(list(maf=maf,coverage=coverage))
}
