#' Search driver pathway using de novo method
#'
#' This function searches driver pathway using de novo method based on mutual exclusivity and coverage
#' It outputs the driver pathway and its p-value as a txt file.
#'
#' @param mutation_matrix Matrix of mutations occur in which patient and which gene.
#' @param maf Mutation Annotation Format, can be whether an R data frame or the path of a txt file.
#' @param threshold Minimum mutation count for genes to be analysed, default to be 5.
#' @param analyse_genes Genes to be analysed.
#' @param exclude_genes Genes excluded from analysis.
#' @param analyse_patients Patrients to be analysed.
#' @param exclude_patients Patrients excluded from analysis.
#' @param background_gene_filter Exclude genes having fewer mutations than expected using a binomial test.
#' @param coverage Coverage file, can be whether an R data frame or the path of a txt file.
#' @param covariate Covariate file, can be whether an R data frame or the path of a txt file.
#' @param original_mutation_rate The default background mutation rate is 1.2e-6, which is used to
#' screen a subset of genes through binomial test.
#' @param ref_genome Reference chromosome, either a folder name or a installed BSgenome package.
#' Default NULL, tries to auto-detect from installed genomes.
#' @param category_num Number of mutation categories, default 4. If the number is 0, categories should exist.
#' @param driver_size Size of output driver gene set, default to be 3.
#' @param pop_size Population size of GA, default to be 200.
#' @param iters Time of iteration of GA, default to be 500.
#' @param permut_time Time of Permutation test, default to be 1000.
#' @param output_file Determine whether to export the preprocessed maf as a txt file.
#' @param quiet Whether to show notes during processing.
#' @examples
#' mutation_matrix <- system.file("extdata", "sample_mutation_matrix.rda",
#' package = "DriverGenePathway")
#' load(mutation_matrix)
#' driverPathway(mutation_matrix, driver_size = 3, pop_size = 30, iters = 100, permut_time = 50)
#'
#' laml_maf <- system.file("extdata", "tcga_laml.maf", package = "DriverGenePathway")
#' driverPathway(maf = laml_maf, driver_size = 3, pop_size = 30, iters = 100, permut_time = 50)
#' @export


driverPathway <- function(mutation_matrix = NULL, maf = NULL, threshold = 5, analyse_genes = NULL,
                          exclude_genes = NULL, analyse_patients = NULL, exclude_patients = NULL,
                          background_gene_filter = FALSE, coverage = NULL, covariate = NULL,
                          original_mutation_rate = 1.2e-6, ref_genome = NULL, category_num = NULL,
                       driver_size=3,
                       pop_size=200,
                       iters=200,
                       permut_time=500,
                       output_file = TRUE,
                       quiet = FALSE) {

    if(!is.null(mutation_matrix)){
        if(!quiet){cat('Loading mutation matrix ... ')}
        if(is.character(mutation_matrix)){
            mutation_matrix <- utils::read.table(mutation_matrix)
        }
        if(!quiet){cat('success!\n')}
    }else if(!is.null(maf)){
        if(!quiet){cat('Transforming maf to mutation matrix ... ')}
        mutation_matrix <- mafToMatrix(maf, threshold, analyse_genes, exclude_genes,
                                analyse_patients, exclude_patients, output_file, quiet = TRUE)
        if(!quiet){cat('success!\n')}
    }else{
        stop("Mutaion matrix or maf file missing, cannot do analysis.")
    }

    if(background_gene_filter){
        if(!quiet){cat('Eliminating genes with fewer mutations than expected using a binomial test... ')}
        backgroundGenes <- backgroundBinomial(maf, coverage, covariate, original_mutation_rate,
                                              ref_genome, category_num)

        mutation_matrix <- mutation_matrix[, which(colnames(mutation_matrix) %in% backgroundGenes)]
        if(!quiet){cat('success!\n')}
    }

    # start with an random population
    total_gene_size <- ncol(mutation_matrix)
    population = matrix(nrow=(pop_size*2), ncol=total_gene_size);
    # fill values
    for (child in 1:(pop_size*2)) {
        population[child,] <- rep(0,total_gene_size)
        population[child,sample(1:total_gene_size,driver_size)] <- 1
    }

    # calculate each object
    evalVals = rep(NA, pop_size*2);
    for (object in 1:(pop_size*2)) {
        evalVals[object] = evalFunc(mutation_matrix,population[object,]);
    }

    order_flag <- order(evalVals,decreasing = F)
    population <- population[order_flag,]
    Record <- 0
    iter <- 1
    pj <- 1
    while( (iter < iters) & (Record < 10) ){
        weight_vector <- evalVals[order_flag]
        t_minweight <- min(weight_vector)
        for (pj in 1:pop_size) {
            selected_index <- select_operator(weight_vector)
            population[pop_size+pj,] <- cross_operator(population[selected_index[1],], population[selected_index[2],],total_gene_size)
            population[pop_size+pj,] <- mutation_operator(mutation_matrix,
                                                          population[pop_size+pj,],total_gene_size,1)
            weight_vector[pop_size+pj] <- evalFunc(mutation_matrix,
                                                   population[pop_size+pj,])
        }
        population <- population[order(weight_vector,decreasing = F),]
        weight_vector <- sort(weight_vector,decreasing = F)

        if(Record == 2){
            temp <- sort(unique(weight_vector),decreasing = F)
            index <- which(weight_vector <= temp[1])
            for (j in min(length(index),1)) {
                population[index[j],] <- mutation_operator(mutation_matrix,
                                                           population[index[j],],
                                                           total_gene_size,
                                                           floor(sqrt(total_gene_size)))
                weight_vector[index[j]] <- evalFunc(mutation_matrix,
                                                    population[index[j],])
            }

        }

        if(Record == 5){
            temp <- sort(unique(weight_vector),decreasing = F)
            index <- which(weight_vector <= temp[1])
            for (j in min(length(index),1)) {
                population[index[j],] <- mutation_operator(mutation_matrix,
                                                           population[index[j],],
                                                           total_gene_size,
                                                           total_gene_size)
                weight_vector[index[j]] <- evalFunc(mutation_matrix,
                                                    population[index[j],])
            }

        }


        min_weight <- min(weight_vector)

        if(min_weight == t_minweight){
            Record <- Record + 1
        }else{
            Record <- 0
        }
        iter <- iter + 1

    }

    order_flag <- order(weight_vector)
    population <- population[order_flag,]
    max_pop <- population[1:pop_size,]

    # delete the duplicate rows
    len_maxpop <- dim(max_pop)[1]
    max_pop <- max_pop[!duplicated(max_pop),]

    if(length(nrow(max_pop) > 1) > 0){
        driver_geneset <- matrix(data = NA, ncol = driver_size, nrow = nrow(max_pop))
        for (i in 1:nrow(max_pop)) {
            driver_geneset[i,] <- colnames(mutation_matrix)[which(max_pop[i,] == 1)]
        }
    }else{
        driver_geneset <- colnames(mutation_matrix)[which(max_pop==1)]

    }

    p_value <- mulExclusive_significance(mutation_matrix,driver_geneset,permut_time = permut_time)
    G <- list(driver_geneset=driver_geneset,p_value=p_value)
    utils::write.table(G,file = "denovoDriverPathway.txt",quote = F,sep = "\t",row.names = F)
    return(G)
}
