driverPathway <- function(maf = NULL, threshold = 5,
                          analyse_genes = NULL,
                       driver_size=3,
                       pop_size=200,
                       iters=200,
                       permut_time=500,
                       output_file = TRUE,
                       quiet = FALSE) {

    # if(!is.null(mutation_matrix)){
    #     if(!quiet){cat('Loading mutation matrix ... ')}
    #     if(is.character(mutation_matrix)){
    #         mutation_matrix <- utils::read.table(mutation_matrix)
    #     }
    #     if(!quiet){cat('success!\n')}

  if(!quiet){cat('Transforming maf to mutation matrix ... ')}
  mutation_matrix <- mafToMatrix(maf, threshold, analyse_genes = NULL, output_file = output_file, quiet = TRUE)
  if(!quiet){cat('success!\n')}

  if(!quiet){cat('Eliminating genes with fewer mutations than expected using a binomial test... ')}
  mutation_matrix <- mutation_matrix[, which(colnames(mutation_matrix) %in% analyse_genes)]
  if(!quiet){cat('success!\n')}




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
    if(output_file){
      utils::write.table(G,file = "denovoDriverPathway.txt",quote = F,sep = "\t",row.names = F)
    }
    return(G)
}
