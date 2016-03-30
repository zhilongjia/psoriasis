


enrichTable <- function(enrichmentTable, CutoffNumGeneset=Inf, CutoffPVal=0.05,
                        orderMethod="max", roundvalue=TRUE) {
    
    # Gene numbers in each clusters, I, II and All.
    NumGeneInCluster <- as.vector( table(gene_cluster_info ) )
    score <- enrichmentTable
    
    # the orderMethod options
    orderMethod <- match.arg(orderMethod, c(rownames(score), "max", "mean"))
    if (orderMethod == "mean") {
        score = score[,order(colMeans(score, na.rm=TRUE), decreasing=TRUE)]
        index_above_cutoffPVal <- which(suppressWarnings(
            apply(score, 2, max, na.rm=TRUE)) > -log2(CutoffPVal))
    } else if (orderMethod == "max") {
        colMax <- function(X) {suppressWarnings(
            apply(X, 2, max, na.rm=TRUE))}
        score = score[,order(colMax(score), decreasing=TRUE)]
        index_above_cutoffPVal <- which(suppressWarnings(
            apply(score, 2, max, na.rm=TRUE)) > -log2(CutoffPVal))
    } else if (orderMethod %in% rownames(score)) {
        score = score[, order(score[orderMethod,], decreasing=TRUE)]
        index_above_cutoffPVal <- 
            which(score[orderMethod,] > -log2(CutoffPVal))
    }
    
    if (length(index_above_cutoffPVal) > CutoffNumGeneset){
        score <- score[,c(1:CutoffNumGeneset)]
    } else if (length(index_above_cutoffPVal) == 0){
        score <- NA
        return (score)
    } else {
        #drop para used as length(index_above_cutoffPVal)==1.
        score <- score[,index_above_cutoffPVal, drop=FALSE]
    }
    
    score[score<=-log2(CutoffPVal)] <- NA
    
    # drop para used as length(index_above_cutoffPVal)==1.
    score <- score[,ncol(score):1, drop=FALSE]
    
    colnames(score) <- tolower(strtrim(colnames(score), 60))
    
    rownames(score) <- paste(rownames(score), as.character(NumGeneInCluster), sep="#")
    
    if (roundvalue){
        score <- round(score,1)
    }
    
    return (score)
    
    
}