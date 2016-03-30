


heatmapTable <- function(enrichmentTable, 
         CutoffNumGeneset=20, CutoffPVal=0.05,
         orderMethod="max", roundvalue=TRUE,
         low="grey", high="red", na.value="white", 
         maintitle=NULL, printGS=FALSE) {

    enrichment_score <- enrichTable(enrichmentTable, CutoffNumGeneset, 
                                   CutoffPVal, orderMethod, roundvalue)
    
    if (length(enrichment_score)==1 && is.na(enrichment_score)){
        return(paste("No enrichment above the cutoff for", method, 
                     "when the number of clusters is", nCluster, 
                     "with the orderMethod:", orderMethod, "!"))
    }
    
    if (printGS==TRUE) {
        cat (rev(colnames(enrichment_score)), sep ="\t")
    }
    
    # cl_color <- sub("#.*", "", rownames(enrichment_score))
    
    enrichment_score <- reshape2::melt(enrichment_score)
    #legend breaks
    cutoff_score <- round(-log2(CutoffPVal), 2)
    max_score <- max(enrichment_score$value, na.rm=TRUE)
    if (is.infinite(max_score)) {max_score <- 1000}
    if (max_score/cutoff_score > 2) {
        breaks <- c(cutoff_score, round( seq(cutoff_score, max_score, length.out=5)[-1]))
    } else if (max_score/cutoff_score >1) {
        breaks <- c(cutoff_score, round( seq(cutoff_score, max_score, length.out=2)[-1]))
    } else {
        breaks <- NULL
    }
    
    Var1=Var2=value=NULL

    ggplot2::ggplot(enrichment_score, aes(as.factor(Var1), Var2)) + 
        geom_tile(aes(fill = value)) + 
        scale_fill_gradient2("score", space="Lab", mid=low, midpoint=4, low=low, 
                             high=high, na.value=na.value, breaks=breaks) +
        geom_text(aes(fill=value, label=value),size=4, na.rm=TRUE) +
        labs(list(title = maintitle, x = "Cluster", y = "Gene set")) +
        theme(axis.text.y = element_text(size = rel(1.5), face="bold")) +
        theme(axis.text.x = element_text(size = rel(1.3), angle=-90, 
                                         # color=cl_color, 
                                         face="bold", vjust=0.5)) +
        theme(panel.grid.major.x = element_line(color = "grey", size = 5),
              panel.grid.major.y = element_blank())
}