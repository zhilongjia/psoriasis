
raw_clEnrich <- function(gene_cluster_info, annofile=NULL, TermFreq=0){
    

    ############################################################################
    # Annotation data
    if (is.null(annofile)) {
        annofile <- system.file("extdata", "c2.cp.kegg.v5.0.symbols.gmt.xz", package="cogena")
    }
    annotation <- gene2set(annofile, names(gene_cluster_info), TermFreq=TermFreq)
    # the background gene gene-sets matrix
    AllGeneSymbols=NULL
    data(AllGeneSymbols, envir = environment())
    annotationGenesPop <- gene2set(annofile, AllGeneSymbols, TermFreq=TermFreq)
    annotationGenesPop <- annotationGenesPop[,colnames(annotation)]
    
    if (ncol(annotationGenesPop) ==0 || ncol(annotation)==0) {
        stop("Error in annotation as ncol equals zero. 
        Maybe lower the TermFreq value.")
    }

    ############################################################################
    # Gene sets enrichment analysis for clusters
    clen <- list()
    
    cluster <- gene_cluster_info
    
    pei <- matrix(NA, nrow=length(unique(cluster)), 
                  ncol=ncol(annotation))
    rownames(pei) <- sort(unique(cluster))
    colnames(pei) <- colnames(annotation)
    for (k in  sort(unique(cluster))) {
        genenames <- names(which(cluster==k))
        pei[as.character(k),] <- cogena::PEI(genenames, annotation=annotation, 
                                             annotationGenesPop=annotationGenesPop)
    }

    # negative log2 p value
    logAdjPEI <- function (pei) {
        # fdr based on pval
        pei.adjust <- matrix(p.adjust(pei, "fdr"), ncol=ncol(pei))
        dimnames(pei.adjust) <- dimnames(pei)
        pei.NeglogPval <- -log2(pei.adjust)
    }
    pei <- logAdjPEI(pei)
    
    # res <- list(gene_cluster_info=gene_cluster_info,
    #             enrichmentTable = pei)
    
    return (pei)
}


