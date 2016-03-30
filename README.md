# Drug repositioning for psoriasis based on cogena

Here is a knitr-powered report, including code, to reproduce all the results in 
the manuscript, **Drug repositioning and drug mode of action discovery based 
on co-expressed gene-set enrichment analysis**.


* __Drug repositioning for psoriasis based on cogena.pdf__ shows all the code and results related with GSE13355.
* __GSE30999.pdf__ shows all the code and results related with GSE30999.
* __src__ contains the Rmd files to generate the two pdf files above, the WGCNA-pathway analysis of the two datasets (GSE13355 and GSE30999).
* __data__ is used to store the raw data from GEO.
* __results__ contains the the input and output of _GSEA_, _CMap_,  _NFFinder_ and _WGCNA_.
* __tmp__ is to store the temporary file.

cogena package is available at [Bioconductor](https://bioconductor.org/packages/cogena/) or [my Github](https://github.com/zhilongjia/cogena).
