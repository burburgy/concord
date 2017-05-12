
coxen <- function( set1, set2, probeset=probe.cor, cutoff=0.2, method="spearman"){
	set1.cor.cor = (cor(t(set1[probeset, ]))) ## correlation among genes within NCI60 drug cell line
	set2.cor.cor = ( cor(t(set2[probeset, ])) )   ## correlation among genes within training set
	n.cor.num <- length(probeset)
	coxen.cor.cor = matrix(1, nrow=n.cor.num, ncol=2)
	for (i in 1:n.cor.num) {
	    #print(i)
	    coxen.cor.cor[i,1] = cor(set1.cor.cor[i,-i], set2.cor.cor[i,-i], method=method)
	    coxen.cor.cor[i,2] = cor.test(set1.cor.cor[i,-i], set2.cor.cor[i,-i], alternative="greater", method=method)$p.value
	}
	result <- list()
	qvalue <- p.adjust(coxen.cor.cor[,2], method="fdr")
	result$probe <- probeset[ which( qvalue < cutoff ) ]
      rownames(coxen.cor.cor) <- probeset
	result$coxen.cor.cor <- coxen.cor.cor
	result$coxen <- coxen.cor.cor[ result$probe, ]
	return(result)
}

