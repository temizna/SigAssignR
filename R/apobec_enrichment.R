#' Calculate APOBEC enrichment
#'
#' @description Calculate APOBEC enrichment
#' @param inputdata a dataframe with columns of twentycontext (20 nucleotide context),
#'  con (trinucleotide context), sample (Sample names), mut_TCW (TCW context mutations), 
#'  MutC (mutations at C) This dataframe is created using either read_nonvcf function or
#'  vcfs_to_inputdata function
#' @return enrichment table with adjusted p-values
#' @export
apobec_enrichment<-function(inputdata){
inputdata$ConTCW=str_count(inputdata$twentycontext,"TCA") + str_count(inputdata$twentycontext,"TCT")+
  str_count(inputdata$twentycontext,"TGA") + str_count(inputdata$twentycontext,"TGT")
inputdata$ConC=str_count(inputdata$twentycontext,"C") + str_count(inputdata$twentycontext,"G")
#Mark TCW and C mutation
tc1<-which((as.character(inputdata$alt)=="A"))
tc21<-which(inputdata$context=="TCA")
tc22<-which(inputdata$context=="TCT")
inputdata$mut_TCW=0
inputdata$mut_TCW[tc21]=1
inputdata$mut_TCW[tc22]=1
inputdata$mut_TCW[tc1]=0
inputdata$MutC=0
c1<-which(as.character(inputdata$ref)=="C")
c12<-which((as.character(inputdata$alt)=="A"))
inputdata$MutC[c1]=1
inputdata$MutC[c12]=0
#count the TCW and C contexts and mutations to create the contingency table

ta<-aggregate(as.numeric(inputdata$mut_TCW),by=list(inputdata$sample), FUN="sum")
ta1<-aggregate(as.numeric(inputdata$mut_TCW),by=list(inputdata$sample), FUN="sum")$x
ta2<-aggregate(as.numeric(inputdata$MutC),by=list(inputdata$sample), FUN="sum")$x
ta3<-aggregate(inputdata$ConTCW,by=list(inputdata$sample), FUN="sum")$x
ta4<-aggregate(inputdata$ConC,by=list(inputdata$sample), FUN="sum")$x

enrich.table=data.frame(mutTCW=ta1, ConTCW=ta3,mutC=ta2, ConC=ta4, row.names = ta$Group.1)
enrich.table$score=(enrich.table$mutTCW/enrich.table$ConTCW)/(enrich.table$mutC/enrich.table$ConC)
enrich.table$MutRatio=enrich.table$mutTCW/enrich.table$ConTCW
enrich.table$ConRatio=enrich.table$mutC/enrich.table$ConC
enrich.table$mutC_TCW=enrich.table$mutC-enrich.table$mutTCW
enrich.table$conC_TCW=enrich.table$ConC-enrich.table$ConTCW

entest<-as.data.frame(enrich.table[,c(1,2,8,9)])
entest2<-as.data.frame(enrich.table[,1:4])


fishers <- t(as.data.frame(apply(entest, 1, exe_fisher)))
fishers <- as.data.frame(fishers)

fishers2 <- t(as.data.frame(apply(entest2, 1, exe_fisher)))
fishers2 <- as.data.frame(fishers2)
fishers2$BH<-p.adjust(fishers2$V1, method="BH")
fishers2$Bon<-p.adjust(fishers2$V1, method="bonferroni")

enrich.table$fisher_pval <- fishers$V1
enrich.table$bh_adj_qval <- p.adjust(enrich.table$fisher_pval, method = "BH")
enrich.table$Bon <- p.adjust(enrich.table$fisher_pval, method = "bonferroni")
#enrich.table=enrich.table[sort(enrich.table[,11]),]
return(enrich.table)
}
