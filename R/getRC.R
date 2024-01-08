#' Get the reverse complement of a DNA string
#'
#' @description et the reverse complement of a DNA string
#' @param s DNA sequence to be converted
#' @return reverse complement sequence string
#' @export
getRC<-function(s){
  out<-character()
  for(i in 1:nchar(s[1])){
    nt<-substr(s,i,i)
    out<-paste(ifelse(nt=="G","C",ifelse(nt=="C","G",ifelse(nt=="A","T","A"))),out,sep="")
  }
  return(out)
}
#' Execute one-sided Fisher's exact test
#'
#' @description Execute one-sided Fisher's exact test
#' @param x data frame with elements of contingency table as columns and samples as rows
#' @return Fisher's exact p-calue
#' @export
exe_fisher <- function(x) {
  m <- matrix(unlist(x), ncol = 2, nrow = 2, byrow = F)
  f <- fisher.test(m, alternative = "greater")
  return(as.data.frame(f$p.value))
}
