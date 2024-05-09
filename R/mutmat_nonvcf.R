#' Make mutation count matrix of 96 trinucleotides from a dataframe with contexts and samplenames 
#' for each mutation (this data frame is created using read_nonvcf)
#'
#' @description Make 96 trinucleotide mutation count matrix
#' @param inputdata Dataframe containing con (context column) and sample column to form the mutation matrix
#' @return 96 mutation count matrix
#' @export
mutmat_nonvcf<-function(inputdata){
  #create the mutation matrix from nonvcf inputdata (see read_nonvcf.R for creating the inputdata from a nonvcf table)
  mutorder=c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T", 
          "G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T",
          "A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T",
          "G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T",
          "A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T",
          "G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T",
          "A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T",
          "G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T",
          "A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T",
          "G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T",
          "A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T",
          "G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
  mut.table<-table(inputdata$con,inputdata$sample)
  missing=which(!(mutorder %in% rownames(mut.table)))
  missdata=as.table(matrix(0, nrow =length(missing), ncol=dim(mut.table)[2]))
  rownames(missdata)=mutorder[missing]
  colnames(missdata)=colnames(mut.table)
  mut.table2=rbind(mut.table, missdata)
  mut.table=mut.table2
  mut.table=mut.table[mutorder,]
  return(mut.table)
}
