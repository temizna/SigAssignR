#' Find optimal nonnegative linear combination of mutation signatures to
#' reconstruct the mutation matrix using nnls package.
#'
#' @param mat mutation count matrix (dimensions: 96 mutation types
#' X n samples)
#' @param known_sig_mat Known Signature matrix (dimensions: 96 mutation types
#' 30 n signatures) 
#' @param cut_off Minimum signature contribution required
#'
#' @return Named list with signature contributions and reconstructed
#' mutation matrix
#'
#' @import nnls
#'
#' @examples
#'
#'
#' ## Get signatures
#' signatures <- get_known_signatures()
#'
#' ## Perform the fitting
#' fit <- assign_signatures(known_sig_mat, mat, cut_off=0.03)
#'

#' @export
assign_signature<-function(known_sig_mat, mat, cut_off=0.03){
  #if(is.na(known_sig_mat)) known_sig_mat=cancer_signatures
  if(dim(known_sig_mat)[1] != dim(mat)[1]){
    stop("Number of mutational channels are not equal")
  }
  assigned.coeff<-matrix(nrow=dim(mat)[2],ncol=dim(known_sig_mat)[2])
  assigned.fit<-matrix(nrow=dim(known_sig_mat)[1],ncol=dim(mat)[2])
  assigned.res<-matrix(nrow=dim(known_sig_mat)[1],ncol=dim(mat)[2])
  assigned.dev=vector(length=dim(mat)[2])
  names(assigned.dev)=colnames(mat)
  rownames(assigned.coeff)<-colnames(mat)
  colnames(assigned.fit)<-colnames(mat)
  colnames(assigned.res)<-colnames(mat)
  colnames(assigned.coeff)=colnames(known_sig_mat)
  rownames(assigned.fit)=rownames(known_sig_mat)
  rownames(assigned.res)=rownames(known_sig_mat)
  
  for (i in colnames(mat)){
    muts=sum(mat[,i])
    assign.sol=nnls(known_sig_mat,mat[,i])
    assign.sol.frac<-assign.sol$x/sum(assign.sol$x)
    sigtouse=which(assign.sol.frac>=cut_off)
    assigned.coeff[i,]=assign.sol$x
    if (length(sigtouse) >= 2){
      assign.sol=nnls(known_sig_mat[,sigtouse],mat[,i])
      assigned.coeff[i,]=0
      assigned.coeff[i,sigtouse]=assign.sol$x
    }
    assigned.fit[,i]=round(assign.sol$fitted[,1])
    assigned.res[,i]=round(assign.sol$residuals[,1])
    assigned.dev[i]=assign.sol$deviance
  }
  results=list(coefficient=floor(assigned.coeff), fitted=floor(assigned.fit),residuals=floor(assigned.res), deviance=assigned.dev)
  return(results)
}

