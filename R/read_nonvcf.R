#' Read a non-vcf mutation file to create the inputdata frame.
#'
#'
#' @param inputdata a dataframe with columns chr, start, stop, ref, alt, sample
#' @param bsg reference BSGenome object to map
#' @param ext extention of context (current only 1 is supported)
#'
#' @return dataframe with columns displaying 20 nucleotide context, 3 nucleotide context
#' @export
read_nonvcf<-function(inputdata, bsg=bsg, ext=1){
  #default bsg = BSgenome.Hsapiens.UCSC.hg19
  cols=c("chr","start","stop","ref","alt", "sample")
  #input file/table should be ordered as "chr", "start", "stop", "ref", "alt", "sample")
  colnames(inputdata)=c("chr","start","stop","ref","alt", "sample")
  chrs=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")
  allchr=which(inputdata$chr %in% chrs)
  #filter ref and alt sof single base substitutions and remove any indels and longer mutations
  inputdata=inputdata[allchr,]
  in1=which(nchar(inputdata$alt)==1)
  inputdata<-inputdata[in1,]
  in2<-which(!(inputdata$alt=="-"))
  inputdata<-inputdata[in2,]
  in1r=which(nchar(inputdata$ref)==1)
  inputdata<-inputdata[in1r,]
  in2r<-which(!(inputdata$ref=="-"))
  inputdata<-inputdata[in2r,]
  
  grlinputdata<-GRanges(seqnames=paste("chr",inputdata$chr,sep=""),
                   ranges=IRanges(start=inputdata$start-ext,end=inputdata$start+ext))
  grlinputdata20<-GRanges(seqnames=paste("chr",inputdata$chr,sep=""),
                     ranges=IRanges(start=inputdata$start-20,end=inputdata$start+20))
  #get 41 bases around each SBS for APOBEC enrichment calculation
  inputdata$twentycontext=BSgenome::getSeq(bsg,grlinputdata20,as.character=T)
  inputdata$context=BSgenome::getSeq(bsg,grlinputdata,as.character=T)
 
  #identify and return the reverse complement for A/G ref positions and their alts 
  temprca<-which(inputdata$ref =="A")         
  temprcg<-which(inputdata$ref == "G") 
  inputdata$ref[temprca]=getRC(inputdata$ref[temprca])
  inputdata$ref[temprcg]=getRC(inputdata$ref[temprcg])
  inputdata$alt[temprca]=getRC(inputdata$alt[temprca])
  inputdata$alt[temprcg]=getRC(inputdata$alt[temprcg])
  inputdata$twentycontext[temprca]=getRC(inputdata$twentycontext[temprca])
  inputdata$twentycontext[temprcg]=getRC(inputdata$twentycontext[temprcg])
  inputdata$context[temprca]=getRC(inputdata$context[temprca])
  inputdata$context[temprcg]=getRC(inputdata$context[temprcg])
  if(ext==1){
    inputdata$con=paste(substr(inputdata$context,1,1),"[",inputdata$ref,">",inputdata$alt,"]",substr(inputdata$context,3,3),sep="")
  } else {
    inputdata$con=paste(substr(inputdata$context,1,2),"[",inputdata$ref,">",inputdata$alt,"]",substr(inputdata$context,4,5),sep="")
  }
 
  inputdata<-unique(inputdata)
  return(inputdata)
}
