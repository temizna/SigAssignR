#' Read a vcf mutation file to create the inputdata frame.
#' Uses read_vcf function from sigminer package
#'
#'
#' @param vcffiles a vcf file list
#' @param bsg reference BSGenome object to map
#' @param ext extention of context (current only 1 is supported)
#' @param samples sample name list
#' @param genome one of "mm10" or "hg38" description
#'
#'
#' @return dataframe with columns displaying 20 nucleotide context, 3 nucleotide context
#' @export
vcfs_to_inputdata<-function(vcffiles, bsg=bsg, ext=1, samples=sample_names, genome=gen){
  mafs<-read_vcf(vcffiles,genome_build = gen)
  mafs2=as.data.frame(mafs@data)
  inputdata=mafs2[,c(2,3,6,4,5,1)]
  #default bsg = BSgenome.Hsapiens.UCSC.hg19
  cols=c("chr","start","stop","ref","alt", "sample")
  #input file/table should be ordered as "chr", "start", "stop", "ref", "alt", "sample")
  colnames(inputdata)=c("chr","start","stop","ref","alt", "sample")
  chrs=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
         "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","X","Y")
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
  
  grlinputdata<-GRanges(seqnames=inputdata$chr,
                        ranges=IRanges(start=inputdata$start-ext,end=inputdata$start+ext))
  grlinputdata20<-GRanges(seqnames=inputdata$chr,
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
