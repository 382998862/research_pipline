#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript


# load library
library('getopt');


#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'help' , 'h', 0, "logical",
  'infile' , 'i', 1, "character",
  'sample' , 's', 1, "character",
  'outfile' , 'o', 1, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);


# define usage function
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
      Usage example:
      1) Rscript simpleBar.r --infile T01,circ,T02.circ --outfile out_simpleBar.circ
      Options:
      --help		NULL 		get this help
      --infile 	character 	the input file [forced]
      --sample  character   split by ","
      --outfile 	character 	the filename for output graph [forced]
      \n")
  q(status=1);
}



# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }


# check non-null args
if ( is.null(opt$infile) )	{ print_usage(spec) }
if ( is.null(opt$outfile) )	{ print_usage(spec) }


#-----------------------------------------------------------------

chooseGeneId <- function(x){
  x<-as.vector(as.matrix(x))
  x=paste(x,collapse=",")
  result <- unlist(strsplit(x,","))
  result <- result[which(result[] != "NA")]
  count <- length(result[which(result[] == "n/a")])
  if(count==length(result)){
    result = 'n/a'
  }else{
    result <- result[which(result[] != "n/a")]
    result <- unique(result)
    result=paste(result,collapse=",")
  }
  return (result)
}


chooseType <- function(x)
{
  x<-as.vector(as.matrix(x))
  if("exon"%in%x)
  {
    return("exon")
    
  }else if("intron"%in%x){
    return("intron")
  }else if("intergenic_region"%in%x){
    return("intergenic_region")
  } else{
    return("n/a")
  }
  
}

chooseStrand <- function(x)
{
  x<-as.vector(as.matrix(x))
  x=x[!is.na(x)]
  x <- unique(x)
  x=paste(x,collapse="")
  return(x)
}
# reading data
#-----------------------------------------------------------------
# reading data
sample<-unlist(strsplit(opt$infile,','))
key<-unlist(strsplit(opt$sample,','))
suffixes <- key
r1 <- read.table(file = sample[1], sep = "\t", header=TRUE)
colnames(r1) <-c("chr","start","end",suffixes[1],"circRNA_type","gene_id","strand")
tmp<-r1
count <-1
for(s in sample)
{

  if(count>1)
  {
    r2 <- read.table(file = s, sep = "\t", header=TRUE)
    colnames(r2) <-c("chr","start","end",suffixes[count],"circRNA_type","gene_id","strand")
    tmp <- merge (tmp, r2, all=TRUE,by=c("chr","start","end"))
    tmp <-data.frame(tmp[,c(1:(ncol(tmp)-7),(ncol(tmp)-3))],circRNA_type=apply(tmp[,c((ncol(tmp)-6),(ncol(tmp)-2))], 1, chooseType),gene_id=apply(tmp[,c((ncol(tmp)-5),(ncol(tmp)-1))], 1, chooseGeneId),strand = apply(tmp[,c((ncol(tmp)-4),ncol(tmp))], 1, chooseStrand),check.names=F)
  }
  count <-count+1
}
tmp[is.na(tmp)]<-0
final1 <- data.frame(circRNA_ID = paste(paste(tmp[,1],tmp[,2],sep=":"),tmp[,3],sep="|"),tmp[,c(4:ncol(tmp))],check.names=F)
final2 <- data.frame(circRNA_ID = paste(paste(tmp[,1],tmp[,2],sep=":"),tmp[,3],sep="|"),tmp[,c(4:ncol(tmp))],check.names=F)
#-----------------------------------------------------------------
# output Result
#-----------------------------------------------------------------
write.table(final2,paste(opt$outfile,"/All_gene_counts.list_tmp",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
write.table(final1,paste(opt$outfile,"/All_gene_counts_detail.xls_tmp",sep=""),sep="\t",quote = FALSE,row.names = FALSE)


