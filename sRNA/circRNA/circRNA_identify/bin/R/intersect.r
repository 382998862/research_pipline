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
  'ciri' , 'c', 1, "character",
  'find_circ','f',0, "character",
  'circexplorer','e',0, "character",
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
      --ciri 	character 	the ciri input file [forced]
      --find_circ 	character 	the find_circ input file
      --circexplorer 	character 	the circexplorer input file
      --outfile 	character 	the filename for output graph [forced]
      \n")
  q(status=1);
}



# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }


# check non-null args
if ( is.null(opt$ciri) )	{ print_usage(spec) }
if ( is.null(opt$outfile) )	{ print_usage(spec) }


#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------
# reading data
rawdata <- list ()
rawdata[[1]] <- read.table(file = opt$ciri, sep = "\t", header=FALSE, comment.char = "!",skip = 1)
tmp <- rawdata[[1]]
file.remove(opt$ciri)
if (!is.null(opt$find_circ))
{
  rawdata[[2]] <- read.table(file = opt$find_circ, sep = "\t", header=FALSE, comment.char = "!",skip = 1)
  file.remove(opt$find_circ)
  rawdata[[2]][,2]<- as.numeric(rawdata[[2]][,2])+1
  head(rawdata[[2]])
  tmp <- merge (tmp, rawdata[[2]], by=c("V1", "V2", "V3"))
  head(tmp)
  tmp<-data.frame(tmp[,1:3],junction=apply(tmp[,c(4,ncol(tmp))],1,mean),tmp[,(ncol(tmp)-3):(ncol(tmp)-1)])
}
if(!is.null(opt$circexplorer))
{
  rawdata[[3]] <- read.table(file = opt$circexplorer, sep = "\t", header=FALSE, comment.char = "!",skip = 1)
  head(rawdata[[3]])
  file.remove(opt$circexplorer)
  rawdata[[3]][,2]<- as.numeric(rawdata[[3]][,2])+1
  tmp <- merge (tmp, rawdata[[3]], by=c("V1", "V2", "V3"))
  head(tmp)
  tmp<-data.frame(tmp[,1:3],junction=apply(tmp[,c(4,ncol(tmp))],1,mean),tmp[,(ncol(tmp)-3):(ncol(tmp)-1)])
}
tmp[,4]<- as.integer(tmp[,4])
colnames(tmp) <-c("chr","start","end","Mean_junction","circRNA_type","gene_id","strand")
#-----------------------------------------------------------------
# output Result
#-----------------------------------------------------------------
write.table(tmp,opt$outfile,sep="\t",quote = FALSE,row.names = FALSE)




