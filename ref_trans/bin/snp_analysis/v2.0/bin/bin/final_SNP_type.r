#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript


# usage function
usage <- function() {
  print("-------------------------------------------------------------------------------")
  print("Usage: Rscript bin/cog_anno_plot.r in.stat out.png")
  print("1) in.stat: the stat file for SNP type")
  print("2) outdir: the dirname for output")
  print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) != 2 ) {
  print(args)
  usage()
  stop("the length of args != 2")
}



# load library
library(ggplot2)
library(grid)


# get args
args <-commandArgs(TRUE)


# reading data
data <- read.delim(args[1], header=TRUE, sep="\t")
colnames(data) <- read.delim(args[1],header=F,check.names =F,stringsAsFactors=F,nrows=1)
head(data)

Sample<-colnames(data[,-1])

# plot
for (i in 1:length(Sample))
{
#  df <- data.frame(Number=data[,i+1], type=factor(data[,1],levels=factor(data[,1][nrow(data):1])))
   df <- data.frame(Number=data[,i+1], type=factor(data[,1],levels=factor(data[,1][1:nrow(data)])))
    p <- ggplot(data=df, aes(x=type, y=Number)) + geom_bar(aes(fill=type), stat="identity")
    p <- p+ theme(axis.text.x = element_text(colour ="black",angle = 45, hjust = 1))
  #p<-p+coord_flip()+guides(fill=FALSE) + theme_bw()
  p<-p+labs(x="")+guides(fill=FALSE)
  p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  # output plot	
  png(filename=paste(args[2],paste(Sample[i],"snp.type.png",sep="."),sep='/'), height = 2500, width = 3000, res = 500, units = "px")
  print(p)
  dev.off()
}









