#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript


# usage function
usage <- function() {
  print("-------------------------------------------------------------------------------")
  print("Usage: Rscript bin/cog_anno_plot.r in.stat out.png")
  print("1) in.stat: the stat file for SNP/InDel Annoation")
  print("2) outdir: the dirname for output")
  print("3) keyword: the keyword for output")
  print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) != 3 ) {
  print(args)
  usage()
  stop("the length of args != 3")
}



# load library
library(ggplot2)
library(grid)


# get args
args <-commandArgs(TRUE)


# reading data
data <- read.delim(args[1], header=TRUE, sep="\t")
colnames(data) <- read.delim(args[1],header=F,check.names =F,stringsAsFactors=F,nrows=1)
data<-data[,-1]
head(data)

Sample<-colnames(data[,-1])

# plot
for (i in 1:length(Sample))
{

  df <- data.frame(Number=data[,i+1], type=factor(data$Type,levels=factor(data$Type[nrow(data):1])))
  p <- ggplot(data=df, aes(x=type, y=Number)) + geom_bar(aes(fill=type), stat="identity")
  p<-p+coord_flip()+guides(fill=FALSE) #+ theme_bw()
  p<-p+geom_text(label=data[,i+1], hjust=-0.5, size=3.5)+ylim(c(0,max(data[,i+1]*1.1)))
  p<-p+labs(x="Type", title="Number of effects by type")
  p <- p+ theme(axis.text.y = element_text(colour ="black"))
  p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  # output plot	
  png(filename=paste(args[2],paste(Sample[i],args[3],"anno.stat.png",sep="."),sep='/'), height = 3000, width = 5000, res = 500, units = "px")
  print(p)
  dev.off()
}









