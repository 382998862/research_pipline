# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript count_rpkm_cor.r in.txt out_prefix")
	print("1) read_count.txt: the read count file for RNA-SEQ")
	print("2) output dir")
	print("3) outfile prefix")
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
require(ggplot2)


# read data
print("read  data ...")
df <- read.table(args[1], header = FALSE)
#colnames(df) <- c("length", "num")
head(df)
print("read data is over")
#df=read.table("count_len.stat",head=FALSE)
#p<-qplot(df[,1], df[,2], xlab="Length", ylab="Number", main="Length Graph", size=I(1.5), geom="path")
sample<-unlist(strsplit(args[3],split="[._]"))
title<-paste("Length Graph","(",sample[1],")",sep="")
p<-qplot(df[,1], df[,2], xlab="Length", ylab="Number", main=title, size=I(1.5),geom="bar" ,stat="identity",fill=factor(df[,1]))
p<-p+theme(legend.position="none",panel.background = element_rect(fill = "transparent", color = "white"),axis.line = element_line(colour = "gray",size = 1))
#p<-p+theme_bw()
#p<-p+theme(panel.grid.major=element_line(colour=NA))
p <- p + scale_x_continuous(breaks=min(df[,1]):max(df[,1]))
## plot
##p <- ggplot(df, aes(x=Site, y=Percent, group=Type))
#p <- ggplot(df, aes(x=length, y=num))
##p <- p + geom_line(aes(colour = Type)) + labs(colour="Type")
#size <- 18
#p <- p + theme(axis.title.x = element_text(face="bold", size=size),
#	axis.text.x  = element_text(face="bold", size=size),
#	axis.title.y = element_text(face="bold", size=size),
#	axis.text.y  = element_text(face="bold", size=size) )
#p <- p + theme(legend.title = element_text(face="bold", size=size),
#	legend.text = element_text(size=size) )
#p <- p + theme(legend.position=c(.7, .85))


# output plot
png(filename=paste(args[2],"/",args[3],".length.png",sep=""), height = 3000, width = 6000, res = 500, units = "px")
print(p)
dev.off()
pdf(file=paste(args[2],"/",args[3],".length.pdf",sep=""), height = 7, width = 14)
print(p)
dev.off()
