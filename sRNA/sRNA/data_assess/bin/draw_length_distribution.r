# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript count_rpkm_cor.r in.txt out_prefix")
	print("1) read_count.txt: the read count file for RNA-SEQ")
	print("2) out_file: output file")
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
require(ggplot2)


# read data
print("read  data ...")
df <- read.table(args[1], header = FALSE)
#colnames(df) <- c("length", "num")
head(df)
print("read data is over")
#df=read.table("count_len.stat",head=FALSE)
#p<-qplot(df[,1], df[,2], xlab="Length", ylab="Number", main="Length Graph", size=I(1.5), geom="path")
p<-qplot(df[,1], df[,2], xlab="Length", ylab="Number", main="Length Graph", size=I(1.5),geom="bar" ,stat="identity",fill=factor(df[,1]))
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
png(filename=args[2], height = 3000, width = 6000, res = 500, units = "px")
print(p)
dev.off()



