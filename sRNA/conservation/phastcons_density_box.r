#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript


# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript phastcons_density_box.r phastcons_score.txt outdir")
	print("1) phastcons_score: the phasecons score file(Contained at least 2 column:mean0 type)")
	print("2) outdir: the dir for output")
	print("3) key: output key")
	print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) != 2 && length(args)!=3) {
	print(args)
	usage()
	stop("the length of args != 2/3")
}

# load library
require(ggplot2)

# read count data
print("read phastcons Score ...")
scores <- read.csv(args[1],header=TRUE, sep="\t")
head(scores)
print("read phastcons score file is over")

# plot fpkm box for all
# plot

library(ggplot2)
p<-ggplot(data=scores,mapping=aes(x=mean0,colour=type)) + geom_density()+xlab("PhastCons Score") + ylab("Density")
p <- p + theme(axis.title.x = element_text(face="bold", size=14),
               axis.text.x  = element_text(face="bold", size=12),
               axis.title.y = element_text(face="bold", size=14),
               axis.text.y  = element_text(face="bold", size=12) )
p <- p + theme(legend.title = element_text(face="bold", size=14),
               legend.text = element_text(size=12)) 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p<-p+scale_x_continuous(limits=c(0,1))

key="RNA"
if(length(args)==3){key=args[3]}

file <- paste(args[2], "/",key,".phastcons.density.png", sep="")
png(filename=file, height = 3000, width = 3400, res = 500, units = "px")
print(p)
dev.off()

file <- paste(args[2], "/",key,".phastcons.density.pdf", sep="")
pdf(file)
print(p)
dev.off()



p<-ggplot(data=scores,mapping=aes(x=type,y=mean0)) + geom_boxplot(aes(fill=type))+xlab("RNA type") + ylab("PhastCons Score")
p<-p+theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title.x = element_text(face="bold", size=14),
               axis.text.x  = element_text(face="bold", size=12),
               axis.title.y = element_text(face="bold", size=14),
               axis.text.y  = element_text(face="bold", size=12) )
p <- p + theme(legend.title = element_text(face="bold", size=14),
               legend.text = element_text(size=12))

file <- paste(args[2], "/",key,".phastcons.boxplot.png", sep="")
png(filename=file, height = 3000, width = 3400, res = 500, units = "px")
print(p)
dev.off()

file <- paste(args[2], "/",key,".phastcons.boxplot.pdf", sep="")
pdf(file)
print(p)
dev.off()



p<-ggplot(data=scores,mapping=aes(x=mean0,colour=type)) + stat_ecdf(size=1,geom ="smooth")+xlab("PhastCons Score") + ylab("Cumulative Frequency")
p <- p + theme(axis.title.x = element_text(face="bold", size=14),
               axis.text.x  = element_text(face="bold", size=12),
               axis.title.y = element_text(face="bold", size=14),
               axis.text.y  = element_text(face="bold", size=12) )
p <- p + theme(legend.title = element_text(face="bold", size=14),
               legend.text = element_text(size=12))
p <- p +theme_classic()
p<-p+scale_x_continuous(limits=c(0,1))

if(length(args)==3){key=args[3]}

file <- paste(args[2], "/",key,".phastcons.cdf.png", sep="")
png(filename=file, height = 3000, width = 3400, res = 500, units = "px")
print(p)
dev.off()

file <- paste(args[2], "/",key,".phastcons.cdf.pdf", sep="")
pdf(file)
print(p)
dev.off()

