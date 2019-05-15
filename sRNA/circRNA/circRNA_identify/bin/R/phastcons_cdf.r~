#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript

#file <- read.delim("/home/songmm/下载/scores1.txt")
# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript boxplot.r phastcons_score.txt outdir")
	print("1) phastcons_score: the phasecons score file")
	print("2) outdir: the dir for output")
	print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)

# check args length
# if( length(args) != 2 ) {
# 	print(args)
# 	usage()
# 	stop("the length of args != 2")
# }

# load library
require(ggplot2)

# read count data
print("read phastcons Score ...")
scores <- read.table(args[1],header=TRUE, quote="\"")
head(scores)
print("read phastcons score file is over")

# plot fpkm box for all
# plot

library(ggplot2)
#p<-ggplot(scores, aes(Score, colour =Type)) + stat_ecdf(size=1,geom ="smooth")+xlab("PhastCons Score") + ylab("Cumulative Frequency")
#p<-ggplot(scores, aes(Score)) + stat_ecdf(size=1,geom ="smooth")+xlab("PhastCons Score") + ylab("Cumulative Frequency")
p<-ggplot(scores[,2]) + stat_ecdf(size=1,geom ="smooth")+xlab("PhastCons Score") + ylab("Cumulative Frequency")
p <- p + theme(axis.title.x = element_text(face="bold", size=14),
               axis.text.x  = element_text(face="bold", size=12),
               axis.title.y = element_text(face="bold", size=14),
               axis.text.y  = element_text(face="bold", size=12) )
p <- p + theme(legend.title = element_text(face="bold", size=14),
               legend.text = element_text(size=12)) 
p <- p +theme_classic()
p<-p+scale_x_continuous(limits=c(0,1))
file <- paste(args[2], "/phastcons.cdf.png", sep="")
png(filename=file, height = 3000, width = 3400, res = 500, units = "px")
print(p)
dev.off()







