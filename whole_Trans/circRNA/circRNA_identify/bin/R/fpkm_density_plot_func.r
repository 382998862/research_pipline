#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript


# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript fpkm_density_plot_func.r read_count.txt outdir")
	print("1) read_count.txt: the read count file for RNA-SEQ")
	print("2) outdir: the dir for output")
	print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)
# load library
require(edgeR)
require(ggplot2)


# read count data
print("read count data ...")
count_data <- read.delim(args[1], row.names = 1, header=TRUE,check.names = FALSE)
head(count_data)
print("read count data is over")
# plot
m <- ggplot(count_data, aes(x=count_data[,1]))
p <- m + geom_density(fill="#CC79A7", size=1, colour="#CC79A7") + xlab("Conservation Score")
p <- p + theme(axis.title.x = element_text(face="bold", size=14),
	axis.text.x  = element_text(face="bold", size=12),
	axis.title.y = element_text(face="bold", size=14),
	axis.text.y  = element_text(face="bold", size=12) )
p <- p + theme(legend.title = element_text(face="bold", size=14),
	legend.text = element_text(size=12) )
# output
file <- paste(args[2],"/phastcons.cdf.png", sep="")
png(filename=file, height = 3000, width = 3400, res = 500, units = "px")
print(p)
dev.off()


## plot fpkm box for all
## plot
#m <- ggplot(log10fpkm, aes(factor(sample), log10fpkm))
#p <- m + geom_boxplot(aes(fill=Sample)) + xlab("Sample") + ylab("log10(FPKM)")
#p <- p + theme(axis.title.x = element_text(face="bold", size=14),
#	axis.text.x  = element_text(face="bold", size=12),
#	axis.title.y = element_text(face="bold", size=14),
#	axis.text.y  = element_text(face="bold", size=12) )
#p <- p + theme(legend.title = element_text(face="bold", size=14),
#	legend.text = element_text(size=12) )
## output
#file <- paste(args[2], "/all", ".fpkm_box.png", sep="")
#png(filename=file, height = 3000, width = 3400, res = 500, units = "px")
#print(p)
#dev.off()







