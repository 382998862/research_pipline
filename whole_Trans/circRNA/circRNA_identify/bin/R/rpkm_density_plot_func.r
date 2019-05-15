#!/share/nas2/genome/biosoft/R/current/bin/Rscript

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

# check args length
if( length(args) != 2 ) {
	print(args)
	usage()
	stop("the length of args != 2")
}


# load library
require(ggplot2)


# read count data
print("read count data ...")
count_data <- read.delim(args[1], row.names = 1, header=TRUE)
head(count_data)
print("read count data is over")




sample <- colnames(count_data)
# init 
all <- NULL
all_sam <- NULL

# iter plot fpkm density
for( i in 1:length(sample) ){
  # fpkm
  c <- count_data[,sample[i]]
  keep <- c > 0
  r <- c[keep]
  log10fpkm <- data.frame(log10fpkm=log10(r))
  
  # update all
  all <- c(all, log10(r))
  all_sam <- c(all_sam, rep(sample[i], length(r)))
}


log10fpkm <- data.frame(log10fpkm=all, sample=all_sam)
Sample <- factor(all_sam)
# output
m <- ggplot(log10fpkm, aes(factor(sample), log10fpkm))
p <- m + geom_boxplot(aes(fill=Sample)) + xlab("Sample") + ylab("log10(FPKM)")
p <- p + theme(axis.title.x = element_text(face="bold", size=14),
               axis.text.x  = element_text(face="bold", size=6),
               axis.title.y = element_text(face="bold", size=14),
               axis.text.y  = element_text(face="bold", size=12) )
p <- p + theme(legend.title = element_text(face="bold", size=14),
               legend.text = element_text(size=12) )
p <- p + theme( panel.background = element_rect(colour="black", size=1, fill="white"),panel.grid = element_blank())
# output
file <- paste(args[2], ".fpkm_box.png", sep="")
png(filename=file, height = 3000, width = 4500, res = 500, units = "px")
print(p)
dev.off()












