#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript


# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript fpkm_density_plot_func.r read_fpkm.txt outdir")
	print("1) fpkm.txt: fpkm file for RNA-SEQ")
	print("2) count.txt: count file for RNA-SEQ")
	print("3) outdir: the dir for output")
	print("4) expression type:FPKM/TPM/SRPBM")
	print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) != 4 ) {
	print(args)
	usage()
	stop("the length of args != 4")
}


# load library
require(ggplot2)

# read count data
print("read fpkm data ...")
fpkm_data <- read.delim(args[1], row.names = 1, header=TRUE,check.names =F)
print("read fpkm data is over")
sam_name <- colnames(fpkm_data)
fpkm <- as.matrix(fpkm_data[ ,sam_name ])
print("calculate FPKM is over")

print("read count data ...")
count_data <- read.delim(args[2], row.names = 1, header=TRUE,check.names =F)
count_data<-as.matrix(count_data)
print("read count data is over")

fpkm<-fpkm[rownames= rownames(count_data),]
log10_fpkm <- fpkm

# init 
all <- NULL
all_sam <- NULL
lab <- paste("log10(",args[4],")",sep = "")
# iter plot fpkm density
for( i in 1:length(sam_name) ){
	c <- count_data[,sam_name[i]]
	keep <- c > 0
	r<-fpkm[keep,sam_name[i]]+0.00001
	log10fpkm <- data.frame(log10fpkm=log10(r))
	
	# update all
	all <- c(all, log10(r))
	all_sam <- c(all_sam, rep(sam_name[i], length(r)))

	# plot
	m <- ggplot(log10fpkm, aes(x=log10fpkm))
	p <- m + geom_density(fill="#CC79A7", size=1, colour="#CC79A7") + xlab(lab)
	p <- p + theme_classic()
	p<- p+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))

	p <- p + theme(axis.title.x = element_text(face="bold", size=14),
		axis.text.x  = element_text(face="bold", size=12),
		axis.title.y = element_text(face="bold", size=14),
		axis.text.y  = element_text(face="bold", size=12) )
	p <- p + theme(legend.title = element_text(face="bold", size=14),
		legend.text = element_text(size=12) )
	if (length(sam_name)>20) {p<- p + guides(fill = guide_legend(nrow = 20))}
	# output plot into file
#	file <- paste(args[3], "/", sam_name[i], ".fpkm_density.png", sep="")
#	png(filename=file, height = 3000, width = 3000, res = 500, units = "px")
#	print(p)
#	dev.off()
}


# plot fpkm density for all
# create data.frame
log10fpkm <- data.frame(log10fpkm=all, sample=all_sam)
Sample <- factor(all_sam)
#write.table(log10fpkm,"all_temp.xls")
# plot
m <- ggplot(log10fpkm, aes(x=log10fpkm))
p <- m + geom_density(aes(fill=Sample, colour=Sample),alpha = 0.2) + xlab(lab)
p <- p + theme(axis.title.x = element_text(face="bold", size=14),
	axis.text.x  = element_text(face="bold", size=12),
	axis.title.y = element_text(face="bold", size=14),
	axis.text.y  = element_text(face="bold", size=12) )
p <- p + theme(legend.title = element_text(face="bold", size=14),
	legend.text = element_text(size=12) )
if (length(sam_name)>20) {p<- p + guides(fill = guide_legend(nrow = 20))}
p <- p + theme_classic()
p<- p+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
# output
file <- paste(args[3], "/all", ".fpkm_density.png", sep="")
png(filename=file, height = 3000, width = 3400, res = 500, units = "px")
print(p)
dev.off()


## plot fpkm box for all
## plot
m <- ggplot(log10fpkm, aes(factor(sample), log10fpkm))

if(length(sam_name)>20){
p <- m + geom_boxplot(aes(fill=Sample),outlier.shape=NA) + xlab("Sample") + ylab(lab)
p <- p + theme_classic()
p<- p+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
p <- p + theme(axis.title.x = element_text(face="bold", size=14),
	axis.text.x  = element_text(face="bold", size=12),
	axis.title.y = element_text(face="bold", size=14),
	axis.text.y  = element_text(face="bold", size=12) )

p <- p + theme(legend.title = element_text(face="bold", size=14),
	legend.text = element_text(size=12) )
if (length(sam_name)>20) {p<- p + guides(fill = guide_legend(nrow = 20))}
# output
file <- paste(args[3], "/all", ".fpkm_box.png", sep="")
png(filename=file, height = 5000, width = 9000, res = 500, units = "px")
print(p)
dev.off()
}else{
p <- m + geom_boxplot(aes(fill=Sample)) + xlab("Sample") + ylab(lab)
p <- p + theme_classic()
p<- p+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
p <- p + theme(axis.title.x = element_text(face="bold", size=14),
        axis.text.x  = element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=14),
        axis.text.y  = element_text(face="bold", size=12) )

p <- p + theme(legend.title = element_text(face="bold", size=14),
        legend.text = element_text(size=12) )
if (length(sam_name)>20) {p<- p + guides(fill = guide_legend(nrow = 20))}
# output
file <- paste(args[3], "/all", ".fpkm_box.png", sep="")
png(filename=file, height = 3000, width = 2400, res = 500, units = "px")
print(p)
dev.off()
}
