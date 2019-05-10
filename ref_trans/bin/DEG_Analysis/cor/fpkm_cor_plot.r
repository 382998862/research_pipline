#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript


# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript fpkm_cor_plot.r read_count.txt read_fpkm.txt outdir")
	print("1) read_count.txt: the read count file for RNA-SEQ")
    print("2) read_fpkm.txt: the read fpkm file for RNA-SEQ")
	print("3) outdir: the dir for output")
	print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) < 3 ) {
	print(args)
	usage()
	stop("the length of args < 3")
}


# load library
require(edgeR)
require(ggplot2)


# read count data
print("read count data ...")
count_data <- read.delim(args[1], row.names = 1, header=TRUE,check.names =F)
fpkm_data <- read.delim(args[2], row.names = 1, header=TRUE,check.names =F)
head(count_data)
head(fpkm_data)
print("read count data is over")


# check geneLength key
if( !("geneLength" %in% colnames(count_data)) ) {
	stop("geneLength is not in read count file")
}


# get sample name
sam_name <- NULL
if(length(args)>=4) {
	sam_name <- unlist( strsplit(args[3], ",") )
} else {
	sam_name <- colnames(count_data)[ !(colnames(count_data)%in%c("geneLength")) ]
}
if(length(sam_name)<2) stop("length(sam_name)<2")


# calculate FPKM
fpkm <- fpkm_data[ ,sam_name ]

# init 
decor <- matrix(1,length(sam_name),length(sam_name))
cor_num <- matrix(0,length(sam_name),length(sam_name))
colnames(decor) <- sam_name
rownames(decor) <- sam_name
colnames(cor_num) <- sam_name
rownames(cor_num) <- sam_name



# iter plot fpkm density
for( i in 1:(length(sam_name)-1) ){
	for( j in (i+1):length(sam_name) ){
		# get sample id
		s1 <- sam_name[i]
		s2 <- sam_name[j]
		print(c(s1, s2))

		# filter
		keep_cor <- rowSums(cpm(count_data[,c(s1,s2)])>1) > 1

		# R2 for fpkm
		c1 <- fpkm[,s1][keep_cor]
		c2 <- fpkm[,s2][keep_cor]
		c_na <- (!is.na(c1)) & (!is.na(c2))
		cor_num[i,j] <- sum(c_na)
		cor_num[j,i] <- sum(c_na)
		c <- cor(c1[c_na],c2[c_na])^2
		c <- round(c, 4)
		decor[i,j] <- c
		decor[j,i] <- c

		# plot fpkm cor
		p <- qplot(log10( fpkm[,s1]+1 ), log10( fpkm[,s2]+1 ), xlab=paste("log10(FPKM+1) of ",s1, sep=""), ylab=paste("log10(FPKM+1) of ",s2, sep=""), main="Correlation Plot", size=I(0.5))
		p <- p + geom_abline(slope=1, size=0.5, colour=6, alpha=0.5)
		xc <- 4/5 * min(log10(fpkm[,s1]+1)) + 1/5* max(log10(fpkm[,s1]+1))
		yc <- 1/6 * min(log10(fpkm[,s2]+1)) + 5/6* max(log10(fpkm[,s2]+1))
		p <- p + annotate("text", label = paste("r^2==",c), x = xc, y = yc, size = 5, colour = "blue", parse=TRUE)
		p <- p + theme_classic()
		p<- p+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
		png(filename=paste(args[3],"/cor_plot/",s1,"_vs_",s2,".cor.png", sep=""), height = 750, width = 750, res = 125, units = "px")
		print(p)
		dev.off()
	}
}



# plot heatmap
library(latticeExtra)
library(lattice)

x <- decor
dd.row <- as.dendrogram(hclust(dist(x)))
row.ord <- order.dendrogram(dd.row)

dd.col <- as.dendrogram(hclust(dist(t(x))))
col.ord <- order.dendrogram(dd.col)

p <- levelplot(x[row.ord, col.ord],
	aspect = "fill",
	scales = list(x = list(rot = 90)),
	colorkey = list(space = "left"),
	legend =
	list(right =
		list(fun = dendrogramGrob,
		args =
		list(x = dd.col, ord = col.ord,
		side = "right",
		size = 10)),
		top =
		list(fun = dendrogramGrob,
		args =
		list(x = dd.row,
		side = "top",
		size = 10))), xlab="Sample",ylab="Sample")

# output
file <- paste(args[3], "/sample_cluster.png", sep="")
png(filename=file, height = 800, width = 800, res = 125, units = "px")
print(p)
dev.off()

file <- paste(args[3], "/sample_cluster.pdf", sep="")
pdf(file)
print(p)
dev.off()


# output cor and cor number
df <- data.frame(decor=decor, cor_num=cor_num)
file <- paste(args[3], "/free_com.cor", sep="")
write.table(df, file=file, sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)



