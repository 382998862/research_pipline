#!/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript

####################################################################
# Copyright 2014, BMK
#
# Author: macx <macx@biomarker.com.cn>
# 
# Function: create pie chart.
# 
######################################################################

#==================================================================
# Get options
#==================================================================

library(getopt)
library(RColorBrewer)

spec <- matrix(c(
	'help', 'h', 0, "logical", "help",
	'input', 'i', 1, "character", "input file, table format:{<key>	<value>}, forced",
	'samName','s',1,  "character","sample Name of pie",
	'output', 'o', 1, "character", "output png file, forced"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

## check options 
if (is.null(opt$input) | is.null(opt$output) | !is.null(opt$help)){
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}
if (is.null(opt$samName)){
	opt$samName="demo";
	
}
#==================================================================
# Main
#==================================================================

## load data 
# The input file format is like this:
#  <label/name>		<frequencies/counts>
data <- read.table(file=opt$input, header=FALSE, comment.char="#", sep="\t")

## preprocess data
slices <- data[,2]
lbls <- data[,1]
lbls <- paste(lbls, slices)
#pct <- round(slices/sum(slices)*100,2)
#lbls <- paste(lbls, pct)  ## add percentage to labels 
#lbls <- paste(lbls, "%", sep="") ## add % to labels 


## plot
mycol<-brewer.pal(8, "Set1")[1:dim(data)[1]]
par(mar=c(5.1,5.1,5.1,6.1))
png(filename=opt$output, width = 3000, height = 3000,res=300)
pie(slices, labels=lbls, col=mycol, border = NA,main=opt$samName,cex.main=2,cex=1.5) #by huangls 2015-2-26
dev.off()

