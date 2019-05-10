#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript
#####################################################################
# Copyright 2015, BMK
#
# Programe name: plot_SNP_stat.R
#
# Author: tengh <tengh@biomarker.com.cn>
#
# Function: draw statics of SNP map
#
# Create Date: 20150803
#
# Modified info:
#####################################################################

# load library
library('getopt');
require(ggplot2)

#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
	'help' , 'h', 0, "logical",
	'infile' , 'i', 1, "character",
	'outfile' , 'o', 1, "character",
	'group.col' , 'g', 1, "integer",
	'x.col' , 'x', 1, "integer",
	'y.col' , 'y', 1, "integer",
	'height' , 'H', 2, "integer",
	'width' , 'W', 2, "integer",
	'group.lab' , 'G', 2, "character",
	'x.lab' , 'X', 2, "character",
	'y.lab' , 'Y', 2, "character",
	'title.lab' , 'T', 2, "character",
	'legend.xpos' , 'a', 2, "double",
	'legend.ypos' , 'b', 2, "double",
	'legend.col' , 'c', 2, "integer",
	'lab.size' , 'l', 2, "integer",
	'axis.size' , 's', 2, "integer",
	'legend.size' , 'd', 2, "integer",
	'no.grid' , 'r', 2, "logical",
	'skip' , 'k', 1, "integer"
	), byrow=TRUE, ncol=4);
opt = getopt(spec);

# define usage function
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
1) Rscript dodgedBar.r --infile in_dodgedBar.data --outfile out_dodgedBar.png \\
	--x.col 1 --group.col 2 --y.col 3 --x.lab \"x lab\" --group.lab \"group lab\" \\
	--y.lab \"y lab\" --title.lab \"title lab\" --legend.col 2
2) Rscript dodgedBar.r --infile in_dodgedBar.data --outfile out_dodgedBar.png \\
	--x.col 1 --group.col 2 --y.col 3 --x.lab \"x lab\" --group.lab \"group lab\" \\
	--y.lab \"y lab\" --title.lab \"title lab\" --height 600 --width 800
3) Rscript dodgedBar.r --infile in_dodgedBar.data --outfile out_dodgedBar.png \\
	--x.col 1 --group.col 2 --y.col 3 --x.lab \"x lab\" --group.lab \"group lab\" \\
	--y.lab \"y lab\" --title.lab \"title lab\" --legend.xpos 0.8 --legend.ypos 0.85

Options: 
--help		-h 	NULL 		get this help
--infile 	-i 	character 	the input file [forced]
--outfile 	-o 	character 	the filename for output graph [forced]
--group.col 	-g 	integer 	the col for group factor [forced]
--x.col 	-x 	integer 	the col for x value [forced]
--y.col 	-y 	integer 	the col for y value [forced]
--height 	-H 	integer 	the height of graph [optional, default: 600]
--width 	-W 	integer 	the width of graph [optional, default: 800]
--group.lab 	-G 	character 	the lab for group factor [forced]
--x.lab 	-X 	character 	the lab for x [forced]
--y.lab 	-Y 	character 	the lab for y [forced]
--title.lab 	-T 	character 	the lab for title [optional, default: NULL]
--legend.xpos 	-a 	double 		the x relative position for legend, (0.0,1.0) [optional, default: NULL]
--legend.ypos 	-b 	double 		the y relative position for legend, (0.0,1.0) [optional, default: NULL]
--legend.col 	-c 	integer 	the col number for legend disply [optional, default: NULL]
--lab.size 	-l 	integer 	the font size of lab [optional, default: 14]
--axis.size 	-s 	integer 	the font size of text for axis [optional, default: 14]
--legend.size 	-d 	integer 	the font size of text for legend [optional, default: 12]
--no.grid	-r 	NULL 		Do not drawing grid [optional, default: FALSE]
--skip 		-k 	integer 	the number of line for skipping [optional, default: 0]
\n")
	q(status=1);
}


# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }


# check non-null args
if ( is.null(opt$infile) )	{ print_usage(spec) }
if ( is.null(opt$outfile) )	{ print_usage(spec) }
if ( is.null(opt$group.col) )	{ print_usage(spec) }
if ( is.null(opt$x.col) )	{ print_usage(spec) }
if ( is.null(opt$y.col) )	{ print_usage(spec) }



#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$group.lab) )  { opt$group.lab="Group" }
if ( is.null(opt$x.lab) )	{ opt$x.lab=""}
if ( is.null(opt$y.lab) )	{opt$y.lab="counts" }
if ( is.null(opt$skip ) )		{ opt$skip = 0 }
if ( is.null(opt$height ) )		{ opt$height = 600 }
if ( is.null(opt$width ) )		{ opt$width = 800 }
if ( is.null(opt$legend.xpos ) )	{ opt$legend.xpos = NULL }
if ( is.null(opt$legend.ypos ) )	{ opt$legend.ypos = NULL }
if ( is.null(opt$legend.col ) )		{ opt$legend.col = NULL }
if ( is.null(opt$lab.size ) )		{ opt$lab.size = 14 }
if ( is.null(opt$axis.size ) )		{ opt$axis.size = 14 }
if ( is.null(opt$legend.size ) )	{ opt$legend.size = 12 }
if ( is.null(opt$title.lab) )		{ opt$title.lab = NULL }
if ( is.null(opt$no.grid) )  	{ opt$no.grid = FALSE }


#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------
# reading data
data <- read.table(as.vector(opt$infile), skip=opt$skip,comment.char = "#",header=F)
header<- as.matrix(read.table(as.vector(opt$infile), skip=opt$skip,comment.char = "",header=F,nrows = 1))
if(grep("^#",x=header[[1]])>0&&ncol(data)==length(header)){
  colnames(data)=header
}
# check dim
data.dim <- dim(data)
if ( is.null(data.dim) ){
	cat("Final Error: the format of infile is error, dim(data) is NULL\n")
	print_usage(spec)
}
# check col size
if ( data.dim[2] < max(opt$x.col, opt$y.col, opt$group.col) ){
	cat("Final Error: max(x.col, y.col, group.col) > the col of infile\n")
	print_usage(spec)
}
# create df
#data<-data[data[,1]!="All",]
df <- data.frame(x=as.numeric(as.factor(data[,opt$x.col])),x.lab=as.factor(data[,opt$x.col]), group=as.factor(data[,opt$group.col]), y=data[,opt$y.col])

#-----------------------------------------------------------------
# plot
#-----------------------------------------------------------------
# mian plot

p <- ggplot(df, aes(x=x, y=y, group=group,colour=group))
p<-p+geom_line( size = 1)
p<-p+geom_point(size=3)
#-----------------------------------------------------------------
# theme
#-----------------------------------------------------------------
# lab

p <- p + labs(fill=opt$group.lab) + xlab(opt$x.lab) + ylab(opt$y.lab) + labs(title=opt$title.lab)
# set lab and axis test size
p <- p + theme(title = element_text(face="bold", size=opt$lab.size), 
	axis.text = element_text(face="bold", size=opt$axis.size))
p <- p + theme(legend.title = element_text(face="bold", size=opt$legend.size),
	legend.text = element_text(size=opt$legend.size) )
# legend position
if( (!is.null(opt$legend.xpos)) && (!is.null(opt$legend.ypos)) ){
	p <- p + theme(legend.position=c(opt$legend.xpos, opt$legend.ypos))
}
x_breaks<-as.factor(levels(as.factor(data[,opt$x.col])))
p <- p + scale_x_continuous(breaks=as.numeric(x_breaks), labels=x_breaks) 
# legend col
#if( is.null(opt$legend.col) ){
#	levels_num <- length( levels(df$group) )
#	if( levels_num%%12==0 )
#		opt$legend.col <- as.integer(levels_num / 12)
#	else
#		opt$legend.col <- as.integer(levels_num / 12) + 1
#}
#p <- p + guides(fill = guide_legend(ncol = opt$legend.col))

if (length(unique(df$group))>20){p<-p + guides(fill = guide_legend(nrow=20))}

# grid and background
if ( opt$no.grid) {
	p <- p + theme( panel.background = element_rect(colour="black", size=1, fill="white"),
		panel.grid = element_blank())
}



#-----------------------------------------------------------------
# output plot
#-----------------------------------------------------------------
#pdf(file=paste(opt$outfile,".pdf",sep=""), height=opt$height*2/1000, width=opt$width*2/1000)
#print(p)
#dev.off()
png(filename =opt$outfile, height=opt$height, width=opt$width, res=100, units="px")
print(p)
dev.off()







