#!/share/nas2/genome/biosoft/R/current/bin/Rscript


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
	'main.col' , 'm', 1, "integer",
	'sub.col' , 's', 1, "integer",
	'value.col' , 'v', 1, "integer",
	'main.lab' , 'M', 1, "character",
	'sub.lab' , 'S', 1, "character",
	'value.lab' , 'V', 1, "character",
	'title.lab' , 'T', 1, "character",
	'legend.xpos' , 'x', 1, "double",
	'legend.ypos' , 'y', 1, "double",
	'skip' , 'k', 1, "integer"
	), byrow=TRUE, ncol=4);
opt = getopt(spec);


# define usage function
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("1) Rscript bar_split.r --infile t.data --outfile out.png --main.col 1 --sub.col 2 \\\n")
	cat("\t--value.col 3 --main.lab mainlab --sub.lab sublab --value.lab valuelab \\\n")
	cat("\t--title.lab titlelab --legend.xpos 0.8 --legend.ypos 0.85\n")
	cat("2) Rscript bar_split.r --infile t.data --outfile out.png --main.col 1 --sub.col 2 \\\n")
	cat("\t--value.col 3 --main.lab mainlab --sub.lab sublab --value.lab valuelab \\\n")
	cat("\t--title.lab titlelab\n")
	cat("3) Rscript bar_split.r --infile t.data --outfile out.png --main.col 1 --sub.col 2 \\\n")
	cat("\t--value.col 3 --main.lab mainlab --sub.lab sublab --value.lab valuelab\n")
	q(status=1);
}


# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }


# check non-null args
if ( is.null(opt$infile) )	{ print_usage(spec) }
if ( is.null(opt$outfile) )	{ print_usage(spec) }
if ( is.null(opt$main.col) )	{ print_usage(spec) }
if ( is.null(opt$sub.col) )	{ print_usage(spec) }
if ( is.null(opt$value.col) )	{ print_usage(spec) }
if ( is.null(opt$main.lab) )	{ print_usage(spec) }
if ( is.null(opt$sub.lab) )	{ print_usage(spec) }
if ( is.null(opt$value.lab) )	{ print_usage(spec) }


#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$skip ) )		{ opt$skip = 0 }
if ( is.null(opt$legend.xpos ) )	{ opt$legend.xpos = 0.7 }
if ( is.null(opt$legend.ypos ) )	{ opt$legend.ypos = 0.85 }
if ( is.null(opt$title.lab) )		{ opt$title.lab = NULL }



#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------
# reading data
data <- read.table(opt$infile, skip=opt$skip)
# check dim
data.dim <- dim(data)
if ( is.null(data.dim) ){
	cat("Final Error: the format of infile is error, dim(data) is NULL\n")
	print_usage(spec)
}
# check col size
if ( data.dim[2] < max(opt$main.col, opt$sub.col, opt$value.col) ){
	cat("Final Error: max(main.col, sub.col, value.col) > the col of infile\n")
	print_usage(spec)
}


#-----------------------------------------------------------------
# plot
#-----------------------------------------------------------------
# plot
p <- ggplot(data, aes(x=as.factor(data[,opt$main.col]), 
	y=as.numeric(data[,opt$value.col]), fill=as.factor(data[,opt$sub.col])))
p <- p + geom_bar(position="dodge",stat="identity") + labs(colour="Type")
# lab
p <- p + labs(x=opt$main.lab) + labs(y=opt$value.lab) + labs(title=opt$title.lab)
p <- p + labs(fill=opt$sub.lab)
# size and style
size <- 18
p <- p + theme(axis.title.x = element_text(face="bold", size=size),
	axis.text.x  = element_text(face="bold", size=size),
	plot.title = element_text(face="bold", size=size),
	axis.title.y = element_text(face="bold", size=size),
	axis.text.y  = element_text(face="bold", size=size) )
p <- p + theme(legend.title = element_text(face="bold", size=size),
	legend.text = element_text(size=size) )
p <- p + theme(legend.position=c(opt$legend.xpos, opt$legend.ypos))


# output plot
png(filename=opt$outfile, height = 3000, width = 6000, res = 500, units = "px")
print(p)
dev.off()






