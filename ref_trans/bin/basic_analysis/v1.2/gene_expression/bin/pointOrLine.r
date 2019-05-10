#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript


# load library
library('getopt');
require(ggplot2)
require(RColorBrewer)


#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
	'help' , 'h', 0, "logical",
	'infile' , 'i', 1, "character",
	'outfile' , 'o', 1, "character",
	'x.col' , 'x', 1, "integer",
	'y.col' , 'y', 1, "integer",
	'height' , 'H', 1, "integer",
	'width' , 'W', 1, "integer",
	'x.lab' , 'X', 1, "character",
	'y.lab' , 'Y', 1, "character",
	'title.lab' , 'T', 1, "character",
	'is.point' , 'a', 0, "logical",
	'shape' , 'b', 1, "integer",
	'point.color' , 'c', 1, "integer",
	'is.line' , 'd', 0, "logical",
	'linetype' , 'e', 1, "integer",
	'line.color' , 'f', 1, "integer",
	'line.size' , 'g', 1, "double",
	'point.size' , 'j', 1, "double",
	'lab.size' , 'l', 1, "integer",
	'axis.size' , 's', 1, "integer",
	'no.grid' , 'r', 0, "logical",
	'skip' , 'k', 1, "integer"
	), byrow=TRUE, ncol=4);
opt = getopt(spec);


# define usage function
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
1) Rscript pointOrLine.r --infile in_pointOrLine.data --outfile out_pointOrLine.png \\
	--x.col 1 --y.col 2 --x.lab \"x lab\" --y.lab \"y lab\" --skip 1 \\
	--is.point --is.line --point.color 5 --line.color 4 --shape 16 --linetype 2
2) Rscript pointOrLine.r --infile in_pointOrLine.data --outfile out_pointOrLine.png \\
	--x.col 1 --y.col 2 --x.lab \"x lab\" --y.lab \"y lab\" --skip 1 \\
	--is.point --is.line 
3) Rscript pointOrLine.r --infile in_pointOrLine.data --outfile out_pointOrLine.png \\
	--x.col 1 --y.col 2 --x.lab \"x lab\" --y.lab \"y lab\" --skip 1 \\
	--is.point --is.line --shape 16 --linetype 2
4) Rscript pointOrLine.r --infile in_pointOrLine.data --outfile out_pointOrLine.png \\
	--x.col 1 --y.col 2 --x.lab \"x lab\" --y.lab \"y lab\" --skip 1 \\
	--is.point --is.line --point.color 5 --line.color 4 --shape 16 --linetype 2 \\
	--point.size 2.4 --line.size 0.9 --no.grid

Options: 
--help		-h 	NULL 		get this help
--infile 	-i 	character 	the input file [forced]
--outfile 	-o 	character 	the filename for output graph [forced]
--x.col 	-x 	integer 	the col for x value [forced]
--y.col 	-y 	integer 	the col for y value [forced]
--height 	-H 	integer 	the height of graph [optional, default: 600]
--width 	-W 	integer 	the width of graph [optional, default: 800]
--x.lab 	-X 	character 	the lab for x [forced]
--y.lab 	-Y 	character 	the lab for y [forced]
--title.lab 	-T 	character 	the lab for title [optional, default: NULL]
--is.point 	-a 	NULL 		drawing point
--shape 	-b 	integer 	the shape of point [1-25] [optional, default: 16]
--point.color 	-c 	integer 	the color of point [1-9] [optional, default: 5]
--point.size 	-j 	double 		the size of point [optional, default: 2]
--is.line 	-d 	NULL 		drawing line
--linetype 	-e 	integer 	the type of line [1-5] [optional, default: 1]
--line.color 	-f 	integer 	the color of line [1-9] [optional, default: 4]
--line.size 	-g 	double 		the size of line [optional, default: 1]
--lab.size 	-l 	integer 	the font size of lab [optional, default: 14]
--axis.size 	-s 	integer 	the font size of text for axis [optional, default: 14]
--no.grid	-r 	NULL 		Do not drawing grid
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
if ( is.null(opt$x.col) )	{ print_usage(spec) }
if ( is.null(opt$y.col) )	{ print_usage(spec) }
if ( is.null(opt$x.lab) )	{ print_usage(spec) }
if ( is.null(opt$y.lab) )	{ print_usage(spec) }


#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$skip ) )		{ opt$skip = 0 }
if ( is.null(opt$height ) )		{ opt$height = 600 }
if ( is.null(opt$width ) )		{ opt$width = 800 }
if ( is.null(opt$lab.size ) )		{ opt$lab.size = 14 }
if ( is.null(opt$axis.size ) )		{ opt$axis.size = 14 }
if ( is.null(opt$title.lab) )		{ opt$title.lab = NULL }
if ( is.null(opt$is.point) )		{ opt$is.point = FALSE }
if ( is.null(opt$point.color) )		{ opt$point.color = 5 }
if ( is.null(opt$shape) )		{ opt$shape = 16 }
if ( is.null(opt$is.line) )		{ opt$is.line = FALSE }
if ( is.null(opt$line.color) )		{ opt$line.color = 4 }
if ( is.null(opt$linetype) )		{ opt$linetype = 1 }
if ( is.null(opt$line.size) )		{ opt$line.size = 0.8 }
if ( is.null(opt$point.size) )		{ opt$point.size = 2 }

# check shape linetype and color
opt$point.color <- ((opt$point.color-1) %% 9) + 1 	# here only 9 type of color
opt$line.color <- ((opt$line.color-1) %% 9) + 1 	# here only 9 type of color
opt$shape <- ((opt$shape-1) %% 25) + 1 			# here only 25 type of shape
opt$linetype <- ((opt$linetype-1) %% 5) + 1 		# here only 5 type of linetype

# check combination non-null args
if ( (!opt$is.point) && (!opt$is.line) ) {
	cat("Final Error: Both is.point and is.line are NULL\n")
	print_usage(spec)
}



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
if ( data.dim[2] < max(opt$x.col, opt$y.col) ){
	cat("Final Error: max(x.col, y.col, group.col) > the col of infile\n")
	print_usage(spec)
}
# create df
df <- data.frame(x=data[,opt$x.col], y=data[,opt$y.col])



#-----------------------------------------------------------------
# plot
#-----------------------------------------------------------------
# mian plot
p <- ggplot(df, aes(x=x, y=y))
# good color
colorSet <- brewer.pal(9,"Set1")
cat("opt$point.color",opt$point.color,"\nopt$line.color",opt$line.color,"\n")
# line
if( opt$is.line ){
	p <- p + geom_line(linetype=opt$linetype, colour=colorSet[opt$line.color], size=opt$line.size)
}
# point
if( opt$is.point ){
	p <- p + geom_point(shape=opt$shape, colour=colorSet[opt$point.color], size=opt$point.size)
}

#-----------------------------------------------------------------
# theme
#-----------------------------------------------------------------
# lab
p <- p + xlab(opt$x.lab) + ylab(opt$y.lab) + labs(title=opt$title.lab)
# set lab and axis test size
p <- p + theme(title = element_text(face="bold", size=opt$lab.size), 
	axis.text = element_text(face="bold", size=opt$axis.size))
# remove legend
p <- p + theme(legend.position = "none")
# grid and background
if ( !is.null(opt$no.grid) ) {
	p <- p + theme( panel.background = element_rect(colour="black", size=1, fill="white"),
		panel.grid = element_blank())
}



#-----------------------------------------------------------------
# output plot
#-----------------------------------------------------------------
#pdf(file=paste(opt$outfile,".pdf",sep=""), height=opt$height*2/1000, width=opt$width*2/1000)
#print(p)
#dev.off()
png(filename=opt$outfile, height=opt$height, width=opt$width, res=100, units="px")
print(p)
dev.off()


