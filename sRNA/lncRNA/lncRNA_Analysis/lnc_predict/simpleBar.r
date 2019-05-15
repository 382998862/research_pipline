#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript


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
	'x.col' , 'x', 1, "integer",
	'y.col' , 'y', 1, "integer",
	'height' , 'H', 1, "integer",
	'width' , 'W', 1, "integer",
	'x.lab' , 'X', 1, "character",
	'y.lab' , 'Y', 1, "character",
	'title.lab' , 'T', 1, "character",
	'lab.size' , 'l', 1, "integer",
	'axis.size' , 's', 1, "integer",
	'no.grid' , 'r', 0, "logical",
	'percent_label' , 'A', 0, "logical",
	'count_label' , 'B', 0, "logical",
	'skip' , 'k', 1, "integer"
	), byrow=TRUE, ncol=4);
opt = getopt(spec);


# define usage function
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
1) Rscript simpleBar.r --infile in_simpleBar.data --outfile out_simpleBar.png \\
	--x.col 1 --y.col 2 --x.lab \"x lab\" --y.lab \"y lab\"
2) Rscript simpleBar.r --infile in_simpleBar.data --outfile out_simpleBar.png \\
	--x.col 1 --y.col 2 --x.lab \"x lab\" --y.lab \"y lab\" \\
	--title.lab \"title lab\" --height 3000 --width 4000
3) Rscript simpleBar.r --infile in_simpleBar.data --outfile out_simpleBar.png \\
	--x.col 1 --y.col 2 --x.lab \"x lab\" --y.lab \"y lab\" \\
	--title.lab \"title lab\" --height 3000 --width 4000 --percent_label

Options: 
--help		NULL 		get this help
--infile 	character 	the input file [forced]
--outfile 	character 	the filename for output graph [forced]
--x.col 	integer 	the col for x value [forced]
--y.col 	integer 	the col for y value [forced]
--height 	integer 	the height of graph [optional, default: 3000]
--width 	integer 	the width of graph [optional, default: 4000]
--x.lab 	character 	the lab for x [forced]
--y.lab 	character 	the lab for y [forced]
--title.lab 	character 	the lab for title [optional, default: NULL]
--lab.size 	integer 	the font size of lab [optional, default: 14]
--axis.size 	integer 	the font size of text for axis [optional, default: 14]
--no.grid	NULL 		Do not drawing grid
--percent_label	NULL 		drawing percent label on the bar
--count_label	NULL 		drawing value label on the bar
--skip 		integer 	the number of line for skipping [optional, default: 0]
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
if ( is.null(opt$height ) )		{ opt$height = 3000 }
if ( is.null(opt$width ) )		{ opt$width = 4000 }
if ( is.null(opt$lab.size ) )		{ opt$lab.size = 14 }
if ( is.null(opt$axis.size ) )		{ opt$axis.size = 14 }
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
if ( data.dim[2] < max(opt$x.col, opt$y.col) ){
	cat("Final Error: max(x.col, y.col) > the col of infile\n")
	print_usage(spec)
}
# create df
df <- data.frame(name=as.factor(data[,opt$x.col]), count=data[,opt$y.col])
df[,1]<-factor(df[,1], levels =df[,1])
rownames(df)<-as.vector(df$name)

#-----------------------------------------------------------------
# plot
#-----------------------------------------------------------------
# mian plot
#p <- qplot(x, y, data=df, geom="bar", stat="identity", fill=x)
p <- ggplot(df, aes(x=name,y=count,fill=type))+geom_bar(stat='identity')
#-----------------------------------------------------------------
# theme
#-----------------------------------------------------------------
# lab
p <- ggplot(df, aes(x=name,y=count,fill=name))+geom_bar(stat='identity')
p <- p + geom_text(aes(label = paste(count,"/",round(count/sum(count)*100,1),"%",sep="")), vjust = -0.2)
p <- p + theme_classic()
p<-p+ xlab("")+ylab("number")
p <- p + theme(title = element_text(face="bold", size=14), axis.text = element_text(face="bold", size=10),legend.title= element_text(size=0))
p<- p+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
#p <- p + xlab(opt$x.lab) + ylab(opt$y.lab) + labs(title=opt$title.lab)
# set lab and axis test size
#p <- p + theme(title = element_text(face="bold", size=opt$lab.size), 
#	axis.text = element_text(face="bold", size=opt$axis.size))
# remove legend
#p<- p +  geom_text(aes(label=y),hjust=0.5, vjust=-0.5,size=3 )
#p <- p + theme(legend.position = "none")
#p <- p + theme_bw()
# percent
if( !is.null(opt$percent_label) ){
	p <- p + geom_text(aes(label = paste(round(y/sum(y)*100,1),"%",sep="")), vjust = -0.2)
}
if( (!is.null(opt$count_label)) && (is.null(opt$percent_label)) ){
	p <- p + geom_text(aes(label = y), vjust = -0.2)
}
# grid and background
if ( !is.null(opt$no.grid) ) {
	p <- p + theme( panel.background = element_rect(colour="black", size=1, fill="white"),
		panel.grid = element_blank())
}



#-----------------------------------------------------------------
# output plot
#-----------------------------------------------------------------
pdf(file=paste(opt$outfile,".pdf",sep=""), height=opt$height*2/1000, width=opt$width*2/1000)
print(p)
dev.off()
png(filename=(opt$outfile), height=opt$height, width=opt$width, res=500, units="px")
print(p)
dev.off()







