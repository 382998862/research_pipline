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
	'sample','sa', 1, "character",
	'x.col' , 'x', 1, "integer",
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
	--x.col 1  --x.lab \"x lab\" --y.lab \"y lab\"
2) Rscript simpleBar.r --infile in_simpleBar.data --outfile out_simpleBar.png \\
	--x.col 1 --x.lab \"x lab\" --y.lab \"y lab\" \\
	--title.lab \"title lab\" --height 3000 --width 4000
3) Rscript simpleBar.r --infile in_simpleBar.data --outfile out_simpleBar.png \\
	--x.col 1 --x.lab \"x lab\" --y.lab \"y lab\" \\
	--title.lab \"title lab\" --height 3000 --width 4000 --percent_label

Options:
--help		NULL 		get this help
--infile 	character 	the input file [forced]
--outfile 	character 	the filename for output graph [forced]
--sample        character       sample name split by ","
--x.col 	integer 	the col for x value [forced]
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
if ( is.null(opt$x.lab) )	{ print_usage(spec) }
if ( is.null(opt$y.lab) )	{ print_usage(spec) }
if ( is.null(opt$sample) )	{ print_usage(spec) }

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

data <- read.table(opt$infile, header=TRUE,sep ="\t",check.names=F)
head(data)
# check dim
data.dim <- dim(data)
if ( is.null(data.dim) ){
	cat("Final Error: the format of infile is error, dim(data) is NULL\n")
	print_usage(spec)
}
sample<-unlist(strsplit(opt$sample,','))
for(s in sample)
{
    sample_res <- data.frame(chr = apply(data,1,function(x){return (unlist(strsplit(x,':'))[1])} ),junction = data[,s])
    head(sample_res)
    df <- aggregate(sample_res[,2], by=list(y=sample_res[,1]), FUN=sum)
    #df <- df[df$x>1000,]
    df <- data.frame(x=as.factor(df$y), y=df$x)
    df <- df[grep("^(chr)?[0-9]+$|^(chr)?X$|^(chr)?Y$", df$x,ignore.case=TRUE),]
#    df$x <- factor(paste(df$x),levels=paste(df$x))
    df$x<-factor(df$x,levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))
    head(df)
#data$group<-factor(data$group,levels=c("Dormant","Reactivating","Active"))
    #-----------------------------------------------------------------
    # plot
    #-----------------------------------------------------------------
    # mian plot
    #p <- qplot(x, y, data=df, geom="bar", stat="identity", fill=x)
    p<- ggplot(df, aes(x=x,y=y,fill=x))+geom_bar(stat='identity')
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
    p <- p + theme(axis.text.x = element_text(angle = 60, hjust = 1))
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
    pdf(file=paste(opt$outfile,s,".chromosome_distrbution.pdf",sep=""), height=opt$height*2/1000, width=opt$width*2/1000)
    print(p)
    dev.off()
    png(filename=paste(opt$outfile,s,".chromosome_distrbution.png",sep=""), height=opt$height, width=opt$width, res=500, units="px")
    print(p)
    dev.off()
}
data <- data[,-ncol(data)]
data <-data.frame(chr = apply(data,1,function(x){return (unlist(strsplit(x,':'))[1])} ),junction=apply(data[,2:ncol(data)],1,sum))
df <- aggregate(data[,2], by=list(y=data[,1]), FUN=sum)
df <- data.frame(x=df$y, y=df$x)
df <- df[grep("^(chr)?[0-9]+$|^(chr)?X$|^(chr)?Y$", df$x,ignore.case=TRUE),]
#df$x <- factor(paste(df$x),levels=paste(df$x))
df$x<-factor(df$x,levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))
head(df)
#-----------------------------------------------------------------
# plot
#-----------------------------------------------------------------
# mian plot
#p <- qplot(x, y, data=df, geom="bar", stat="identity", fill=x)
p<- ggplot(df, aes(x=x,y=y,fill=x))+geom_bar(stat='identity')
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
p <- p + theme(axis.text.x = element_text(angle = 60, hjust = 1))
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
pdf(file=paste(opt$outfile,"all.chromosome_distrbution.pdf",sep=""), height=opt$height*2/1000, width=opt$width*2/1000)
print(p)
dev.off()
png(filename=paste(opt$outfile,"all.chromosome_distrbution.png",sep=""), height=opt$height, width=opt$width, res=500, units="px")
print(p)
dev.off()




