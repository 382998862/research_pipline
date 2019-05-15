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
if ( is.null(opt$lab.size ) )		{ opt$lab.size = 18 }
if ( is.null(opt$axis.size ) )		{ opt$axis.size = 12 }
if ( is.null(opt$title.lab) )		{ opt$title.lab = NULL }


#### first
# data <- read.delim(opt$infile, skip=opt$skip,sep ="\t")
# df <- data.frame(x=as.factor(data[,opt$x.col]), group=data[,opt$y.col])
# #gg_normal <-  ggplot(data = df, aes(x = x, fill = group))
# p <-  ggplot(data = df, aes(x=reorder(x,rep(1,length(x)),sum), fill = group))
# p <- p + geom_bar(position = "stack")+ ggtitle("position = stack")
# p <- p + xlab(opt$x.lab) + ylab(opt$y.lab) + labs(title=opt$title)
# p<-p+theme(axis.title.x = element_text(face="bold", size=opt$lab.size),
#            axis.text.x  = element_text(face="bold",color = "black",size=opt$axis.size,angle = 90, hjust = 1),
#            axis.title.y = element_text(face="bold", size=opt$lab.size),
#            axis.text.y  = element_text(face="bold", size=opt$axis.size,color = "black"),
#            title  = element_text(face="bold", size=20,color = "black"))
# p<-p+theme(legend.title = element_text(face="bold", size=opt$lab.size),
#            legend.text = element_text(size=opt$axis.size) )
##############second
# data <- read.delim(opt$infile, skip=opt$skip,sep ="\t")
# df <- data.frame(x=as.factor(data[,opt$x.col]), group=data[,opt$y.col])
# #gg_normal <-  ggplot(data = df, aes(x = x, fill = group))
# p <-  ggplot(data = df, aes(x=reorder(x,rep(1,length(x)),sum), fill = group))
# p <- p + geom_bar(position = "dodge")+ ggtitle("position = dodge")
# p <- p + xlab(opt$x.lab) + ylab(opt$y.lab) + labs(title=opt$title)
# p<-p+theme(axis.title.x = element_text(face="bold", size=opt$lab.size),
#            axis.text.x  = element_text(face="bold",color = "black",size=opt$axis.size,angle = 60, hjust = 1),
#            axis.title.y = element_text(face="bold", size=opt$lab.size),
#            axis.text.y  = element_text(face="bold", size=opt$axis.size,color = "black"),
#            title  = element_text(face="bold", size=20,color = "black"))
# p<-p+theme(legend.title = element_text(face="bold", size=opt$lab.size),
#            legend.text = element_text(size=opt$axis.size) )
####################third
data <- read.delim(opt$infile, skip=opt$skip,sep ="\t")
limit <- max(data[,opt$x.col])
if(limit >20000)
{
  limit<-20000
}
#w=(trunc((3.5*sd(data[,opt$x.col])) / (length(data[,opt$x.col])^(1/3))))/10
#breaks <- pretty(range(data[,3]), n = nclass.FD(data[,3]), min.n = 1)
#w <- breaks[2]-breaks[1]
df <- data.frame(x=data[,opt$x.col], group=data[,opt$y.col])
p <-  ggplot(data = df, aes(x=x, fill = group))
#p <- p +  geom_histogram(position = 'stack')+ ggtitle("position = stack")+scale_x_continuous(breaks = seq(0,limit,by=floor(limit/100)))
p <- p +  geom_histogram(position = 'dodge')+ ggtitle("position = dodge")+scale_x_continuous(limits = c(0,limit),breaks = seq(0,limit,by=ceiling(limit/50)))
p <- p + xlab(opt$x.lab) + ylab(opt$y.lab) + labs(title=opt$title)
p<-p+theme(axis.title.x = element_text(face="bold", size=opt$lab.size),
           axis.text.x  = element_text(face="bold",color = "black",size=opt$axis.size,angle = 90, hjust = 1),
           axis.title.y = element_text(face="bold", size=opt$lab.size),
           axis.text.y  = element_text(face="bold", size=opt$axis.size,color = "black"),
           title  = element_text(face="bold", size=20,color = "black"))
p<-p+theme(legend.title = element_text(face="bold", size=opt$lab.size),
           legend.text = element_text(size=opt$axis.size) )
p <- p + theme( panel.background = element_rect(colour="black", size=1, fill="white"),panel.grid = element_blank())
#-----------------------------------------------------------------
# output plot
#-----------------------------------------------------------------
pdf(file=paste(opt$outfile,".pdf",sep=""), height=opt$height*2/400, width=opt$width*2/400)
print(p)
dev.off()
png(filename=paste(opt$outfile,".png",sep=""), height=opt$height, width=opt$width, res=230, units="px")
print(p)
dev.off()
