#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript


# load library
library('getopt');
#source("/share/nas1/tengh/research/pipelines/source/pieplot.R")
#opt<-data.frame(outfile="pieplot",infile="E:/R_workplace/20150624pieplot/AllSample.data.stat",title="",height=3000,width=4000,legend.pos="bottom",col="Accent",lab.size=5,legend.size=10,legend.ncol=2,bg=NA)

#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'help' , 'h', 0, "logical",
  'infile' , 'i', 1, "character",
  'outfile' , 'o', 0, "character",
  'key' , 'k', 0, "character",
  'title' , 't', 0, "character",
  'height' , 'e', 0, "integer",
  'width' , 'w', 0, "integer",
  'legend.pos' , 'P', 0, "character",
  'legend.ncol' , 'N', 0, "integer",
  'lab.size' , 'l',0, "integer",
  'legend.size' , 'S', 0, "integer",
  'col','c',0,"character",
  'bg' , 'b', 0, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

# define usage function
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
      Usage example: 
      1) Rscript pie_plot.R --infile in_Pie.data --outfile Pie \\
      --col Accent
      2) Rscript pie_plot.R --infile in_Pie.data --outfile Pie \\
      --legend.pos bottom --legend.ncol 2
      
      Options: 
      --help		-h 	NULL 		get this help
      --infile 	-i 	character 	the input file [forced]
      --outfile 	-o 	character 	the filename for output graph [optional]
      --key 	-k 	character 	the  key for output graph [optional]
      --height 	-h 	integer 	the height of graph [optional, default: 3000]
      --width 	-w 	integer 	the width of graph [optional, default: 4000]
      --legend.pos 	-P 	character 		Set legend position[optional:right]
      --legend.ncol 	-N 	integer 	the col number for legend disply [optional, default: 1]
      --title 	-t 	character 	the lab for title [optional, default: '']
      --lab.size 	-l 	integer 	the font size of lab [optional, default: 5]
      --legend.size 	-S 	integer 	the font size of text for legend [optional, default: 10]
      --bg	-b 	character 		set the backgroud [optional, default: NA]
      --col 		-c 	character 	choose the colour set(Set1,Set2,Set3,Accent,Dark2,Pastel1,Pastel2,Paired)or set colour splited by , [optional, default: Set3]
      \n")
  q(status=1);
}


# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }


# check non-null args
if ( is.null(opt$infile) )	{ print_usage(spec) }else{opt$infile<-gsub("\\\\",replacement = "/",x = opt$infile)}

#set some reasonable defaults for the options that are needed,
#but were not specified.

if ( is.null(opt$outfile) )  { opt$outfile =paste(c(getwd(),"/","pie_plot"),collapse = "")}
if ( is.null(opt$height ) )		{ opt$height = 3000 }
if ( is.null(opt$width ) )		{ opt$width = 4000 }
if ( is.null(opt$title) )		{ opt$title="" }
if ( is.null(opt$legend.pos) )		{ opt$legend.pos ="right" }
if ( is.null(opt$lab.size ) )		{ opt$lab.size = 5}
if ( is.null(opt$legend.size ) )	{ opt$legend.size = 10 }
if ( is.null(opt$legend.ncol ) )	{ opt$legend.ncol = 1 }
if ( is.null(opt$col ) )  { opt$col ="Accent" }
if ( is.null(opt$bg ) )  { opt$bg= NA}
col=strsplit(x =as.vector(opt$col),split = "," )[[1]]

#factor(a$name, levels=unique(a$name))
#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------
# reading data

data <- read.delim(opt$infile,sep="\t",header = T)
if ( is.null(opt$key ) )  {
	if (is.null(dim(data)) ){
		 sample<-as.vector(data[1])
		 pie.data<-as.matrix(data[2:5],nrow=1)  
	}else{
		sample<-as.vector(data[,1])
		pie.data<-as.matrix(data[,2:5])   
	}	
}else{
	sample<-read.delim(opt$key,sep="\t",header = F)
	sample<-as.vector(t(sample[,1]))
	pie.data<-as.matrix(data[,2:5])
}

# check dim
#if (is.null(dim(data)) ){
#  sample<-as.vector(data[2])
#  pie.data<-as.matrix(data[3:11],nrow=1)  
#}else{
#  sample<-as.vector(data[,2])
#  pie.data<-as.matrix(data[,3:11])   
#}
pie.data<-data.frame(Adapter=pie.data[,grepl("^Adapter",x = colnames(pie.data))],Low_quality=pie.data[,grepl("^Inferior",x = colnames(pie.data))])
pie.data$Clean<-100-pie.data$Adapter-pie.data$Low_quality
pie.data<-as.matrix(pie.data,ncol=3)
colnames(pie.data)<-c('Adapter Related',"Low Quality","Clean Reads")
row.names(pie.data)<-sample
#plot pie charts
pieplot<-function(x,legend=names(x),labels=NULL,title="",outfile="pie.plot",height=3000,width=4000,legend.pos="right",col="Set1",is.axis=F,lab.size=5,legend.size=10,no.grid=T,legend.ncol=1,bg=NA){  
  library(ggplot2)
  require(RColorBrewer)
  library(grid)
  df<-data.frame(y=x,x=legend)
  coln<-c(8,8,12,9,8,9,8,12)
  names(coln)<-c("Accent","Dark2","Paired","Pastel1","Pastel2","Set1","Set2","Set3")
  if(length(col)>1&&length(x)!=length(col))col="Set3"
  if(length(col)==1){
    if(is.na(col)||is.na(coln[col])){
      col="Set3"
      message("Select colourSet 'Set3'!")
    }    
    col<-rep(brewer.pal(coln[col],col),length.out =length(x))
  }
   col<-c("#ff0000","#0000ff","#339900","#ffcc00","#00ccff")
  
  
  #-----------------------------------------------------------------
  # plot
  #-----------------------------------------------------------------
  # mian plot
  if((!is.null(labels))&&length(labels)==length(x)){
    j=0;
    xat=numeric(0)
    yat=numeric(0)
    for(i in 1:length(x)){
      if(x[i]<=0){
        xat=c(xat,1.6)
        yat=c(yat,(df$y/2 + c(0, cumsum(df$y)[-length(df$y)]))[i])
        labels[i]="";
      }else if(x[i]/sum(x)<0.1){
        if(i==0&&j==0)j=-min(cumsum(df$y)[-length(df$y)])
        xat=c(xat,1.6)
        yat=c(yat,(df$y/2 + c(0, cumsum(df$y)[-length(df$y)]))[i]+j)
        j=j+nchar(labels[i])/2*lab.size*0.3;
      }else{ 
        xat=c(xat,0.9)
        yat=c(yat,(df$y/2 + c(0, cumsum(df$y)[-length(df$y)]))[i])
        j=0
      }
    }    
  }
  df$xat=xat
  df$yat=yat
  df$labels=labels
  p <- ggplot(df, aes(x="", y=y, fill=x,)) + geom_bar(width = 1, stat="identity",alpha=0.5)
  p <- p + scale_fill_manual(values = col)
  
  if((!is.null(labels))&&length(labels)==length(x)){
    p=p + geom_text(aes(x=xat, y =yat, label=labels), size=lab.size)
  }
  
  # coord polar
  p <- p + coord_polar(theta = "y")
  
  #-----------------------------------------------------------------
  # theme
  #-----------------------------------------------------------------
  # lab
  p <- p +labs(title=title,x="",y="")+xlab("")
  # set lab and axis test size
  p <- p + theme(title = element_text(face="bold", size=lab.size), 
                 axis.text = element_text(face="bold", size=lab.size))
  p <- p + theme(legend.title = element_text(face="bold", size=legend.size),
                 legend.text = element_text(size=legend.size) ,legend.position=legend.pos)
  p
  #delete
  if( !is.axis ){
    p <- p + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                   panel.grid  = element_blank())
  }
  
  # legend col
  if( is.null(legend.ncol) ){
    legend.ncol=1
  }
  p <- p + guides(fill = guide_legend(ncol = legend.ncol,title=""))
  # grid and background
  if ( !is.na(bg)) {
    p <- p + theme( panel.background = element_rect(colour=bg, size=1, fill=bg),
                    panel.grid = element_blank())
  }
  p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
 p<-p+ theme(panel.background = element_blank(),panel.border = element_blank(),panel.grid = element_blank(),axis.text = element_blank())

  #-----------------------------------------------------------------
  # output plot
  #-----------------------------------------------------------------
  #pdf(file=paste(outfile,".pdf",sep=""), height=opt$height*2/1000, width=opt$width*2/1000)
  #print(p)
  #dev.off()
  png(filename=paste(outfile,".png",sep=""), height=opt$height, width=opt$width, res=500, units="px")
  print(p)
  dev.off()
}

for(i in 1:nrow(pie.data)){
  x<-pie.data[i,]
  legend<-colnames(pie.data)
  outfile<-paste(as.vector(opt$outfile),sample[i],sep="_")
  labels=paste(x,"%",sep="")
  pieplot(x = x,legend = legend,lab.size=4,labels = labels,title =as.vector(opt$title),outfile =outfile,height = opt$height,width = opt$width,legend.pos = as.vector(opt$legend.pos),col = col,legend.ncol = opt$legend.ncol,legend.size = opt$legend.size,bg = opt$bg)
}






