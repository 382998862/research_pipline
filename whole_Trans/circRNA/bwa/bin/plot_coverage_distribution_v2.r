#!/share/nas2/genome/biosoft/R/current/bin/Rscript

#####################################################################
# Copyright 2015, BMK
#
# Author:tengh <tengh@biomarker.com.cn>
#
# Function: draw genomewide cytosine coverage distribution map
#
# Modify date: 20150807
#####################################################################

library(getopt)
library(R.utils)

#+--------------------
# get options
#+--------------------
#+--------------------
spec <- matrix(c(
  'help', 'h', 0, "logical", "help",
  'input1', 'i', 1, "character", "input file with preprocessed cytosine coverage depth info, forced.",
  'input2','I',1,"character", "input file with preprocessed cytosine coverage depth info, forced.",
  'output', 'o', 1, "character", "output png file, forced.",
  'window', 'w', 2, "integer", "window size, forced",
  'nchro', 'n', 2, "integer", "Choose the top n chromosome to draw, default [NA]",
  'x.title', 'x', 2, "character", "x title, default [Chromosome position]",
  'y.title', 'y', 2, "character", "y title, default [Median reads density(log2)]",
  'title', 't', 2, "character", "graph title, default [Genomewide distribution of base coverage depth]",
  'col','c',2,"character", "line col,comma separated, or any two of color(Re,Bl,Gr,Pu,Or,Ye,Br,Pi,Gy)combined together, default [BlGr]",
  'bg','b',2,"logical","plot background and gid or not[TRUE]",
  'width','W',2,"integer","graph width, default[according to chromosome names]",
  'height','H',2,"integer","graph height, default[according to chromosome number]",
  'size','s',2,'integer',"fontsize of lab and title text,default[16]",
  'chro','C',2,"character", "select chromosome to draw, default [NA]"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

#+--------------------
# check options
#+--------------------
if ( !is.null(opt$help) | is.null(opt$input1) | is.null(opt$input2)|is.null(opt$output) | is.null(opt$window)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

input1<-as.vector(opt$input1)
if(file.exists(input1)){
  input1=getAbsolutePath(input1)
}else{
  stop("Error:input1 file not exist!")
}
input2<-as.vector(opt$input2)
if(file.exists(input2)){
  input2=getAbsolutePath(input2)
}else{
  stop("Error:input2 file not exist!")
}

#+--------------------
# some default options
#+--------------------
if ( is.null(opt$nchro) ) opt$nchro <- 0
if ( is.null(opt$x.title) ) opt$x.title <- "Chromosome position"
if ( is.null(opt$y.title) ) opt$y.title <- "Median reads density(log2)"
if ( is.null(opt$title) ) opt$title <- "Genomewide distribution of base coverage depth"
if ( is.null(opt$col) ) opt$col <- "BlGr"
if ( is.null(opt$bg) ) opt$bg <- FALSE
if ( is.null(opt$width) ) opt$width <- 5000
if ( is.null(opt$height) ) opt$height <- 4000
if ( is.null(opt$size) ) opt$size <- 16
if ( is.null(opt$chro) ) chrolist <- NA else chrolist=as.vector(opt$chro)
nchro<-as.numeric(opt$nchro )


#+--------------------
# Main
#+--------------------

## load ggplot2 
library(ggplot2)
library(grid)

## get data
reads1 <- read.table(input1, header=FALSE, sep="\t", comment.char = "#")
if (ncol(reads1) != 3) stop("Input file1 format error: must be <chromosome> <position> <coverage>")
reads1<-reads1[reads1[,3]>0,]
reads1<-data.frame(reads1,flag=1,min=0)
colnames(reads1) <- c("chromosome", "position", "coverage","flag","min")


reads2<-read.table(input2, header=FALSE, sep="\t", comment.char = "#")
if (ncol(reads2) != 3) stop("Input file format error: must be <chromosome> <position> <coverage>")
reads2<-reads2[reads2[,3]>0,]
reads2<-data.frame(reads2,flag=2,min=0)
colnames(reads2) <- c("chromosome", "position", "coverage","flag","min")

reads2$min=-reads2$coverage
reads2$coverage<-0
reads<-rbind(reads1,reads2)
chromosome<-as.vector(reads$chromosome)
chro<-unique(chromosome)
#gregexpr("(\\d*.*\\d+)$", "Maize_23newGene_8", perl = TRUE)
rowname_parts <- function(x) {
  x <-as.character(x)
  m <- regexec("(\\d*\\.*\\d+)$", x)
  parts_1<-as.vector(do.call(rbind,lapply(x,function(x)gsub("(\\d*\\.*\\d+)$","",x))))
  parts_2 <-as.numeric(as.vector(do.call(rbind,lapply(regmatches(x, m), `[`, 1L))))
  result<-data.frame(x,name=parts_1,id=parts_2) 
}
chro<-rowname_parts(chro)
chro<-chro[order(chro[,2],chro[,3]),]
chro<-as.vector(chro[,1])
chro1<-unique(as.vector(reads1[,1]))
chro2<-unique(as.vector(reads2[,1]))

if(is.na(chrolist)||!file.exists(chrolist)){
  #check whether the input1 input2 chromosome is identity or not
  if(!setequal(chro1,chro2)){
    if(sum(!is.element(chro1,chro2))>0)message(paste(c(setdiff(chro1,chro2),"not exist in input file2:",input1,"!"),collapse = "\t"))
    if(sum(!is.element(chro2,chro1))>0)message(paste(c(setdiff(chro2,chro1),"not exist in input file1:",input2,"!"),collapse = "\t"))
    #chro=intersect(chro1,chro2)
    #reads<-reads[is.element(chro1,chro),]
    #chro1=chro
    #chromosome<-chromosome[is.element(chromosome,chro)]
  }
  ## subset data frame if chromosome in scaffold level and specified the top_n option
  if(nchro<1){
    nchro=ifelse(length(chro)>20,10,length(chro))
  }
  if (nchro>0) {
    if(nchro>length(chro)){
      message("Coverage distribution parameter: Set the nchro to the numbr of chromosome!")
      nchro=lenght(chro)
    }
    if(opt$window<1000){
      chro2<-tapply(reads[,"position"],reads[,"chromosome"],length)
      chro2<-names(chro2)[order(-chro2)]
      chro_kept<-intersect(chro,chro2[1:nchro])
    }else{
      chro_kept <-chro[1:nchro]      
    }
    message("chro list:")
    message(paste(chro_kept,collapse=";"))
    reads <- reads[chromosome %in% chro_kept,]
    reads$chromosome<-factor(reads$chromosome,levels = chro_kept)
    reads<-reads[order(reads$chromosome),]
  }
}else{
  #输入文件列表
  chrolist<-as.matrix(read.table(chrolist, header=FALSE, sep="\t", comment.char = "#"))[,1]
  chro<-as.vector(chrolist)
  if(sum(!is.element(chro,chro1))>0)message(paste(c(setdiff(chro,chro1),"not exist in input file1!\n",input1," "),collapse = ""))
  if(sum(!is.element(chro,chro2))>0)message(paste(c(setdiff(chro,chro2),"not exist in input file2!\n",input2," "),collapse = ""))
  reads<-reads[is.element(chromosome,chro),]
  nchro=length(chro)  
}
scale=ifelse(opt$window>1000,1e6,1e3)
names(scale)<-ifelse(opt$window>1000,"M","K")
chro<-unique(as.vector(reads$chromosome))
#plot the alignment map
reads$flag=factor(reads$flag,levels = (c("1","2")))
reads$chromosome=factor(reads$chromosome,levels =chro)
## plot 
p<-ggplot(reads, aes(x=position*as.numeric(opt$window)/scale, ymax=coverage,ymin=min,colour=flag))+geom_linerange()## create plot
p <- p + facet_grid(chromosome~., space="fixed") ## facet by chromosome
p <- p + theme(legend.position="none")
x_max <- round(max(reads$position * as.numeric(opt$window)/scale), digits=0)
x_step <- round((round(x_max / 5) / 10+0.5)) * 10
x_breaks <- seq(0, x_max, x_step)
p <- p + scale_x_continuous(breaks=x_breaks, labels=paste(x_breaks, names(scale), sep=""))  ## set x and y axis text 
ymax<-round(0.8*max(reads$coverage))
p <- p + scale_y_continuous(breaks=c(-ymax,0,ymax), labels=c(-ymax,0,ymax))  ## set x and y axis text 

p <- p + theme(strip.background=element_rect(fill="white", colour="white"), strip.text.y=element_text(angle=360, vjust=0.5), 
               axis.line=element_line(colour="black"),text = element_text(size=opt$size, family="serif",face="bold"),
               axis.text = element_text(colour="black",size=10,face="plain")) ## set other theme opts 
p <- p + labs(title=opt$title, x=opt$x.title, y=opt$y.title) ## set x and y labs

#set col
if(!is.na(opt$col)&&opt$col!=""){
  opt$col=unlist(strsplit(opt$col,split = ","))
  if(length(opt$col)==1){
    col<-c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999")
    names(col)<-c("Re","Bl","Gr","Pu","Or","Ye","Br","Pi","Gy")
    nchar<-nchar(opt$col)
    if(nchar>=4){
      col1=substr(opt$col,1,2)
      col2=substr(opt$col,3,4)
      if(sum(!is.element(c(col1,col2),names(col)))>0){
        col=col[c("Bl","Gr")]
      }else{
        col=col[c(col1,col2)]
      }
    }else{
      col=col[c("Bl","Gr")]
    }
  }else{
    col=opt$col[1:2]
  }
  names(col)<-c("1","2")
  p<-p+scale_color_manual(values=col)
}

if(!opt$bg){
  p<-p+theme(axis.line = element_line(colour = "black"),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.border = element_blank(),
             panel.background = element_rect(fill="white", colour="white")) 
}
if(is.na(opt$height)){
  opt$height=ifelse(nchro*500+500>4500,4500,nchro*220+400)
}
if(is.na(opt$width)){
  opt$width=4000+ifelse(max(nchar(chro))*33>1000,1000,max(nchar(chro))*33)
}
## output to png file 
ggsave(filename=opt$output, p,height =opt$height/500, width = opt$width/500, dpi = 500, units = "in")
