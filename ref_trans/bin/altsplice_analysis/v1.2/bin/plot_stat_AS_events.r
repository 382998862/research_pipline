
#!/share/nas2/genome/biosoft/R/current/bin/Rscript

#####################################sub function
facetAdjust <- function(x, pos = c("up", "down")){
  pos <- match.arg(pos)
  p <- ggplot_build(x)
  gtable <- ggplot_gtable(p); dev.off()
  dims <- apply(p$panel$layout[2:3], 2, max)
  nrow <- dims[1]
  ncol <- dims[2]
  panels <- sum(grepl("panel", names(gtable$grobs)))
  space <- ncol * nrow
  n <- space - panels
  if(panels != space){
    idx <- (space - ncol - n + 1):(space - ncol)
    gtable$grobs[paste0("axis_b",idx)] <- list(gtable$grobs[[paste0("axis_b",panels)]])
    if(pos == "down"){
      rows <- grep(paste0("axis_b\\-[", idx[1], "-", idx[n], "]"), 
                   gtable$layout$name)
      lastAxis <- grep(paste0("axis_b\\-", panels), gtable$layout$name)
      gtable$layout[rows, c("t","b")] <- gtable$layout[lastAxis, c("t")]
    }
  }
  class(gtable) <- c("facetAdjust", "gtable", "ggplot"); gtable
}

print.facetAdjust <- function(x, newpage = is.null(vp), vp = NULL) {
  if(newpage)grid.newpage()
  if(is.null(vp)){
    grid.draw(x)
  } else {
    if (is.character(vp)) 
      seekViewport(vp)
    else pushViewport(vp)
    grid.draw(x)
    upViewport()
  }
  invisible(x)
}

barMap<-function(bardata,bg=F,size,x.title,y.title,title,output,width,height,ncol=3){
  
  ## plot 
  p<-ggplot(bardata, aes(x=ASType, y=number,fill=ASType))+geom_bar(stat="identity")## create plot
  p<-p+coord_flip()
  p<-p+facet_wrap( ~ sample, ncol = 3) 
  p <- p + theme(legend.position="none") 
  p <- p + theme(text = element_text(size=size, family="serif"),
                 axis.text = element_text(colour="black",size=10),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank()) ## set other theme opts 
  p <- p + labs(title=title, x=x.title, y=y.title) ## set labs 
  if(bg){
    p<-p+theme(axis.line=element_line(colour="black"),
               strip.background=element_rect(fill="white", colour="white"),
               strip.text.y=element_text(angle=360, vjust=0.5),
               axis.line = element_line(colour = "black"),
               panel.border = element_blank(),
               panel.background = element_blank(),
               panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }
  p <- facetAdjust(p,pos="down")
  ggsave(filename=output, p,height = height/500, width = width/500, dpi = 500, units = "in")
}
###################################################
#####################################################################
# Copyright 2014, BMK
#
# Author: tengh <tengh@biomarker.com.cn>
#
# Function: draw statics of AS events map
#
#####################################################################

library(getopt)
library(R.utils)
#+--------------------
# get options
#+--------------------
spec <- matrix(c(
  'help', 'h', 0, "logical", "help",
  'input','i',1,"character", "input file with preprocessed cytosine coverage depth info, forced.",
  'output', 'o', 1, "character", "output png file, forced.",
  'x.title', 'x', 2, "character", "x title, default [Number].",
  'y.title', 'y', 2, "character", "y title, default [AS type].",
  'title', 't', 2, "character", "graph title, default [Statistics of AS events]",
  'size','s',2,'integer',"fontsize of lab and title text,default[16]",
  'bg','b',2,"logical","plot background or not[FALSE]",
  'width','W',2,"integer","graph width default[4000]",
  'height','H',2,"integer","graph height default[3000]"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)
#+--------------------
# check options
#+--------------------
if ( !is.null(opt$help) | is.null(opt$input) |is.null(opt$output)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

input<-as.vector(opt$input)
if(file.exists(input)){
  input=getAbsolutePath(input)
}else{
  stop("Error:input file not exist!")
}


#+--------------------
# some default options
#+--------------------

if ( is.null(opt$x.title) ) opt$x.title <- "AS type"
if ( is.null(opt$y.title) ) opt$y.title <- "Number"
if ( is.null(opt$title) ) opt$title <- ""
if ( is.null(opt$bg) ) opt$bg <- FALSE
if ( is.null(opt$width) ) opt$width <- 4000
if ( is.null(opt$height) ) opt$height <- 3000
if ( is.null(opt$size) ) opt$size <- 15
#+--------------------
# Main
#+--------------------

## load ggplot2 
library(ggplot2)
library(grid)
#opt<-data.frame(input="As_Result.stat" ,output="As_event_stat",stringsAsFactors=F)
## get data
dat <- read.table(input, header=T, sep="\t", comment.char = "",check.names =F)
colnames(dat)<-read.delim(input,header=F,check.names =F,stringsAsFactors=F,nrows=1)
rownames(dat)<-as.vector(dat[,1])
#dat<-as.matrix(dat[,-1])
sample<-colnames(dat)[-1]
bardat<-data.frame(sample=rep(colnames(dat)[-1],each=nrow(dat)),ASType=rep(rownames(dat),ncol(dat)-1),number=as.vector(as.matrix(dat[,-1])))
ASType<-rownames(dat)
bardat$ASType=factor(bardat$ASType,levels = ASType)
for(i in seq(1,length(sample),by=6)){
  samID<-ifelse(i+5>length(sample),length(sample),i+5)
  if(samID-i<2){
    ncol=samID-i+1
  }else if(samID-i==3){
    ncol=2
  }else{
    ncol=3
  }  
  output=paste(opt$output,i,"_",samID,".png",sep="")
  samID<-sample[i:samID]
  temp<-bardat[is.element(bardat$sample,samID),]
  temp$sample<-factor(temp$sample,levels=samID)
  png(file = output) # NOT create unnecessary Rplots.pdf file
  barMap(bardata = temp,bg=opt$bg,size=opt$size,x.title=opt$x.title,y.title=opt$y.title,title=opt$title,output=output,width=opt$width,height=opt$height,ncol=ncol)  
}


