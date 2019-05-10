#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript

gplot=function(infile,outfile="plot.png",sep="",skip=0,geom1="geom_point",geom2=NULL,
               x=1,y=2,width=120,height=90,col=brewer.pal(9,"Set1")[1],
               xlim=NULL,ylim=NULL,xlab=NULL,ylab=NULL,title=NULL,bg=F){
# #----------------------------------------------------------------------------------------------------
#       Descript: simple general ggplot2 
#       
#       Options:
#       infile:     character   input file.  [forced]
#       skip:       integer     the number of line for skipping [optional, default: 0]
#       sep         character 	the field separator character for infile. default is ¡®white space¡¯
#       outfile:    character   filename for output graph [forced]
#       geom1:      character   the graph function,see ggplot2 package. default is geom_point. more counld be
#       geom_line,geom_smooth, geom_path.
#       geom2:     character    the same as geom1. default geom2=NULL
#       col:       character     the colour of geom1
#       x:         integer      the column for x value. default x=1
#       y:         integer      the column for x value. default y=2
#       width:      integer     the width of graph, the unit is mm [optional, default:120]
#       height:     integer     the height of graph,the unit is mm [optional, default: 90]
#       xlim:       vector       set the limits of the x axis.
#       ylim:       vector       set the limits of the y axis.
#       xlab:       character    set the  x axis legend titles.
#       ylab:       character    set the  y axis legend titles.
#       title:      character    set the main title of graph
#       bg          logical      whether use background, the value is T or F, T use background, F not [default: F]
#------------------------------------------------------------------------------------------------------------------
    library(ggplot2)
    library(grid)
    require(RColorBrewer)
  
  #read data----------------------------------
  data=read.table(file=infile,sep=sep,skip=skip,as.is=T,comment="")
  
  #plot par 
  plotpar=list()
  if(grepl("point",geom1)){
    plotpar$size=3*width/210
  }else{
    plotpar$size=1.3*width/210
  }
  
  plotpar$colour=col
  
  if(grepl("smooth",geom1)){
    plotpar$se=F
    plotpar$level=0.99
    plotpar$span=0.1
  }
  
  
  #theme
    if(bg){
      mytheme=theme_get()
    }else{
      mytheme=theme_bw()
      mytheme$panel.grid.minor$colour="white"
      mytheme$panel.grid.major$colour="white"
      mytheme$panel.border$colour="black"
    }

  size_ratio=width/210
  mytheme$axis.title.x$size=14*size_ratio
  mytheme$axis.title.y$size=14*size_ratio
  mytheme$axis.text.x$size=12*size_ratio
  mytheme$axis.text.y$size=12*size_ratio
  theme_set(mytheme)
  
  #-------------------------------------
  colnames(data)[c(x,y)]=c("x","y")
  p=ggplot(data,aes(x,y))+
    do.call(geom1,plotpar)
  
  if(!is.null(geom2)){p=p+ do.call(geom2,list())}
  
  if(!is.null(xlim)){p=p+xlim(xlim)}
  if(!is.null(ylim)){p=p+ylim(ylim)}
  if(!is.null(xlab)){p=p+xlab(xlab)}
  if(!is.null(ylab)){p=p+ylab(ylab)}
  if(!is.null(title)){p=p+ggtitle(title)}
  ggsave(filename=outfile,p,width=width,height=height,units="mm",dpi=600)
}
#*************************************section 1 usage

usage=function(){
  cat("
      Descript: plot smooth insert size graph
      Email:    wangmin@biomarker.com.cn
      Options:
      infile:     character   input file whith tab and no header.  [forced]
      outfile:    character   filename for output graph [default: Inset_Size_v1.png]
      col:        character     the colour of line {default: red}
      width:      integer     the width of graph, the unit is mm [optional, default:120]
      height:     integer     the height of graph,the unit is mm [optional, default: 90]
      bg          logical      whether use background, the value is T or F, T use background, F not [default: F]
      example:
      Rscript plot_insertsize_v1.R   infile=T05.insertSize.r.list
      "
  )
  q(status=1)
}

#----------------------------------------------------------------------------------------------------

mergelist=function(...){
  if (missing(...)) return(invisible(NULL))
  dots=list(...)
  List=list()
  for(i in 1:length(dots)){
    for(j in names(dots[[i]])){
      List[[j]]=dots[[i]][[j]]
    }
  }
  List
}

#************************section 3 call main function*********************************************** 

arg <- commandArgs(T)
#
if(length(arg)==0) usage()

# ##eval character options to R variable
arg=read.table(text=arg,sep="=",row.names=1,as.is=T)
arg=data.frame(t(arg),stringsAsFactors=F)
arg=as.list(arg)
for(i in names(arg)){arg[[i]]=type.convert(arg[[i]],as.is=T)}

#-------------------------------------------------------------
if(!all(c("infile")%in% names(arg))) {usage()}

default=list(outfile="Inset_Size_V1.png",sep="\t",skip=0,geom1="geom_smooth",geom2=NULL,
x=1,y=2,width=120,height=90,
xlim=NULL,ylim=NULL,xlab="Inset Size(bp)",ylab="Number of Reads",title=NULL,bg=F)

arg=mergelist(default,arg)

do.call(gplot,arg)
