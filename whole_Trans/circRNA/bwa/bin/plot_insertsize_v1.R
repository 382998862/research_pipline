#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript

source("/share/nas1/wangm/yf/plot/gplot/gplot_v2.R")

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
