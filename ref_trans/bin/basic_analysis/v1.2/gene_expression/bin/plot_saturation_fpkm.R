#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript
#library(data.table)
require(RColorBrewer)
library(ggplot2)
library(reshape2)

usage=function(){
  cat("
      Descript: plot pie graph
      Email:    wangmin@biomarker.com.cn
      usage: Rscript plot_sample_fpkm.R  infile=infile.txt outfile=out.png width=140 height=140
      Options:
      -h:          help       NULL
      infile:     character   input file.  [forced]
      outfile:    character   filename for output graph [optional,default:out.png]
      width:      integer     the width of graph, the unit is mm [optional, default:140]
      height:     integer     the height of graph,the unit is mm [optional, default:140]
      title:      character    set the main title of graph
      "
  )
  q(status=1)
}

######---------------------------------------------------------------
#run function
#####----------------------------------------------------------------
myfun=function(infile,outfile="out.png",
               width=140,height=140,
               title=NULL
){
  #-------------------------------------------------
  #read data  and select the final fpkm value >0.1
  #-------------------------------------------------
  d=read.table(infile,header=T)
  d=as.data.frame(d)
  d=d[,-1]
  d=d[d[,10]>0.1,]
  
  #--------------------------------------------------------------
  #calculate the absolution difference percent out of final fpkm
  #--------------------------------------------------------------
  ch=apply(d,1,function(x){
    abs((x-x[10])/x[10])
  })
  ch=t(ch)
  
  #------------------------------------------------------------------
  #the gene number with  diff percnet <0.15 in ervery interval fpkm
  #-----------------------------------------------------------------
  fpkm=d[,10]
  interval=apply(ch,2,function(x){
    fm0=fpkm[x<0.15]   
    
    c(length(which(fm0<=1)),
      length(which(fm0<=5&fm0>1)),
      length(which(fm0<=10&fm0>5)),
      length(which(fm0<=50&fm0>10)),
      length(which(fm0>50))
    )
    
  })
  
  #-------------------------------------------------------
  #the percent gene out of final value
  #-------------------------------------------------------
  interval=apply(interval,1,function(x){x/x[10]})
  colnames(interval)=c("0-1","1-5","5-10","10-50",">50")
  data=data.frame(precent=seq(10,100,10),interval,check.names =F)
  data=melt(data,id=1)
  
  #----------------------------------------------------------------
  #ggplot 
  #---------------------------------------------------------------
  
  size_ratio=width/210
  mytheme=theme_get()
  mytheme$axis.title.x$size=14*size_ratio
  mytheme$axis.title.y$size=14*size_ratio
  mytheme$axis.text.x$size=12*size_ratio
  mytheme$axis.text.y$size=12*size_ratio
  mytheme$legend.title$size=11*size_ratio
  mytheme$legend.text$size=9*size_ratio
  theme_set(mytheme)
  
  p=ggplot(data,aes(x=precent,y=value,colour=variable))+
    geom_line()+
    geom_point(size=2*width/210)+
    ylab("Percent of genes with 15% of final fpkm")+
    xlab("Mapped Reads(%)")+
    scale_colour_manual(values =brewer.pal(5,"Dark2"),name="fpkm interval")
    
    ggsave(filename=outfile,width=width,height=height,unit="mm",dpi=300)
  
}

###################################################################

arg <- commandArgs(T)
#
if(length(arg)==0){
  usage()
}

# ##eval character options to R variable
arg=read.table(text=arg,sep="=",row.names=1,as.is=T)
arg=data.frame(t(arg),stringsAsFactors=F)
arg=as.list(arg)
for(i in names(arg)){arg[[i]]=type.convert(arg[[i]],as.is=T)}

#-----------------------------------------
if("infile" %in% names(arg)){
  do.call(myfun,arg)
}else{
  usage()
}



