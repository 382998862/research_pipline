#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript


#*************************************section 1 usage

usage=function(){
  cat("
      Descript: plot smooth randcheck graph
      Email:    wangmin@biomarker.com.cn
      Options:
      infile:     character    input file.  [forced]
      col:        character     the colour of line [default:red]
      width:      integer      the width of graph, the unit is mm [optional, default:120]
      height:     integer      the height of graph,the unit is mm [optional, default: 60]
      bg          logical      whether use background, the value is T or F, T use background, F not [default: T]
      example
      Rscript plot_random.R infile=Total.randcheck.list
      Rscript plot_random.R infile=Total.randcheck.list  bg=F \n"
  )
  q(status=1)
}

#*************************************section 2 main function**************************************
plot_random=function(infile,width=120,height=60,col="red",bg=T){
  library(ggplot2)
  library(grid)
  require(RColorBrewer)
  data=read.table(infile,header=T,sep="\t",check=F)
  colnames(data)=c("Sample", "x", "y")
  #ylim
  ylim=c(0,1.5)
  if(max(data$y)>1.45){ylim=c(0,2)}
  
  #theme
  if(bg){
    mytheme=theme_grey()
  }else{
    mytheme=theme_bw()
    mytheme$panel.grid.minor$colour="white"
    mytheme$panel.grid.major$colour="white"
    mytheme$panel.border$colour="gray25"
  }
  size_ratio=width/210
  mytheme$axis.title.x$size=12*size_ratio
  mytheme$axis.title.y$size=12*size_ratio
  mytheme$axis.text.x$size=10*size_ratio
  mytheme$axis.text.y$size=10*size_ratio
  mytheme=mytheme+ theme(legend.title = element_text(size=8*size_ratio),
                         legend.text = element_text(size=8*size_ratio),
                         legend.key.size=unit(width*0.03,"mm"),
                         legend.margin=unit(0,"mm"),
                         plot.margin=unit(c(0.5,0,0,0),"lines"))
  theme_set(mytheme)
  
  p=ggplot(data, aes(x=x, y=y, colour=Sample))+
    geom_smooth(size=0.4,se=F,level=0.99,span=0.1)+
    labs(colour="Sample") + 
    xlab("Relative Position in Genes(5'-3')")+
    ylab("Percent of Reads")+
    ylim(ylim)
  ggsave(filename="Total.randcheck.png",p,width=width,height=height,units="mm",dpi=300)
  
  
  data$Sample=as.character(data$Sample)
  
  for(i in unique(data$Sample)){
    d=subset(data,Sample==i)
    p=ggplot(data, aes(x=x, y=y))+
      geom_smooth(size=0.6,se=F,level=0.99,span=0.1,colour=col)+
      labs(colour="Sample") + 
      xlab("Relative Position in Genes(5'-3')")+
      ylab("Percent of Reads")+
      ylim(ylim)
    ggsave(filename=paste(i,".randcheck.png",sep=""),p,width=width,height=height,units="mm",dpi=300)
  }
  
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

do.call(plot_random,arg)




