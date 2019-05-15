#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript

usage=function(){
  cat("
  Descript: plot pie graph
  Email:    wangmin@biomarker.com.cn
  usage: Rscript pie.R  infile=value1 outfile=value2 value.col=value3 legend.col=value4 title=value5
  Options:
        -h:          help       NULL
        infile:     character   the input file,Note:default header=F [forced]
        skip:       integer     the number of line for skipping [optional, default: 0]
        sep         character 	the field separator character for infile,tab is symbol as t. default is blank
        outfile     character 	the filename for output graph [default: pie.png]
        legend.col: integer     the column which is used to display legend,default is 1
        value.col:  integer     the column which is used to display pie chart, default is 2              
        title:      charcter    the title of the graph, default is blank
        color:      charcter    A palette name pass to function brewer.pal() in package RColorBrewer, default is Set1
 example: 
  Rscript pie.R  infile=test.txt outfile=out_pie.png legend.col=1  value.col=3 skip=1 sep=t   
    "
  )
  q(status=1)
}
mypie=function(infile,skip=0,sep=" ",outfile="pie.png",legend.col=1,value.col=2,title="",color="Set1"){
  
  library(RColorBrewer)
  data=read.table(infile,sep=sep,skip=skip,as.is=T,comment="")
  
for(i in value.col){
  #data=data[order(data[,i], decreasing =T),]
  data[,i]<-as.numeric(data[,i])
  labels=data[,i]/sum(data[,i])*100
  labels=paste( round(labels,digits=1),"%",sep="")
  #legend=paste(data[,legend.col],": ",data[,i],sep="")
  legend=data[,legend.col]
  rmar=(max(nchar(legend)))*0.4
  if(rmar<0.5){rmar=0.5}
 # ofl=paste(outfile,colnames(data)[i],"png",sep=".")
  
  png(filename=outfile,width=115*(1+rmar*0.05),height=100,units="mm",res=600)
  par(mar=c(0.5,0.5,0.5,rmar))
  col=brewer.pal(nrow(data),color)
  
  pie(data[,i],labels=labels,col=col,
      border="gray",main=title,init.angle=-20)
  
  legend(1.05,0,legend,fill=col,col=col,box.col ="transparent",bg="transparent",
         xpd=TRUE,x.intersp=0.05,yjust=0.5)
  dev.off()
}
  
  
}
##start run
arg <- commandArgs(T)

if(length(arg)==0){usage()}

arg=read.table(text=arg,sep="=",row.names=1,as.is=T)
arg=data.frame(t(arg),stringsAsFactors=F)
arg=as.list(arg)

if(!"infile" %in% names(arg)){cat("\n Error: infile must be specified \n");q(status=1)}



for(i in names(arg)){arg[[i]]=type.convert(arg[[i]],as.is=T)}

if("sep" %in% names(arg)){
  if(arg[["sep"]]=="t"){arg[["sep"]]="\t"}
}

do.call(mypie,arg)


