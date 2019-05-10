#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript

#*************************************section 1 usage***************************************

# define usage function
usage=function(){
  
  cat("
      discription:
      plot actgn distribution.
      Email:    wangmin@biomarker.com.cn
      Options: 
          infile 	 	character 	the input file [forced]
          outfile  	character 	the filename for output graph [ default: <infile>.png]
          height 	 	integer 	the height of graph,the unit is mm [optional, default: 80]
          width 	 	integer 	the width of graph, the unit is mm [optional, default: 100]
      example:
          Rscript plot_acgtn.R infile=T01.acgtn
      \n")
  q(status=1);
}


getarg=function(usage=function(){cat("No argument! \n")}){
  #-------------------------------------------------------------------------
  #Description: getopt is primarily intended to be used with “Rscript”. get arguments from Rscript.
  #             return a list data structure containing names of the flags
  #author:     wangmin@biomarker.com.cn
  #---------------------------------------------------------------------------
  arg <- commandArgs(T)
  #
  if(length(arg)==0){
    do.call(usage,list())
    return(NULL)
  }
  arg=read.table(text=arg,sep="=",row.names=1,as.is=T,colClasses ="character")
  arg=data.frame(t(arg),stringsAsFactors=F)
  arg=as.list(arg)
  # ##eval character options to R variable
  for(i in names(arg)){arg[[i]]=type.convert(arg[[i]],as.is=T)}
  arg
}

#******************************************section 2 main function********************************
plot_acgtn=function(infile,outfile=paste(infile,".png",sep=""),height=80,width=100){
  library(RColorBrewer)
  data=read.delim(infile,check=F)
  
  
  png(filename=outfile,width=width,height=height,units="mm",res=600)
  par(mar=c(3,3,2,0.5),mgp=c(1.5,0.01,0),tck=-0.005, cex=0.7,cex.lab=0.8,cex.main=0.8,
      xaxs="i",font.lab=2,font.axis=2)
  
  matplot(x=data[,1],y=data[,-1],type="l",lty=1,col=c("#008B45","#CD1076","#4169E1","#FF7F24","#636363"),axes=F,
          xlim=c(0,nrow(data)),main="Base Distribution",ylab="Percentage",xlab="",lwd=0.4,font.lab=2,las=2)
  
  
  if(nrow(data)>200){
    out=cbind(data[,1],rep(1:(nrow(data)/2),times=2))
    read=max(out[,1])/2
    abline(v=read,col="darkgray",lwd=0.8)
    pos=seq(0,nrow(data)/2,50)
    pos[1]=1
    at=out[,2]%in%pos
    axis(1,labels=rep(pos,times=2),at=out[at,1],col="gray30",cex.axis=0.6,col.axis="gray25",lwd=0.5)
    axis(2,cex.axis=0.6,col.axis="gray25",lwd=0.5,las=2)
    mtext("Read1",side=1,line=1.5,at=read/2,font=2,cex=0.6)
    mtext("Read2",side=1,line=1.5,at=read*1.5,font=2,cex=0.6)
  }else{
    pos=seq(0,nrow(data),25)
    pos[1]=1
    axis(1,labels=pos,at=pos,col="gray30",cex.axis=0.6,col.axis="gray25",lwd=0.5)
    axis(2,col="gray30",cex.axis=0.6,col.axis="gray25",lwd=0.5,las=2)
    mtext("Cycle Number",side=1,line=1.5,at=max(data[,1])/2,font=2,cex=0.6)
  }
  
  lgd=sub("\\(%\\)","",colnames(data)[-1])
  legend( "topright",lgd,col=c("#008B45","#CD1076","#4169E1","#FF7F24","#636363"),box.lty=1,bty="n",lwd=1,cex=0.5)

  box(col="gray30",lwd=0.5)
  dev.off()
  
  
}


#************************section 3 call main function*********************************************** 
arg=getarg(usage())

if(!all(c("infile")%in% names(arg))) {usage()}

do.call(plot_acgtn,arg)
