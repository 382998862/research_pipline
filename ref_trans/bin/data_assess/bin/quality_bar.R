#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript

#*************************************section 1 usage***************************************

# define usage function
print_usage=function(spec=NULL){
 
  cat("
	discription:
		plot reads average error rate from quality file.see test file.
	Email:    wangmin@biomarker.com.cn
	Options: 
		infile 	 	character 	the input file [forced]
		outfile  	character 	the filename for output graph [optional,quality_bar.png]
		height 	 	integer 	the height of graph,the unit is mm [optional, default: 120]
		width 	 	integer 	the width of graph, the unit is mm [optional, default: 180]
		col 	    character [optional, default: blue]
	example:
		Rscript quality_bar.R infile=T08.quality outfile=T08.quality.png
      \n")
  q(status=1);
}

#******************************************section 2 main function********************************
quality_bar=function(infile,outfile="quality_bar.png",height=80,width=120,col="#4169E1"){
  
  data=read.table(infile,comment.char ="",header=T,check=F,row.names=1)
  #convert Q value to error p value
  Q=as.numeric(colnames(data))
  p=10^((-Q)/10)
  

  #calculate average error p value for every cycle
  err=apply(data,1,function(x){
    mean(x*p/100)
  })
  #y labels
  y.at=seq(0,max(err),by=0.0005)
  if(max(err)>max(y.at)){
    y.at=c(y.at,max(y.at)+0.0005)
  }
  labels=y.at
  labels[2]=0.0005
  labels=formatC(labels,format="f",digits=4)
  labels[1]="0"
  png(filename=outfile,width=width,height=height,units="mm",res=600)
  
  par(mar=c(3,3,2,0.5),mgp=c(1.7,0.1,0),tck=-0.005, cex=0.7,cex.lab=0.8,cex.main=0.8,
      xaxs="i",font.lab=1,font.axis=1)
  out=barplot(err,plot=F,col=col,border="white",axes =F,axisnames=F,space=0)
  
  barplot(err,col=col,border=NA,axes =F,axisnames=F,space=0,
          xlab="",ylab="Average Error Rate",
          main="Reads Average Error Rate",ylim=c(0,max(y.at)),xlim=c(0,max(out[,1])+1),font.lab=2)
  abline(v=1:length(err),col="white",lwd=0.1)
  
  if(length(err)>126){
    out=cbind(out,rep(1:(nrow(data)/2),times=2))
    read=max(out[,1])/2
    segments(read,0,read,5,col="darkgray",lwd=0.5)
    pos=seq(0,nrow(data)/2,50)
    pos[1]=1
    at=out[,2]%in%pos
    axis(1,labels=rep(pos,times=2),at=out[at,1],col="gray30",cex.axis=0.7,font.axis=2,col.axis="gray25",lwd=0.5)
    axis(2,at=y.at,labels=labels,col="gray30",cex.axis=0.6,col.axis="gray25",font.axis=2,lwd=0.5,las=2)
    mtext("Read1",side=1,line=1.5,at=read/2,font=2,cex=0.5)
    mtext("Read2",side=1,line=1.5,at=read*1.5,font=2,cex=0.5)
  }else{
    pos=seq(0,length(err),25)
    pos[1]=1
    axis(1,labels=pos,at=out[pos,1],col="gray30",cex.axis=0.7,col.axis="gray25",font.axis=2,lwd=0.5)
    axis(2,at=y.at,labels=labels,col="gray30",cex.axis=0.6,col.axis="gray25",font.axis=2,lwd=0.5,las=2)
    mtext("Cycle Number",side=1,line=1.5,at=max(out[,1])/2,font=2,cex=0.5)
  }
  
  box(col="gray30",lwd=0.5)
  dev.off()
  
}



#************************section 3 call main function*********************************************** 

arg <- commandArgs(T)
#
if(length(arg)==0) print_usage(spec) 

# ##eval character options to R variable
arg=read.table(text=arg,sep="=",row.names=1,as.is=T)
arg=data.frame(t(arg),stringsAsFactors=F)
arg=as.list(arg)
for(i in names(arg)){arg[[i]]=type.convert(arg[[i]],as.is=T)}

#-------------------------------------------------------------
if(!all(c("infile")%in% names(arg))) {print_usage(spec) }

do.call(quality_bar,arg)
