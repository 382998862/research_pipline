args<-commandArgs(TRUE)
usage <- function() {
        print("-------------------------------------------------------------------------------")
        print("Usage: Rscript vn.r  input outpath keyword")
        print("1) input: multiple files seperated by comma or one file contained multiple columns")
        print("2) outpath: output path")
        print("3) keyword: out png/pdf keyword")
        print("-------------------------------------------------------------------------------")
        q(status=1);
}

if(length(args)!=3){
        usage()
}

library(Vennerable)
library(venn,lib.loc='/share/nas1/liubo/qing/app_pipeline/Metabolome_analysis/v1.1/bin/Diff_analysis/venn_lib')
library(grid)
library(VennDiagram)
lists<-vector(mode="list")
if(file.exists(args[1])==TRUE){
	read.csv(args[1],sep="\t",header=T)->data
	for(i in 1:dim(data)[2]){
		lists[[i]]<-unique(sort(as.character(data[,i])))
	}
	Vstem <- Venn(lists,SetNames=colnames(data))
	num=length(colnames(data))
}

zero<-0
if(file.exists(args[1])==FALSE){
	files=strsplit(args[1],split=",")[[1]]
	bases<-vector(length=length(files))
	for(i in 1:length(files)){
		tmp<-strsplit(files[i],split="/")[[1]]
		read.csv(files[i],sep="\t",header=T,row.names=1)->f
		genes<-row.names(f)
		len <-length(genes)
		if(len==0){zero<-1}
		lists[[i]]<-genes
		bases[i]<-strsplit(tmp[length(tmp)],split="[.]")[[1]][1]
	}
	
	Vstem <- Venn(lists,SetNames=bases)
	num<-length(files)
}
print(zero)
names(lists)<-bases

if(length(bases)>=2 & length(bases)<=3){
	print("2-3 sets:Vennerable plot")
	png(paste(args[2],"/",args[3],".png",sep=""),width=800,height=800)
	plot(Vstem, type="circles",doEuler=F,doWeight=F)
	dev.off()
	pdf(paste(args[2],"/",args[3],".pdf",sep=""),width=16,height=16)
	plot(Vstem, type="circles",doEuler=F,doWeight=F)
	dev.off()
}

if(length(bases)>=4 && length(bases)<=5){
	if(zero==2){
		print("4-5 sets:venn plot")
		png(paste(args[2],"/",args[3],".png",sep=""))
		venn(lists,zcolor='style')
		dev.off()
		pdf(paste(args[2],"/",args[3],".pdf",sep=""))
		venn(lists,zcolor='style')
		dev.off()
	}else{
		print("4-5 sets:VennDiagram plot")
		if(length(bases)==4){
			col=c('cornflowerblue','green','yellow','darkorchid1')
		}else{
			col=c('cornflowerblue','green','yellow','darkorchid1','red')
		}
		pic_png=paste(args[2],"/",args[3],".png",sep="")
		pic_pdf=paste(args[2],"/",args[3],".pdf",sep="")
		pic_svg=paste(args[2],"/",args[3],".svg",sep="")
		venn.diagram(lists,filename=pic_png,imagetype="png",margin=0.25,height=900,width=900,resolution=200,units='px',lwd=1,fill=col,cex=0.6,cat.cex=0.8,scaled=0)
		venn.diagram(lists,filename=pic_svg,imagetype="svg",margin=0.25,height=12,width=12,resolution=500,lwd=1,fill=col,cex=1,cat.cex=1,scaled=0)
		shell_cmd<-paste("convert ",pic_svg," ",pic_pdf)
		print(shell_cmd)
		system(shell_cmd)
	}
}
if(length(bases)>7){
	print("more tahn 7 sets:Vennerable plot")
	png(paste(args[2],"/",args[3],".png",sep=""))
	plot(Vstem, type = "battle")
	dev.off()
	pdf(paste(args[2],"/",args[3],".pdf",sep=""))
	plot(Vstem, type = "battle")
	dev.off()
}

