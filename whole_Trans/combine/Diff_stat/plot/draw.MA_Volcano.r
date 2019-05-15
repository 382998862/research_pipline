library(getopt)
spec=matrix(c(
	'help','h',0,"logical",
	'input','i',1,"character",
	'outpath','o',1,"character",
	'key','k',1,"character",
	'type','t',1,"character"
	),byrow=TRUE,ncol=4
)
opt=getopt(spec)
print_usage<-function(spec=NULL){
	cat(getopt(spec,usage=TRUE));
	cat("Usage example: \n")
	cat("
		Contact: Yanhua Wen <wenyh@biomarker>	
		Options:
		--help		NULL		get this help
		--input		character	the input exp path,must be absolute path
		--outpath	character	the input methy path,must be absolute path
		--key		character	out put path
		--type		character	FDR/PValue
		\n")

	q(status=1);
}
if(is.null(opt$input)) print_usage(spec)
if(is.null(opt$type))	opt$type="PValue"
opt$FC=2
opt$p=0.01
library(ggplot2)
data<- read.csv(opt$input,header=T,sep="\t")

Volpng<-paste(opt$outpath,"/",opt$key,".Volcano.png",sep="")
Volpdf<-paste(opt$outpath,"/",opt$key,".Volcano.pdf",sep="")
MApng<-paste(opt$outpath,"/",opt$key,".MA.png",sep="")
MApdf<-paste(opt$outpath,"/",opt$key,".MA.pdf",sep="")
which(data$log2FC != "Inf" & data$log2FC !="-Inf")->N

fc<-data$log2FC[N]
m<-max(abs(fc))

df<-data.frame(log2FC=data$log2FC,pvalue=-log(data[,5],10),lab=factor(data$Sig,levels=c("mRNA","lncRNA","circRNA","miRNA","unchange")))
p=ggplot(df,mapping=aes(x=log2FC, y=pvalue,color=lab))
p=p+geom_point(size=1)
p=p+xlim(-m,m)
p=p+xlab("log2(FC)")+ylab(paste("-log10(",opt$type,")",sep=""))+ggtitle("Volcano plot")+theme(plot.title=element_text(face="bold",size=14))

p=p+geom_hline(yintercept=-log10(opt$p),linetype="longdash",size=0.2,colour="blue")+geom_vline(colour="blue",size=0.2,xintercept=c(log2(opt$FC),log2(1/opt$FC)),linetype="longdash")
p=p+theme_classic()
p=p+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
p=p+scale_color_manual(name="Significant",values=c("mRNA"="red","lncRNA"="blue","circRNA"="green","miRNA"="orange","unchange"="black"))
png(Volpng,width=3000,height=3000,res=300)
print(p)
dev.off()





p=ggplot(data,mapping=aes(x=log2(FPKM), y=log2FC,color=Sig))
p=p+geom_point(size=1)
p=p+ylim(-m,m)
p=p+xlab("log2(FPKM)")+ylab("log2(FC)")+ggtitle("MA plot")+theme(plot.title=element_text(face="bold",size=14))
p=p+theme_classic()
p<- p+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
p=p+scale_color_manual(name="Sig",values=c("mRNA"="red","lncRNA"="blue","circRNA"="green","miRNA"="orange","unchange"="black"))
png(MApng,width=3000,height=3000,res=300)
print(p)
dev.off()


