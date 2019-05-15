library(getopt)
spec=matrix(c(
	'help','h',0,"logical",
	'input','i',1,"character",
	'outpath','o',1,"character",
	'key','k',1,"character"
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
		--input		character	the input file,forced
						Contained three column
						group	RNA_type	Number
						W1	mRNA		20
						W2      lncRNA            20
		--outpath	character	output path,forced
		--key		character	key word,forced

		\n")

	q(status=1);
}
if(is.null(opt$input)) print_usage(spec)
if(is.null(opt$outpath)) print_usage(spec)
if(is.null(opt$key)) print_usage(spec)



library(ggplot2)
library(RColorBrewer)

read.csv(opt$input,header=T,sep="\t")->data
pngfile<-paste(opt$outpath,"/barplot_",opt$key,".png",sep="")
pdffile<-paste(opt$outpath,"/barplot_",opt$key,".pdf",sep="")

len<-length(unique(as.character(data$RNA_type)))
mycol<-brewer.pal(8, "Set1")
col<-mycol[1:len]
p=ggplot(data=data,aes(x=group,y=Number,fill=RNA_type))+geom_bar(position = "fill",stat="identity")+scale_fill_manual(values = col)
p=p+xlab("Difference Group")+ylab("Ratio") + theme(axis.title.y = element_text(size = 25,face = "bold"), axis.title.x = element_text(size = 25,face = "bold"))
p=p+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+scale_y_continuous(expand = c(0,0))

png(pngfile,width=3000,height=3000,res=300)
print(p)
dev.off()
pdf(pdffile)
print(p)
dev.off()









