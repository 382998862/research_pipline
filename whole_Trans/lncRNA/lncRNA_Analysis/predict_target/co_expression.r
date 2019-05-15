library(getopt)
spec=matrix(c(
	'help','h',0,"logical",
	'rna1','i',1,"character",
	'rna2','r',1,"character",
	'outpath','o',1,"character",
	'cor','c',1,"double",
        'pvalue','p',1,"double",
        'top','t',1,"double",
	'key','k',1,"character",
	'method','m',1,"character"
	),byrow=TRUE,ncol=4
)
opt=getopt(spec)
print_usage<-function(spec=NULL){
	cat(getopt(spec,usage=TRUE));
	cat("Usage example: \n")
	cat("
		Description:
		Contact: Yanhua Wen <wenyh@biomarker>	
		Usage example:
			
		Options:
		--help		NULL		get this help
		--rna1		character	input rna1 file
		--rna2		character	input rna2 file, default rna1
		--outpath	character	output path
		--key		character	output keyword,default relation
		--method	character	cor-relation method(pearson, kendall, spearman),default pearson
		--cor		double		default 0
		--pvalue        double		default 1
						if cor and pvalue not defined, will produce all pairs
		--top           double		default
		\n")

	q(status=1);
}
if(is.null(opt$rna1))		print_usage(spec)
if(is.null(opt$rna2))		opt$rna2=opt$rna1
if(is.null(opt$outpath))	print_usage(spec)

if(is.null(opt$method))		opt$method="pearson"
if(is.null(opt$cor)) 		opt$cor=0
if(is.null(opt$pvalue)) 	opt$pvalue=1
if(is.null(opt$key))		opt$key="relation"

opt$cor=as.numeric(opt$cor)
opt$pvalue=as.numeric(opt$pvalue)
print(opt$cor)
print(opt$pvalue)
read.csv(opt$rna1,sep="\t",header=T)->data1
read.csv(opt$rna2,sep="\t",header=T)->data2

gene1<-as.character(data1[,1])
gene2<-as.character(data2[,1])
data1<-as.matrix(data1[,2:dim(data1)[2]])
data2<-as.matrix(data2[,2:dim(data2)[2]])

outsig<-paste(opt$outpath,"/co_expression.",opt$key,".Sig.xls",sep="")
out<-t(c("RNA1","RNA2","coefficient","pvalue"))
write.table(out,outsig,row.names=F,col.names=F,quote=F,sep="\t")

for(i in 1:dim(data1)[1]){
	for(j in 1:dim(data2)[1]){
		if(gene1[i] == gene2[j])	next			

		cor<-cor.test(as.numeric(as.character(data1[i,])),as.numeric(as.character(data2[j,])),method=opt$method)
		p<-cor$p.value
		co<-cor$estimate
		out<-t(c(gene1[i],gene2[j],co,p))
		if(is.na(p)||is.na(co)){
			next
		}
		if(abs(co)>=opt$cor && p<=opt$pvalue){
			write.table(out,outsig,row.names=F,col.names=F,quote=F,sep="\t",append=T)
		}
	}
}







