library(getopt)
spec=matrix(c(
	'help','h',0,"logical",
	'infile','i',1,"character",
	'outfile','o',1,"character",
	'totalnum','t',1,"integer"
	),byrow=TRUE,ncol=4
)
opt=getopt(spec)
print_usage<-function(spec=NULL){
	cat(getopt(spec,usage=TRUE));
	cat("Usage example: \n")
	cat("
		Description: adjust pvalue 
		Contact: Yanhua Wen <wenyh@biomarker>	
		Usage example:
			Rscript padjust_fdr.r -i pvalue.txt -o adjust_pvalue.txt
		Options:
		--help		NULL		get this help
		--infile	character	the input file contained five column divided by tab
						cerna1	cerna2	cerna_num1	cerna_num2	over_num
		--outfile	character	output file,forced
		--totalnum	character	total background num
		\n")

	q(status=1);
}
if(is.null(opt$infile)) print_usage(spec)
if(is.null(opt$outfile)) print_usage(spec)
if(is.null(opt$totalnum)) print_usage(spec)


read.csv(opt$infile,sep="\t",header=F)->data
data<-as.matrix(data[,3:5])
pvalue<-vector(length=0)
for(i in 1:dim(data)[1]){
#	print(i)
	pvalue[i]<-1-phyper(data[i,3]-1,data[i,1],opt$totalnum-data[i,1],data[i,2])
}

write.table(pvalue,opt$outfile,row.names=F,col.names=F,quote=F,sep="\t")

