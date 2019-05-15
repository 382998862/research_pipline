library(getopt)
spec=matrix(c(
	'help','h',0,"logical",
	'infile','i',1,"character",
	'outfile','o',1,"character",
	'method','m',1,"character"
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
			Rscript p_adjust.r -i pvalue.txt -o adjust_pvalue.txt
		Options:
		--help		NULL		get this help
		--infile	character	the input file contained pvalue without header and rownames
		--outfile	character	output file,forced
		--method	character	Optional method,holm,hochberg,hommel,bonferroni,BH,BY,fdr,default fdr
		\n")

	q(status=1);
}
if(is.null(opt$infile)) print_usage(spec)
if(is.null(opt$outfile)) print_usage(spec)

if(!is.null(opt$method)) method<-opt$method
if(is.null(opt$method)) method<-"fdr"

read.csv(opt$infile,sep="\t",header=F)->pvalue
pvalue<-pvalue[,1]
adjust<-p.adjust(pvalue,method=method,n=length(pvalue))
out<-data.frame(pvalue,adjust)
write.table(out,opt$outfile,row.names=F,col.names=F,quote=F,sep="\t")
