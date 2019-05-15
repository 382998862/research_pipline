library(getopt)
spec=matrix(c(
	'help','h',0,"logical",
	'network','n',1,"character",
	'outpath','o',1,"character",
	'key','k',1,"character",
	'dir','d',1,"character",
	'from','f',1,"integer",
	'to','t',1,"integer"
	)
	,byrow=TRUE,ncol=4
)
opt=getopt(spec)
print_usage<-function(spec=NULL){
	cat(getopt(spec,usage=TRUE));
	cat("Usage example: \n")
	cat("
		Description: Random walk for network
		Contact: Yanhua Wen <wenyh@biomarker>	
		Usage example:

		Options:
		--help		NULL		get this help
		--network	character	input network
		--outpath	character	output path
		--key		character	keyword
		--from		int		default, 1
		--to		int		default, 2
		--dir		character	whether Directed network
						TRUE/FALSE, default FALSE
		\n")

	q(status=1);
}
if(is.null(opt$network))	print_usage(spec)
if(is.null(opt$outpath))	print_usage(spec)
if(is.null(opt$from))		opt$from=1
if(is.null(opt$to))		opt$to=2
if(is.null(opt$dir))		opt$dir="FALSE"
if(is.null(opt$key))            opt$key="Net"

library(igraph)
read.csv(opt$network,sep="\t",header=T)->data
edge<-data.frame(data[,opt$from],data[,opt$to])
g <- graph.data.frame(edge, directed = opt$dir)

rank<-page.rank(g,algo=c("prpack","arpack","power"),vids=V(g),directed=opt$dir,damping=0.85,personalized=NULL,weights=NULL,options=NULL)
score<-sort(rank$vector,decreasing = TRUE)
out<-paste(opt$outpath,"/",opt$key,"_score.txt",sep="")
write.table(score,out,quote=F,sep="\t",row.names=T,col.names=F)

