#!/share/nas2/genome/biosoft/R/3.3.2/bin/Rscript
library(getopt)
#library('/home/niepy/R/x86_64-unknown-linux-gnu-library/3.1/getopt/R/getopt')
spec=matrix(c(
	'help','h',0,"logical",
	'go','g',1,"character",
	'kegg','k',1,"character",
	'deg','d',1,"character",
	'all','a',1,"character",
	'other','t',1,"character",
	'od','o',1,"character",
	'prefix','p',1,"character",
	'enrichn','n',1,"integer",
	"gseap",'P',1,"double",
	"gseaq",'Q',1,"double",
	"color",'c',1,"character",
	"column",'C',1,"integer",
	"len","l","1","integer"
	),byrow=TRUE,ncol=4
)
opt=getopt(spec)
print_usage<-function(spec=NULL){
	cat(getopt(spec,usage=TRUE));
	cat("Usage example: \n")
	cat("
		Description: DO enrichment analysis for given gene set
		Contact: Yanhua Wen <wenyh@biomarker>	
		Usage example:
			/share/nas2/genome/biosoft/R/3.3.2/bin/Rscript enrich.r --inpath /share/nas1/wenyh/B237/DMR/J01_vs_J02 --outpath ./out  
			
		Options:
		--help		NULL		get this help
		--deg		character	input DEG_final file or Interested gene sets, file should contained header
						This file would be used for GO/KEGG enrichment analysis.
		--all		character	input DEG.all file, this file should containd column log2FC. Gene sets are all genes.
						This file would be used for GO/KEGG GSEA.
		--od		character	output dir,default current path
		--prefix	character	output key,default Genes

		--color		character	One of pvalue p.adjust qvalue, default p.adjust
		--column	interger	for deg file which column is gene id, default 1
		--kegg		character	KEGG info contained three column
						eg:ko_id	pathway_name	gene_id
		--go		character	GO info contained three column
						eg:go_id	GO_term name	gene_id
		--other		character	other interested gene set function, Satisfy the the column format
						eg:Disease_id	Disease name	gene_id

						Note: One of kegg/go/other must exist.
		--enrichn	int		For enrichment analysis, the first enrichn term would be plotted, default 10
						
		--len		int		The length of GO/KEGG term would be Truncated. default 50
		--gseap		double		GSEA cutoff pvalue, default 1
		--gseaq		double		GSEA cutoff qvalue, default 1
						only terms satisfy the cutoff would be plot.
						if both gseap and gseaq not exist, default gseap=0.001 gseaq=0.05

		\n")

	q(status=1);
}

if(is.null(opt$deg) && is.null(opt$all))	print_usage(spec)
if(is.null(opt$od))	opt$od="./"
if(is.null(opt$prefix)) opt$prefix="Genes"
if(is.null(opt$color))	opt$color="p.adjust"
if(is.null(opt$column))	opt$column=1
if(is.null(opt$len))	opt$len=50
dir.create(opt$od)
#############GO/KEGG/OTHER
if(is.null(opt$go) && is.null(opt$kegg) && is.null(opt$other)) print_usage(spec)
if(!is.null(opt$go)){
	read.csv(opt$go,sep="\t",header=F)->GO
}

if(!is.null(opt$kegg)){
	read.csv(opt$kegg,sep="\t",header=F)->KEGG
}

if(!is.null(opt$other)){
	read.csv(opt$other,sep="\t",header=F)->OTHER
}

#############enrich parameter
if(is.null(opt$gseap) && is.null(opt$gseaq)){
	opt$gseap=0.001
	opt$gseaq=0.05
}
if(is.null(opt$gseap))	{opt$gseap=1}
if(is.null(opt$gseaq))	{opt$gseaq=1}


#############begin enrich
library(clusterProfiler)

#############defined sub function
myenrich<-function(genes,type,info){

	TERM2GENE<-data.frame(info[,1],info[,3])
	TERM2NAME<-data.frame(info[,1],info[,2])
	TERM2NAME<-unique(TERM2NAME[order(TERM2NAME[,1]),])

	enrich<-enricher(genes,TERM2GENE=TERM2GENE,TERM2NAME=TERM2NAME,pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod = "fdr")
	if(is.null(enrich)){
		print("null")
		return(1)
	}
	GeneRatio<-as.numeric(lapply(strsplit(enrich$GeneRatio,split="/"),function(x) as.numeric(x[1])/as.numeric(x[2])))
	BgRatio<-as.numeric(lapply(strsplit(enrich$BgRatio,split="/"),function(x) as.numeric(x[1])/as.numeric(x[2])  ))
	enrich_factor<-GeneRatio/BgRatio
	type<-chartr(" ","_",type)
        out<-data.frame(enrich$ID,enrich$Description,enrich$GeneRatio,enrich$BgRatio,round(enrich_factor,2),enrich$pvalue,enrich$qvalue,enrich$geneID)
        outfile<-paste(opt$od,"/",opt$prefix,"_",type,"_enrich.list",sep="")
        header<-c("ID","Description","GeneRatio","BgRatio","enrich_factor","pvalue","qvalue","geneID")
        write.table(t(header),outfile,row.names=F,col.names=F,quote=F,sep="\t")
        write.table(out,outfile,row.names=F,col.names=F,quote=F,sep="\t",append=T)
	substr<-function(x){
		y<-x
		if(nchar(x)>=opt$len){y=paste(substring(x,1,opt$len),"..",sep="")}
		return(y)
	}
	info[,2]<-sapply(as.character(info[,2]),substr)
	TERM2NAME<-data.frame(info[,1],info[,2])
	TERM2NAME<-unique(TERM2NAME[order(TERM2NAME[,1]),])
	enrich<-enricher(genes,TERM2GENE=TERM2GENE,TERM2NAME=TERM2NAME,pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod = "fdr")
	###plot begin
		bar<-barplot(enrich,showCategory=10,title=type,colorBy=opt$color)
		png(paste(opt$od,"/",opt$prefix,"_",type,"_enrich_barplot.png",sep=""),width = 720)
		print(bar)
		dev.off()

		pdf(paste(opt$od,"/",opt$prefix,"_",type,"_enrich_barplot.pdf",sep=""),width=11)
		print(bar)
		dev.off()

		dot<-dotplot(enrich,x="geneRatio",colorBy=opt$color,showCategory=10,font.size=12,title=type)
		png(paste(opt$od,"/",opt$prefix,"_",type,"_enrich_dotplot.png",sep=""),width = 720)
		print(dot)
		dev.off()

		pdf(paste(opt$od,"/",opt$prefix,"_",type,"_enrich_dotplot.pdf",sep=""),width=11)
       		print(dot)
		dev.off()
		
	###plot done

}

myGSEA<-function(genes,type,info){
	type=chartr(" ","_",type)
	TERM2GENE<-data.frame(info[,1],info[,3])
	TERM2NAME<-data.frame(info[,1],info[,2])
	TERM2NAME<-unique(TERM2NAME[order(TERM2NAME[,1]),])
	gsea<-GSEA(genes,TERM2GENE=TERM2GENE,TERM2NAME=TERM2NAME,pvalueCutoff=1,nPerm=1000)
	
	##plot
	gsea1<-data.frame(gsea)
	which(colnames(gsea1)=="qvalues")->n
	gsea2<-gsea1[order(gsea1[,n]),]
	row<-nrow(gsea2)-1
	which(gsea1$pvalue<opt$gseap)->s1
	which(gsea1$qvalue<opt$gseaq)->s2
	s<-intersect(s1,s2)
	if(row>5){num<-5}else{num<-row}
	if(length(s)>0){if(length(s)<=num){num<-length(s)}}	#the number of (p<0.001 q<0.05) 
#	else{num<-0}
	print(num)	#the number of picture will be plotted
#	if(num>0){	
		for(i in 1:num){
			runningScore<-gseaplot(gsea,gsea2$ID[i],title=gsea2$Description[i],by="runningScore")
			id<-chartr(":","_",gsea2$ID[i])
			preranked<-gseaplot(gsea,gsea2$ID[i],title=gsea2$Description[i],by="preranked")
			png(paste(opt$od,"/GSEA/",opt$prefix,"_",type,"_",id,".runningScore.png",sep=""),height=240)
			print(runningScore)
			dev.off()
			pdf(paste(opt$od,"/GSEA/",opt$prefix,"_",type,"_",id,".runningScore.pdf",sep=""),height=3.5)
			print(runningScore)
			dev.off()
			png(paste(opt$od,"/GSEA/",opt$prefix,"_",type,"_",id,".preranked.png",sep=""),height=240)
			print(preranked)
			dev.off()
			pdf(paste(opt$od,"/GSEA/",opt$prefix,"_",type,"_",id,".preranked.pdf",sep=""),height=3.5)
			print(preranked)
			dev.off()
		}
#	}
	##plot
	out<-data.frame(gsea2$ID,gsea2$Description,gsea2$setSize,gsea2$enrichmentScore,gsea2$NES,gsea2$pvalue,gsea2$qvalues,gsea2$rank,gsea2$core_enrichment)
	outfile<-paste(opt$od,"/",opt$prefix,"_",type,"_GSEA.xls",sep="")
	header<-c("ID","Description","setSize","enrichmentScore","NES","pvalue","qvalue","rank","core_enrichment")
	write.table(t(header),outfile,row.names=F,col.names=F,quote=F,sep="\t")
	write.table(out,outfile,row.names=F,col.names=F,quote=F,sep="\t",append=T)
}

if(!is.null(opt$deg)){
	read.csv(opt$deg,sep="\t",header=T)->deg
	genes<-unique(sort(as.character(deg[,opt$column])))
	if(length(genes)<1){
		print("No deg exists!")
		q()
	}
	if(!is.null(opt$go)){
		types<-c("Biological Process","Molecular Function","Cellular Component")	
		for(i in 1:length(types)){
			which(GO[,4]==types[i])->s
			if(length(s)>0){
				info<-GO[s,1:3]
				myenrich(genes,types[i],info)
			}
		}
	}
	if(!is.null(opt$kegg)){
		myenrich(genes,"KEGG pathway",KEGG)
	}	
	if(!is.null(opt$other)){
		myenrich(genes,"OTHER",OTHER)
	}
}

if(!is.null(opt$all)){
	read.csv(opt$all,sep="\t",header=T)->all
	#ids<-as.character(all[,1])

	which(colnames(all)=="log2FC")->s
	if(length(s)==0){
		print(paste(opt$all," does not have log2FC column",sep=""))
		q()
	}
	dir.create(paste(opt$od,"/GSEA",sep=""))
	df1<-all[which(all$log2FC=="Inf"|all$log2FC=="-Inf"),]
	ids1<-as.character(df1[,1])
	df2<-all[which(all$log2FC!="Inf" & all$log2FC!="-Inf"),]
	ids2<-as.character(df2[,1])
	lists1<-df1$log2FC
	lists2<-df2$log2FC
	a=abs(max(lists2))
	b=abs(min(lists2))
	c=max(a,b)+1
	lists1[lists1=="Inf"] <- c
	lists1[lists1=="-Inf"]<- -c
	ids<-c(ids1,ids2)
	lists<-c(lists1,lists2)
	names(lists)<-ids
	genes<-sort(lists, decreasing = TRUE)
	if(!is.null(opt$go)){
		types<-c("Biological Process","Molecular Function","Cellular Component")	
		for(i in 1:length(types)){
			which(GO[,4]==types[i])->s		
			if(length(s)>0){
				info<-GO[s,1:3]
				myGSEA(genes,types[i],info)	
			}
		}
	}
	if(!is.null(opt$kegg)){
		myGSEA(genes,"KEGG pathway",KEGG)
	}	
	if(!is.null(opt$other)){
		type<-basename(opt$other)
		myGSEA(genes,type,OTHER)
	}
}

