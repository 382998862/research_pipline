args<-commandArgs(TRUE)
usage <- function() {
  print("-------------------------------------------------------------------------------")
  print("Usage: Rscript TFBS_seqlogo.R TFBSres_allseqs geneID outpath")
  print("1) TFBSres_allseqs: the sequences of TFBS predict result")
  print("2) geneID: gene name")
  print("3) outpath: output path")
  print("-------------------------------------------------------------------------------")
  q(status=1);
}

if(length(args)!=3){
  usage()
}

#######require R-package
library(TFBSTools)
library(Biostrings)

sitesSeqs <- as.matrix(read.table(args[1],header=F,sep="\t"))

if(length(nrow(sitesSeqs))!=0){
	tfseq <- vector()
	name <- vector()
	for (i in 1:nrow(sitesSeqs)) { 
	   tf <- as.character(sitesSeqs[i,2])
	   id <- as.character(sitesSeqs[i,1])
	   tfseq <- c(tfseq,tf)
	   name <- c(name,id)
	}	
	names(tfseq) <- name

	countMatrix <- consensusMatrix(tfseq)
	pfm <- PFMatrix(ID=args[2], name=args[2], profileMatrix=countMatrix)
	icm <- toICM(pfm, pseudocounts=0.8)

	########Plot
	pdffile = paste(args[3],"/",args[2],"_seqlogo.pdf",sep="")
	pdf(pdffile)
	seqLogo(icm)
	dev.off()
	#png(paste(args[3],"/",args[2],"_seqlogo.png",sep=""))
	#seqLogo(icm)
	#dev.off()
}else{
	stop("There is no sequence of TFBS, cannot plot.") 
}
