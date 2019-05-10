#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript


#************************section 1 usage and some init  function*********************************************** 

usage=function(){
  cat("
      Descript: 
      Email:    wangmin@biomarker.com.cn
    
      Options:
            infile:  character       the path.txt file maked by bam_count.R  [forced]
            group:   character       the group information, like T01_T02_vs_T03_T04 .  will auto match count_list.txt [forced]
            od:      character       the output dir [forced]
            FDR:     float           the FDR cutoff to HTMLReport [default: 0.01]
        
     usage: 
          Rscript DEU_DEXSeq.R infile=path.txt group=T01_T02_vs_T03_T04 od=./Result
      \n"
  )
  q(status=1)
}

getarg=function(usage=function(){cat("No argument! \n")}){
  #-------------------------------------------------------------------------
  #Description: getopt is primarily intended to be used with “Rscript”. get arguments from Rscript.
  #             return a list data structure containing names of the flags
  #author:     wangmin@biomarker.com.cn
  #---------------------------------------------------------------------------
  arg <- commandArgs(T)
  #
  arg=sub("\b"," ",arg)
  if(length(arg)==0){
    do.call(usage,list())
    return(NULL)
  }
  arg=read.table(text=arg,sep="=",row.names=1,as.is=T,colClasses ="character")
  arg=data.frame(t(arg),stringsAsFactors=F)
  arg=as.list(arg)
  # ##eval character options to R variable
  for(i in names(arg)){arg[[i]]=type.convert(arg[[i]],as.is=T)}
  arg
}

###############################################################
getcmd=function(){
  #----------------------------------------
  # get running command from Rscript
  #----------------------------------------
  r=commandArgs()[1]
  r=regmatches(r,regexpr("^.+/bin/",r))
  r=paste(r,"Rscript",sep="")
  prog <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])
  if(grepl("/",prog)){
    path=sub("[^/]+$","",prog)
  }else{
    path=getwd()
    prog=paste(path,prog,sep='/')
  }
  cmd=c(r,prog,commandArgs(T))
  paste(cmd,collapse=" ")
}


################################
start.run=function(){
  #---------------------------------------------
  #print running start
  #---------------------------------------------
  time=format(Sys.time(), "%a %b %d %X %Y")
  cat("\n ",
      time, "   the script start running",
      "\n the run command:\n\t\t",
      getcmd(),
      "\n\n"
  )
  time
}

##########################################################################
end.run=function(start.time=NULL){
  time=format(Sys.time(), "%a %b %d %X %Y")
  cat("\n ",
      time,"  the script run over",
      "\n running from ",start.time," to ", time,
      "\n"
  )
}




#************************section 2 main function*********************************************** 

DEU_DEXSeq=function(infile,group,od="./",FDR=0.01,small=F){
  start.time=start.run()
  
  library(DEXSeq)
  #---od exist or not  and create od dir
  if(substr(od,nchar(od),nchar(od))!="/"){od=paste(od,"/",sep="")}
  if(!file.exists(od)){dir.create(od,recursive =T)}
  
  #--------------------------------------------
  fls=read.table(infile,as.is=T,sep="\t")[,1]
  gff=fls[2]
  
  #------------------split group----------------------------
  group1=strsplit(group,"_vs_")[[1]][1]
  group1=strsplit(group1,"_")[[1]]
  
  group2=strsplit(group,"_vs_")[[1]][2]
  group2=strsplit(group2,"_")[[1]]
  group=c(group1,group2)
  
  #--------------read count file and match group------------------------------
  count.fl=read.table(fls[1],as.is=T,sep="\t")[,1]
  
  count.flm=NULL
  for(i in 1:length(group)){
    count.flm[i]=count.fl[grepl(group[i],count.fl)]
  }
 
  #---------------------------make sample info--------------------------------
  cat("start dexseq and the sample info is :\n")
  condition=c(rep("control",length(group1)),rep("case",length(group2)))
  sampleTable = data.frame(row.names =group,condition =condition)
  sampleinfo=data.frame(file=count.flm,sampleTable)
  
  print(sampleinfo)
  cat("\n")
  
  #---------------------start dexseq--------------------------------
  dxd = DEXSeqDataSetFromHTSeq(
    count.flm,
    sampleData=sampleTable,
    flattenedfile=gff
  )
  dxd$condition <- relevel(dxd$condition, "control")
  
  #extract small data to test.  just used in test
  if(small){
    cat("extract small data to test. just used in test \n")
    warning("extract small data to test. just used in test \n")
    dxd = dxd[geneIDs( dxd ) %in% sample(geneIDs( dxd ),100),]
  }
  
  cat("end read data, start DEU\n")
  dxr=DEXSeq(dxd,BPPARAM=MulticoreParam())
  sigdeu <- any(dxr$padj<FDR, na.rm=T)
  if(sigdeu==F){
    unlink(od, recursive = TRUE)
    cat("There are no significant results in the group \n")
    return (0) 
  }
  cat("end DEU, start HTMLreport\n")
  DEXSeqHTML(dxr,FDR=FDR,path=paste(od,"DEXSeqReport",sep=""),BPPARAM=MulticoreParam())
  
  d=as.data.frame(dxr)
  
  #
  write.table(d, file=paste(od,"DEU_Result_All.xls",sep=""),
              row.names = F,sep="\t",quote=F)
  
  #filter result
  d=subset(d,padj<=FDR)
  
  d=d[,c("groupID","featureID","log2fold_case_control","pvalue","padj")]
  colnames(d)=c("geneID","exonID","log2(FC)","pvalue","FDR")
  
  write.table(d, file=paste(od,"DEU_Result_Final.xls",sep=""),
              row.names = F,sep="\t",quote=F)
  
  end.run(start.time)

}

#************************section 3 call main function*********************************************** 

arg=getarg(usage())
#-------------------------------------------------------------
if(!all(c("infile","group","od")%in% names(arg))) {usage()}

do.call(DEU_DEXSeq,arg)
