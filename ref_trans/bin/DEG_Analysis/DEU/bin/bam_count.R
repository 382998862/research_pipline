#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript


#************************section 1 usage*********************************************** 

usage=function(){
  cat("
      Descript: count the number of reads for exon by DEXseq. see more DEXseq package manual
      Email:    wangmin@biomarker.com.cn
      Options:
            infile:      character       the file  contains sam or bam files, sam or bam file must sorted by position  [forced]
               gtf:      character       the gtf file [forced]
                od:      character       the output dir [forced]
              type:      character       the type of file: bam or sam. [default: sam]
	     strand:	 character	 the libtype :frstrand [default: yes]
             python:    character	 the python dir  :dir [default: /share/nas2/genome/biosoft/Python/2.7.8/bin/python]
             qsub:       character	 the python dir:dir [default: /share/nas2/genome/biosoft/Python/2.7.8/bin/python]
             queue:      character       the qsub queue [default: general.q]
               cpu:       integer        the qsub cpu    [default: auto select the number of sample ]
                vf:       character      the qsub memory      [default: 5G]
      value:
            the path.txt in dir od  is the input of DEU_DEXSeq.R 
      usage:
          Rscript bam_count.R  infile=bam.list  gtf=genom.gtf od=./
      \n"
  )
  q(status=1)
}

#************************section 2 main function*********************************************** 

bam_count=function(
                   infile,gtf,od,type="sam",strand="yes",
                   queue,qsub,python,cpu=NULL,vf="15G"
                   ){
  library(DEXSeq)
  #python="/share/nas2/genome/biosoft/Python/2.7.8/bin/python"
  if(!file.exists(od)){dir.create(od,recursive =T)}
  if(substr(od,nchar(od),nchar(od))!="/"){od=paste(od,"/",sep="")}
  setwd(od)
  
  #--------------get python scripts in library DEXSeq---------------------------------------------------------------
  an.py=list.files(system.file( "python_scripts", package="DEXSeq" ),pattern="annotation",full.names=T)
  count.py=list.files(system.file( "python_scripts", package="DEXSeq" ),pattern="count",full.names=T)
  
  #--------annotation by python---------------------------------
  cat("start annotation by python\n")
  cat("gtf file: " ,gtf,"\n")
  gff=regmatches(gtf,regexpr("[^/]+$", gtf))
  gff=sub("gtf$","gff",gff)
  gff=paste(od,gff,sep="")
  an.cd=paste(python,an.py,gtf,gff," -r no")
  print(an.cd)
  system(an.cd)
  
  cat("end annotation \n gff file: ",gff,"\n")
  
  #-----count by count.py-----------------------------------
  
  samfile=read.table(infile,as.is=T)[,1]
  
  cat("start count by python.
      sam files: ",samfile,sep="\n")
  
  #ID=regmatches(samfile,regexpr("[A-Za-z0-9]+", samfile))
  ID_tem=strsplit(samfile,"/")
  ID_len=length(ID_tem[[1]])
  ID_tem2=t(as.data.frame(ID_tem))
  ID=as.vector(ID_tem2[,ID_len-1])
  count.fl=paste(od,ID,".count",sep="")
  ct.cd=paste(python,count.py,gff,samfile, count.fl,"-p yes -r pos -s ",strand, " -f",type)
  shfile=paste(od,"count.sh",sep="")
  cat(ct.cd,file=shfile,sep="\n")
  
  #--------------make qsub command and run----------------------------------------------------------------
  if(is.null(cpu)){cpu=length(ID)}
  shcm=paste("sh ", qsub ," ", shfile,
             " --queue ",queue ," --maxproc ",cpu," --resource vf=",vf ,
             " --independent",sep="")
  system(shcm)
  #-------------------------------------------------------------------------------------------------------
  cntfl=paste(od,"count_list.txt",sep="")
  cat(count.fl,file=cntfl,sep="\n")
  
  cat(cntfl,gff,file=paste(od,"path.txt",sep=""),sep="\n")
  cat("python count bam file finished\n")
}

#************************section 3 call main function*********************************************** 

arg <- commandArgs(T)
#
if(length(arg)==0) usage()

# ##eval character options to R variable
arg=read.table(text=arg,sep="=",row.names=1,as.is=T)
arg=data.frame(t(arg),stringsAsFactors=F)
arg=as.list(arg)
for(i in names(arg)){arg[[i]]=type.convert(arg[[i]],as.is=T)}

#-------------------------------------------------------------
if(!all(c("infile","gtf","od")%in% names(arg))) {usage()}

do.call(bam_count,arg)



