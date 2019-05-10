


#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript

#******************************************************************part 1 usage ***************************************

usage=function(){
  cat("
      Descript: 
      Email:    wangmin@biomarker.com.cn
      Options:
  		  bam:     character       input bam file.  [forced]
		  gtf:     character       the merge gtf file [forced]
		  key:     character       the key word for output [forced]
		  od:      character       the  output dir [forced]
     lib:      character       library-type pass to cuffquant --library-type  [default: fr-unstranded]
     run:      T or F          if F, just create sh file, do not qsub sh file. [default: T] 
      usage: 
		Rscript fpkm_saturation_v0.1.R  bam=accepted_hits.bam gtf=merged.gtf od=./saturation key=T01
		  \n"
  )
  q(status=1)
}


#***********************************************part 2 some useful function. you can skip them****************************

getarg=function(usage=function(){cat("No argument! \n"); q(status=1)}){
  #---------------------------------------------------------------------------------
  #Description: getarg is primarily intended to be used with “Rscript”. get arguments from Rscript.
  #             return a list data structure containing names of the flags
  #             the style of = and --(-) are supported
  #             {code} supported pass R code.
  #
  #Example:    Rscript arg.R str=A B str2=c-d log1=T num=1 code={1:10}
  #            Rscript arg.R --str A B -str2 c-d --log1 --num 1 --code {1:10}
  #author:     wangmin@biomarker.com.cn
  #----------------------------------------------------------------------------------
  arg=commandArgs(T)
  if(length(arg)==0){
    usage()
    stop()
  }
  #---------------------------------------------------------
  #match the variable and value
  arg <- paste(arg,collapse=" ")
  if(grepl("^-+",arg)){
    m="(^| )-+[^- ]+"
    style="-+| "
  }else if(grepl("=",arg)){
    m="(^| )[^= ]+="
    style="=| "
  }else{
    usage()
  }
  var=regmatches(arg, gregexpr(m, arg))[[1]]
  var=gsub(style,"",var)
  
  value=regmatches(arg, gregexpr(m, arg),invert =T)[[1]][-1]
  value=gsub("^ | $","",value)
  value[value==""]="T"
  
  #check var and value
  if(length(var)!=length(value)){
    usage()
    stop()
  }
  
  #eval character options to R variable
  arg=list()
  for(i in 1:length(var)){
    if(grepl("\\{.+?\\}",value[i])){
      arg[[var[i]]]=eval(parse(text=gsub("^\\{|\\}$","",value[i])))
    }else{
      arg[[var[i]]]=type.convert(value[i],as.is=T)
    }
  }
  arg
}

##############################################################################################
getRfl=function(){
  #----------------------------------------------------------------------------
  #return the file name and path
  prog <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])
  if(!grepl("/",prog)){prog=file.path(getwd(),prog)}
  prog
}

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
timer=function(sec){
  #------------------------
  #convert second to format hour:minute:second
  sec=ceiling(sec)
  h=sec%/%3600
  m=sec%%3600%/%60
  s=sec%%3600%%60
  paste(h,m,s,sep=":")
}

#Time function
Time=function()format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")

################################
start.run=function(cmd=NULL){
  #---------------------------------------------
  #print running start
  #---------------------------------------------
  if(is.null(cmd)){cmd=getcmd()}
  cat(qq(
    "\n{Time()}  Command start running. 
    the command is:
    \t\t{cmd} \n"
  ))
  proc.time()
}
####################################################
end.run=function(start.time=NULL){
  cat(qq(
    "\n{Time()}  Command run over. 
    The time used (hour:minute:second): {timer(sum(proc.time()-start.time,na.rm=T))} \n"
  ))
}

#-----------------------------------------------------------------------------
qq=function(x,m="\\{code\\}",
            enter=F,blank=F,tab=F,head=F,
            envir=parent.frame()
){
  #--------------------------------------------------
  #variable insert into string like perl
  #argument:
  #        x    the sring in which variables are marked with certain rules
  #        m    the variables match pattern. default '\\{code\\}'
  #   enter:    the '\n' sub to blank is or not. default: F
  #   blank:    the multi-blank sub to single space is or not. default:F
  #     tab:    replace blank to tab
  #   head:     no space at start of every lines
  #   envir:    the enviroment where to find those variables.
  #example:
  #       chrom=1:10
  #       qq(chr={chrom})
  #some code refered by:
  #  https://github.com/jokergoo/GetoptLong/blob/master/R/qq.R
  #--------------------------------------------------
  if(length(x) != 1) {
    stop("Now only support text with length of 1.\n")
  }
  
  if(enter){x=gsub("[\n\r]+"," ",x)}
  if(blank){x=gsub(" +"," ",x)}
  if(tab){x=gsub(" +","\t",x)}
  if(head){x=gsub("\n[\t ]+","\n",x)}
  
  if(!is.null(envir)) {
    if(!is.environment(envir)) {e = as.environment(envir)}
  } 
  
  m2=gsub("code", ".+?", m)
  var=regmatches(x,gregexpr(m2,x))[[1]]
  edge = strsplit(m, "code")[[1]]
  var=gsub(paste("^", edge[1], "|", edge[2], "$", sep=""), "", var)
  str=strsplit(x,m2)[[1]]
  res=""
  if(length(var)>0){
    for(i in 1:length(var)){
      res=paste(res,str[i],eval(parse(text=var[i]),envir = envir),sep="")
    }
  }
  i=length(var)
  if(length(str)>i){res=paste(res,str[i+1],sep="")}
  res 
}

###################################################################

qsub=function(file,text=NULL,queue="general.q",cpu=6,vf="6G",reqsub="--reqsub",run=T,wait=T,
              sh="sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh"){
  #-------------------------------------------------------------------------
  #description: run qsub for sh file
  #argument:
  #         file    character    the file name contains sh command.[forced]
  #-------------------------------------------------------------------------
  if(missing(file)){
    stop("\nsh file path must be specificed\n")
  }
  
  if(!is.null(text)){cat(text,file=file,sep="\n")}
  
  cmd=qq("{sh} {file} --queue {queue} --maxproc {cpu} --resource vf={vf}
         --independent {reqsub}",enter=T,blank=T)
  if(run){
    system(cmd,wait=wait)
    if(!wait){cat("\nqsub submit successfully. please qstat see more information.\n")}
  }
  cmd
}



#******************************************************************part 3 the main function***************************

sample_quant=function(bam,gtf=NULL,key,od=NULL,lib="fr-unstranded",run=T){
  
  ##################################################################################
  downsample=" java -XX:ParallelGCThreads=5 -Xmx4g -jar  /share/nas2/genome/biosoft/picard-tools/1.94/DownsampleSam.jar"
  quant="/share/nas2/genome/biosoft/cufflinks/2.2.1/cuffquant"
  norm="/share/nas2/genome/biosoft/cufflinks/2.2.1/cuffnorm"
  Rscript="/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript"
  #################################################################################
  
  # id=regmatches(bam,regexpr("[TL][0-9]{1,2}", bam))
  
  if(is.null(od)){
    od=sub("accepted_hits.bam$","",bam)
    od=sub("/Tophat/","/sample_cuffquant/",od)
  }
  
  if(substr(od,nchar(od),nchar(od))!="/"){od=paste(od,"/",sep="")}
  
  if(is.null(gtf)){
    gtf=sub("Tophat/[^/]+/accepted_hits.bam","Cuffmerge/merged.gtf",bam)
  }
  
  ######################################################################
  smpath=qq("{od}sample/")
  
  if(!file.exists(smpath)){dir.create(smpath,recursive =T)}
  
  rate=seq(0.1,0.9,0.1)
  
  smfl=qq("{smpath}{rate}.bam")
  
  smcd=qq(
    "{downsample}  I={bam} O={smfl} P={rate}"
    )
  
  #-----------------------------------------
  #make cuffquant sh
  #-----------------------------------------
  smfl=c(smfl,bam)
  rate=seq(0.1,1,0.1)
  qout=qq("{od}cuffquant/p{rate}")
  
  qcd=qq("{quant} -p 6 -no-update-check -q --library-type {lib} -o {qout} {gtf} {smfl} ")
  
  
  #------------------------------------------
  #make cuffnorm sh
  #------------------------------------------
  nmout=paste(od,"cuffnorm",sep="")
  lbl=paste(paste("p",seq(0.1,1,0.1),sep="") ,collapse=",")
  innm=paste(qout,"/abundances.cxb",sep="")
  innm=paste(innm,collapse=" ")
  nmcd=qq("{norm} -p 6 -no-update-check -q --library-type {lib} -L {lbl} -o {nmout} {gtf} {innm} ")
  
  
  #--------------------------------------------------------------------------------------
  #plot 
  #--------------------------------------------------------------------------------------
  bin=dirname(getRfl())
  rcd=qq("{bin}/plot_saturation_fpkm.R")
  plotcd=qq("{Rscript} {rcd} infile={nmout}/genes.fpkm_table outfile={od}Saturation_{key}.png")
  
  #--------------------------------------------------------------------
  #output all work sh
  #--------------------------------------------------------------------
  odsh=qq("{od}worksh/")
  if(!file.exists(odsh)){dir.create(odsh,recursive =T)}
  
  shfiles=c(qq("{odsh}sample_1.sh"),qq("{odsh}cuffquant_2.sh"),
            qq("{odsh}cuffnorm_3.sh"),qq("{odsh}plot_4.sh"))
  
  cat(smcd,file=shfiles[1],sep="\n")
  cat(qcd,file=shfiles[2],sep="\n")
  cat(nmcd,file=shfiles[3],sep="\n")
  cat(plotcd,file=shfiles[4],sep="\n")
  
  #----------------------------------------------------------------------
  #run qsub
  #-----------------------------------------------------------------------
  cmd=qsub(shfiles,cpu=c(5,5,1,1),run=F)
  cmd=paste(cmd,collapse =" && ")
  
  if(run){system(cmd,wait=T)}else{cmd}
}

#*************************************************************call the main function********************************************

stm=start.run()
arg=getarg(usage())
do.call(sample_quant,arg)
end.run(stm)



