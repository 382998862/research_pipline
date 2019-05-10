#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($cfg, $od, $step, $log, $oneStepOnly);
GetOptions(
				"help|?" =>\&USAGE,
				"cfg:s"  =>\$cfg,
				"od:s"   =>\$od,
				"s:s"    =>\$step,
				"oneStepOnly:s"=>\$oneStepOnly,
				) or &USAGE;
&USAGE unless ($cfg and $od) ;
################################

my $notename = `hostname`; chomp $notename;

$cfg = &ABSOLUTE_DIR($cfg);    &MKDIR($od);
$od  = &ABSOLUTE_DIR($od);     &MKDIR("$od/work_sh");

$step = $step || 1;

#
# log file 
#
my $startTime = GetTime();
my $user = `whoami`;  chomp $user;
my $workDir = `pwd`;  chomp $workDir;
my $task = "perl $Bin/$Script ".join("\t", @Original_ARGV);

open ($log,">", "$od/Tophat_cufflinks.".time().".log") or die $!;

print $log "######################################\n";
print $log "$startTime\n";
print $log "user:	$user\n";
print $log "work directory:	$workDir\n";
print $log "job:	$task\n";
print $log "######################################\n\n";

print $log "config file:  $cfg\n";
print $log "output dir:  $od\n";
print $log "start from step: $step\n\n";

#==================================================================
# bins 
#==================================================================

my $GFFREAD_BIN     = "/share/nas2/genome/biosoft/cufflinks/current/gffread";
my $TOPHAT_BIN      = "/share/nas2/genome/biosoft/tophat/current/tophat2";
## cufflinks suite
my $CUFFLINKS_BIN   = "/share/nas2/genome/biosoft/cufflinks/current/cufflinks";
my $CUFFCOMPARE_BIN = "/share/nas2/genome/biosoft/cufflinks/current/cuffcompare";
my $CUFFMERGE_BIN   = "/share/nas2/genome/biosoft/cufflinks/current/cuffmerge";
my $CUFFDIFF_BIN    = "/share/nas2/genome/biosoft/cufflinks/current/cuffdiff";

#==================================================================
# load config file 
#==================================================================

my %total_read;
my %para;
my %sample;

open (IN,"$cfg") || die "$!";
while (<IN>) {
	chomp;
	s/\r$//;s/^\s+//;s/\s+$//;
	next if (/^\#/ || /^$/);
	
	my @tmp=split /\s+/,$_;
	if ($tmp[0]=~m/Sample/) {
		my $fq1=<IN>;  chomp $fq1;
		my $fq2=<IN>;  chomp $fq2;
		my @fq_1=split /\s+/,$fq1;
		$sample{$tmp[1]}{FQ1}=$fq_1[1];
		if (!-f "$od/totalRead.stat.xls") {
			my $total_line=`less -S $sample{$tmp[1]}{FQ1} | wc -l `;chomp $total_line;
			$total_read{$tmp[1]}=$total_line/2;
		}
		my @fq_2=split /\s+/,$fq2;
		$sample{$tmp[1]}{FQ2}=$fq_2[1];
	}
	$para{$tmp[0]}=$tmp[1];
}
close IN;

if (!-f "$od/totalRead.stat.xls") {
	open OUT,">$od/totalRead.stat.xls" || die;
	foreach my $sam (sort keys %total_read) {
		print OUT "$sam\t$total_read{$sam}\n";
	}
	close OUT;
}
else {
	open IN,"$od/totalRead.stat.xls" || die;
	while (<IN>) {
		chomp;
		next if (/^$/);
		my @tmp=split/\s+/,$_;
		$total_read{$tmp[0]}=$tmp[1];
	}
}

#==================================================================
# pipeline 
#==================================================================

#######################################
#
# step 1: check and build bowtie index
#
######

my $genome;
my $gtf;
my $gff;

if ($step!=1) {
	$genome=basename ($para{Genome_fa});
	$gff = basename ($para{GFF});
	my $gff_name=$gff;
	$gff_name =~ s/\.gff$//ig;
	$gtf = "$gff_name.gtf";
}
if ($step==1) {
	print STDOUT "=== check and build bowtie index ===\n";
	print $log  "=== check and build bowtie index ===\n";

	&MKDIR("$od/Ref_Genome");
	chdir "$od/Ref_Genome";
	$genome=basename ($para{Genome_fa});
	if (!-f "$genome.1.ebwt") {
		system "ln -s $para{Genome_fa} ./";
		system "bowtie-build $genome $genome";
	}

	################################## gff2gtf
	my $gff_name=basename ($para{GFF});
	$gff = $gff_name;
	$gff_name=~s/\.gff$//i;
	system "ln -s $para{GFF} ./";
	system "$GFFREAD_BIN $gff -T -o $gff_name.gtf";
	$gtf = "$gff_name.gtf";
	chdir "../";

	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";

}

#################################### 
#
# step 2: Align the RNA-seq read to genome using Tophat2 + bowtie1
#
#########

if ($step==2) {
	print STDOUT "=== Align the RNA-seq read to genome using Tophat2 + bowtie1  ===\n";
	print STDOUT "Tophat mapping shell file: $od/work_sh/Tophat.sh\n";
	print $log "=== Align the RNA-seq read to genome using Tophat2 + bowtie1  ===\n";
	print $log "Tophat mapping shell file: $od/work_sh/Tophat.sh\n";
	
	#
	# write shell
	#
	open OUT,">$od/work_sh/Tophat.sh" || die;
	&MKDIR("$od/Tophat");
	foreach my $sam (sort keys %sample) {
		&MKDIR("$od/Tophat/$sam");
		print OUT "cd $od/Tophat/$sam && ";
		print OUT "$TOPHAT_BIN -p 6 -G $od/Ref_Genome/$gtf --bowtie1 -N $para{Mis} -r $para{Mate} --library-type $para{Lib} --fusion-search -o ./ $od/Ref_Genome/$genome $sample{$sam}{FQ1} $sample{$sam}{FQ2} && ";
		print OUT "samtools index accepted_hits.bam accepted_hits.bam.bai &&\n";
	}
	close OUT;
	#
	# qsub
	#
	$para{Mem}||="15G";
	$para{Host}||="general.q";
	&Cut_shell_qsub("$od/work_sh/Tophat.sh",6,"$para{Mem}","$para{Host}");
	&Check_qsub_error("$od/work_sh/Tophat.sh");
	

	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}

#####################################
#
# step 3: Statistic bam files
#
########

if ($step==3) {
	print STDOUT "=== Statistic bam files  ===\n";
	print STDOUT "shell file: $od/work_sh/Tophat_bam_stat.sh\n";
	print STDOUT "shell file: $od/work_sh/genome_bam2depth.sh\n";
	print STDOUT "shell file: $od/work_sh/genome_Checkgraph.sh\n";
	print $log "=== Statistic bam files  ===\n";
	print $log "shell file: $od/work_sh/Tophat_bam_stat.sh\n";
	print $log "shell file: $od/work_sh/genome_bam2depth.sh\n";
	print $log "shell file: $od/work_sh/genome_Checkgraph.sh\n";
	
	#
	# write shell 
	#
	open OUT1, ">$od/work_sh/Tophat_bam_stat.sh"   || die;
	open OUT2, ">$od/work_sh/genome_bam2depth.sh"  || die;
	open OUT3, ">$od/work_sh/genome_Checkgraph.sh" || die;
	open OUT4, ">$od/work_sh/plot_ReadDensity.sh"  || die;
	&MKDIR("$od/Map_Stat");
	my @str_type_stat;
	foreach my $sam (sort keys %sample) {
#		&MKDIR("$od/Map_Stat/$sam");
		push @str_type_stat, "$od/Map_Stat/$sam.type.stat";
		print OUT1 "perl $Bin/bin/bam2map_stat.pl -i $sam -bam $od/Tophat/$sam/accepted_hits.bam -totalRead $total_read{$sam} -od $od/Map_Stat\n";
		print OUT2 "samtools depth $od/Tophat/$sam/accepted_hits.bam >$od/Tophat/$sam/$sam.sort.bam.depth &&\n";
		print OUT3 "perl $Bin/bin/get_percent_of_exon_intro_inter_by_gff_v1.4.pl -gff $od/Ref_Genome/$gff -i $od/Tophat/$sam/$sam.sort.bam.depth -od $od/Map_Stat -index $sam &&\n";
#		print OUT4 "perl $Bin/bin/plotReadDensity.pl -i $od/Tophat/$sam/$sam.sort.bam.depth -o $od/Map_Stat -k $sam -png -n 10 &&\n";
        print OUT4 "perl $Bin/bin/plotReadDensity.pl -i $od/Tophat/$sam/$sam.sort.bam.depth -o $od/Map_Stat/ -k $sam -p &&\n";
	}
	close OUT1;
	close OUT2;
	close OUT3;
	close OUT4;

	#
	# qsub
	#
	#&Cut_shell_qsub("$od/work_sh/Tophat_bam_stat.sh",   6, "20G", "general.q");
	`sh $od/work_sh/Tophat_bam_stat.sh`;
	&Cut_shell_qsub("$od/work_sh/genome_bam2depth.sh",  6, "15G", "general.q");
	&Cut_shell_qsub("$od/work_sh/genome_Checkgraph.sh", 6, "5G",  "general.q");
	&Cut_shell_qsub("$od/work_sh/plot_ReadDensity.sh",  6, "5G",  "general.q");
	&Check_qsub_error("$od/work_sh/Tophat_bam_stat.sh");
	
	my $str_type_stat = join " ", @str_type_stat;
	`perl $Bin/bin/map_total_type.pl -i $str_type_stat -o $od/Map_Stat/Total_Type.png`;

	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}

#################################### 
#
# step 4: Cufflinks Assembly Analysis
#
########

if ($step==4) {
	print STDOUT "=== Cufflinks Assembly Analysis ===\n";
	print STDOUT "cufflinks assembly shell file: $od/work_sh/Cufflinks.sh\n";
	print $log "=== Cufflinks Assembly Analysis ===\n";
	print $log "cufflinks assembly shell file: $od/work_sh/Cufflinks.sh\n"; 
	
	#
	# write shell 
	#
	open OUT,">$od/work_sh/Cufflinks.sh" || die;
	&MKDIR ("$od/Cufflinks");
	foreach my $sam (sort keys %sample) {
		&MKDIR("$od/Cufflinks/$sam");
		print OUT "cd $od/Cufflinks/$sam && ";
		print OUT "$CUFFLINKS_BIN -o ./ -p 6 -g $od/Ref_Genome/$gtf --library-type $para{Lib} -u -L $sam $od/Tophat/$sam/accepted_hits.bam && \n";
	}
	close OUT;
	#
	# qsub
	#
	$para{Mem}||="15G";
	$para{Host}||="general.q";
	&Cut_shell_qsub("$od/work_sh/Cufflinks.sh",6,"$para{Mem}","$para{Host}");
	&Check_qsub_error("$od/work_sh/Cufflinks.sh");
	
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}

#################################### 
#
# step 5 & 6: Cuffcompare & Cuffmerge Analysis
#
########

if ($step==5) {
	print STDOUT "=== Cuffcompare: compare gtf of each sample with known gene model ===\n";
	print STDOUT "cuffcompare shell file: $od/work_sh/Cuffcompare.sh\n";
	print $log "=== Cuffcompare: compare gtf of each sample with known gene model ===\n";	
	print $log "cuffcompare shell file: $od/work_sh/Cuffcompare.sh\n";
	
	#
	# write shell 
	#
	open OUT1,">$od/work_sh/Cuffcompare.sh" || die;
	&MKDIR("$od/Cuffcompare");
	open OUT2,">$od/assembly_GTF_list.txt" || die;
	foreach my $sam (sort keys %sample) {
		&MKDIR("$od/Cuffcompare/$sam");
		print OUT1 "cd $od/Cuffcompare/$sam && ";
		print OUT1 "$CUFFCOMPARE_BIN -r $od/Ref_Genome/$gtf $od/Cufflinks/$sam/transcripts.gtf \n";
		print OUT2 "$od/Cufflinks/$sam/transcripts.gtf\n";
	}
	close OUT1;
	close OUT2;
	
	#
	# qsub
	#
	&Cut_shell_qsub("$od/work_sh/Cuffcompare.sh",6,"5G","general.q");
	&Check_qsub_error("$od/work_sh/Cuffcompare.sh");
	
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}

if ($step==6) {
	print STDOUT "=== Cuffmerge: merge gtf files of all samples ===\n";
	print STDOUT "cuffmerge shell file: $od/work_sh/Cuffmerge.sh\n";
	print $log "=== Cuffmerge: merge gtf files of all samples ===\n";
	print $log "cuffmerge shell file: $od/work_sh/Cuffmerge.sh\n";

	open OUT3,">$od/work_sh/Cuffmerge.sh" || die;
	&MKDIR("$od/Cuffmerge");
	#print OUT3 "$CUFFMERGE_BIN -o $od/Cuffmerge -g $od/Ref_Genome/$gtf -s $od/Ref_Genome/$genome -p 6 $od/assembly_GTF_list.txt\n";
	print OUT3 "$CUFFMERGE_BIN -o $od/Cuffmerge -g $od/Ref_Genome/$gtf  -p 6 $od/assembly_GTF_list.txt\n";
	close OUT3;

	#system "$CUFFMERGE_BIN -o $od/Cuffmerge -g $od/Ref_Genome/$gtf -s $od/Ref_Genome/$genome -p 6 $od/assembly_GTF_list.txt\n";
	system "$CUFFMERGE_BIN -o $od/Cuffmerge -g $od/Ref_Genome/$gtf  -p 6 $od/assembly_GTF_list.txt\n";

	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}

##################################### 
#
# step 7: Cuffdiff Analysis
#
#########

if ($step==7) {
	print STDOUT "=== Cuffdiff: different expression and splicing analysis(only use gene FPKM result) ===\n";
	print STDOUT "Cuffdiff shell file: $od/work_sh/Cuffdiff.sh\n";
	print $log "=== Cuffdiff: different expression and splicing analysis(only use gene FPKM result) ===\n";
	print $log "Cuffdiff shell file: $od/work_sh/Cuffdiff.sh\n";
	
	#
	# write shell
	#
	open OUT1,">$od/work_sh/Cuffdiff.sh" || die;
	open OUT2,">$od/work_sh/bam2sam.sh" || die;
	&MKDIR("$od/Cuffdiff");
	print OUT1 "$CUFFDIFF_BIN --library-type $para{Lib} -o $od/Cuffdiff -p 5 $od/Cuffmerge/merged.gtf ";
	print OUT1 "-L ";
	my $lable;
	my $sam_str;
	foreach my $sam (sort keys %sample) {
		print OUT2 "samtools view -h $od/Tophat/$sam/accepted_hits.bam -o $od/Tophat/$sam/accepted_hits.sam && ";
		print OUT2 "grep -v \"XF:\" $od/Tophat/$sam/accepted_hits.sam > $od/Tophat/$sam/accepted_hits.sam.for_AS \n";
		$lable.="$sam".",";
		$sam_str.="$od/Tophat/$sam/accepted_hits.sam ";
	}
	$lable=~s/,$//;
	$sam_str=~s/\s+$//;
	print OUT1 "$lable $sam_str\n";
	close OUT1;
	close OUT2;
	
	#
	# qsub 
	#
	&Cut_shell_qsub("$od/work_sh/bam2sam.sh",6,"15G","general.q");
	&Check_qsub_error("$od/work_sh/bam2sam.sh");
	&Cut_shell_qsub("$od/work_sh/Cuffdiff.sh",6,"30G","great.q");
	&Check_qsub_error("$od/work_sh/Cuffdiff.sh");
	
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}

###################################### Extract Gene Expression For DEG Analysis
#if ($step==8) {
#	open OUT,">$od/work_sh/GeneExpression.sh" || die;
#	&MKDIR("$od/geneExpression");
#	print OUT "perl $Bin/bin/FPKM_extract.pl -list $od/NewGene/$para{SPE}.Cuffmerge.ID.list -id $od/Cuffdiff -od $od/geneExpression -i $para{SPE}\n";
#	close OUT;
#
#	system "perl $Bin/bin/FPKM_extract.pl -list $od/NewGene/$para{SPE}.Cuffmerge.ID.list -id $od/Cuffdiff -od $od/geneExpression -i $para{SPE}";
#	$step=9;
#}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
print $log  "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
close($log);
#==================================================================
# subs 
#==================================================================
sub LOAD_PARA {
	my $para_file= shift;
	my $para= shift;

	my $error_status = 0;
	open IN,$para_file || die "fail open: $para_file";
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r$//;
		next if (/^$/ or /^\#/) ;
		my ($para_key,$para_value) = split(/\s+/,$_);
		$para->{$para_key} = $para_value;
		if (!defined $para_value) {
			warn "Non-exist: $para_key\n";
			$error_status = 1;
		}
	}
	close IN;
	die "\nExit due to error Undefined parameter\n" if ($error_status) ;
}

sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	my $line = `less -S $shell |wc -l `;
	if ($line<=1000) {
		if ($notename=~/cluster/) {
			system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
		}
		else
		{
			system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
		}
	}
	if ($line>1000) {
		my @div=glob "$shell.div*";
		foreach (@div) {
			if (-e $_) {
				system "rm $_";
			}
		}
		@div=();
		my $div_index=1;
		my $line_num=1;
		open IN,"$shell" || die;
		while (<IN>) {
			chomp;
			open OUT,">>$shell.div.$div_index" || die;
			if ($line_num<1000) {
				print OUT "$_\n";
				$line_num++;
			}
			else {
				print OUT "$_\n";
				$div_index++;
				$line_num=1;
				close OUT;
			}
		}
		if ($line_num!=0) {
			close OUT;
		}
		@div=glob "$shell.div*";
		foreach my $div_file (@div) {
			if ($notename=~/cluster/) {
				system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
			}
			else
			{
				system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
			}
		}
	}
}

sub Check_qsub_error {#
	# Check The qsub process if error happend 
	my $sh=shift;
	my @Check_file=glob "$sh*.qsub/*.Check";
	my @sh_file=glob "$sh*.qsub/*.sh";

	if ($#sh_file!=$#Check_file) {
		print "Their Some Error Happend in $sh qsub, Please Check..\n";
		die;
	}
	else {
		print "$sh qsub is Done!\n";
	}
}

sub GetTMR {#
	#Get Total Mapped Reads from file line which contain string "Mapped Reads\s"
	my $fStat=shift;
	open (IN,"<",$fStat) or die $!;
	while (<IN>) {
		if (/^Mapped Reads\s(\d+)/) {
			close (IN) ;
			return $1;
		}
	}
	close (IN) ;
	die "Error Reads Stat file.\n";
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;chomp $return;
	}
	else
	{
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: Tophat&Cufflinks_Analysis Procedure
Version: $version
Contact: Meng Fei <mengf\@biomarker.com.cn>

Description:
	This program is a Procedure deal with Multiple Samples RNA_Analysis
	Tophat+Cufflinks Combination£ºDesigned for RNA Analysis with a Reference Genome

	The program will calculate the Junctions, Transcripts(RABT) Assembly, FPKM of genes & isoforms

Usage:
	-cfg              cfg file              must be given
	-od               output dir            must be given
	-s                step of the program   option,default 1;
                            1  Check the index of Genome & build index for alignment
                            2  run Tophat analysis
                            3  stat Tophat bam out with Genome
                            4  run Cufflinks assembly Transcripts & FPKM
                            5  run Cuffcompare analysis Alt splice
                            6  run Cuffmerge to merge Transcripts
                            7  run Cuffdiff
	-oneStepOnly      perform one of the pipeline steps above mentioned
USAGE
	print $usage;
	exit;
}
