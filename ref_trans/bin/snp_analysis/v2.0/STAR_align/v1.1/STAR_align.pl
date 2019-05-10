#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.1.0";
my @Original_ARGV=@ARGV;
#######################################################################################

# ==============================================================
# Get Options
# ==============================================================
my ($cfg,$dOut);

GetOptions(
				"help|?" =>\&USAGE,
				"cfg:s"=>\$cfg,
				"od:s"=>\$dOut,
				) or &USAGE;
&USAGE unless ($cfg);
my $step=0;
#===============================================================
# Default optional value 
#===============================================================
my $_STAR_		=	"$Bin/bin/v2.4.0j/STAR";
my $_reflen_	=	"$Bin/bin/ref_GC_len.pl";
$dOut||=abs_path."/STAR_Result";&MKDIR($dOut);
$dOut=abs_path($dOut);
if (-e "$dOut/work_sh") {
	die("work shell exists!\n");
}
&MKDIR("$dOut/work_sh");
my $notename = `hostname`; chomp $notename;

#===============================================================
# Load Para 
#===============================================================
my %para;
my %sample;
&read_config($cfg);
my $Thread=6;
my $cpu=10;
my $queue=$para{'queue'};
my $vf=$para{'vf'};
$vf||="50g";
$queue||="middle.q";
$para{readlength}=101 if(!exists $para{readlength} or $para{readlength} !~ /\d/g);
my $readlength=$para{readlength};
my $genomeFile=$para{'genome'};

#=================================================================
#  Generate Genome DIR
#=================================================================
my $lenfile="$dOut/tmp/".basename($genomeFile).".len";
`perl $_reflen_ -ref $genomeFile -od $dOut/tmp` if(!-e $lenfile);
my $chrNum=`grep -v \"#\" $lenfile|wc -l `;
my $genomeLength=`tail -1 $lenfile|cut -f2`;
chomp $chrNum;
chomp $genomeLength;
my $SAindex=int(Min(14, log2($genomeLength)/2-2));
my $ChrBin =int(Min(18, log2($genomeLength/$chrNum)));
my $outGenomeDir	=	$dOut."/Genome";
&MKDIR("$outGenomeDir");
&GenerateGenome();
#=================================================================
#  align to Genome
#=================================================================
&MKDIR("$dOut/Alignment");
foreach my $sample (sort keys(%sample)) {
	&alignToGenome($sample,$sample{$sample}{FQ1},$sample{$sample}{FQ2});
}
my @shell=glob("$dOut/work_sh/STAR*.sh");
foreach my $shell (@shell) {
	$shell=~/STAR(\d).sh$/;
	next if($1<$step);
	print "doing $shell\n";
	&Cut_shell_qsub($shell,$cpu,$vf,$queue);
	&Check_qsub_error($shell);
}

#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ==============================================================
# sub function
# ==============================================================
sub GenerateGenome{
	my $starcmd1		=	"$_STAR_ --runMode genomeGenerate --genomeDir $outGenomeDir --genomeFastaFiles $genomeFile  --runThreadN $Thread ";
	   $starcmd1		.=	" --genomeChrBinNbits $ChrBin --genomeSAindexNbases $SAindex";
	   $starcmd1		.=	"\n";
	open (STAR1,">>$dOut/work_sh/STAR1.sh") or die $!;
	print STAR1 $starcmd1;
	close(STAR1);
}

sub alignToGenome{
	my ($sample,$fq1,$fq2) = @_;
	my $runDir		=	"$dOut/Alignment/$sample/";
	&MKDIR($runDir);
	#step1 generate genome
	# step2 mapping genome
	my $sjdbover=$readlength-1;
	my $starcmd2		=	"cd $runDir && $_STAR_ --genomeDir $outGenomeDir --readFilesIn $fq1 $fq2 --runThreadN $Thread --twopass1readsN 10000000000 --sjdbOverhang  $sjdbover ";
	if ($fq1 =~/\.gz$/ and $fq2 =~/\.gz$/) {
		$starcmd2.=" --readFilesCommand zcat ";
	}
	if ($para{'forCuff'} eq "Y") {
		$starcmd2		.=	" --outSAMstrandField intronMotif ";
	}
	$starcmd2.=" --alignSoftClipAtReferenceEnds No ";
	$starcmd2.=" --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM  90000000000";
	$starcmd2.="\n";
	open (STAR2,">>$dOut/work_sh/STAR2.sh") or die $!;
	print STAR2 $starcmd2;
	close(STAR2);
}
sub log2{
	my($a)=@_;
	my $v=log($a)/log(2);
	return $v;
}
sub Min{
	my($a,$b)=@_;
	if ($a<$b) {
		return ($a);
	}else{
		return ($b);
	}
}

sub read_config {#
	my($cfg)=@_;
	open (IN,"$cfg") || die "$!";
	while (<IN>) {
		chomp;
		s/\r$//;
		s/\s+$//;
		next if (/\#/ || /^$/);
		my @tmp=split /\s+/,$_;
		if ($tmp[0]=~m/Sample/) {
			my $fq1=<IN>;chomp $fq1;
			my $fq2=<IN>;chomp $fq2;
			my @fq_1=split /\s+/,$fq1;
			$sample{$tmp[1]}{FQ1}=$fq_1[1];
			my @fq_2=split /\s+/,$fq2;
			$sample{$tmp[1]}{FQ2}=$fq_2[1];
		}
		$para{$tmp[0]}=$tmp[1];
	}
	close IN;
}
sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	mkdir($dir) if(!-d $dir);
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
sub check_file {#
	my ($file) = @_;
	if (-e $file) {
		return 1;
	}else{
		return 0;
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
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Shi Tongwei <shitw\@biomarker.com.cn> 
Discription:
===============================================================================
config file template: $Bin/starconf.cfg
		20150203 v1.1 update: change sam output to bam  
		                      merge 2-pass pipeline in one step
===============================================================================
Usage:
  Options:
  -cfg	<str>	config file,forced
  -od	<str>	Directory where output file produced,optional,default [./STAR_Result]
  -h		Help

USAGE
	print $usage;
	exit;
}
