#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.0";
my %config=%{readconf("$Bin/../config/db_file.cfg")}; 

my @Original_ARGV=@ARGV;
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($cfg1, $cfg2, $od, $step, $oneStepOnly);
GetOptions(
				"help|?" =>\&USAGE,
				"cfg1:s"  =>\$cfg1,
				"cfg2:s"  =>\$cfg2,
				"od:s"   =>\$od,
				"s:s"    =>\$step,
				"oneStepOnly:s"=>\$oneStepOnly,
				) or &USAGE;
&USAGE unless ($cfg1 and $cfg2 and $od) ;
################################

my $notename = `hostname`; chomp $notename;

$cfg1 = &ABSOLUTE_DIR($cfg1);
$cfg2 = &ABSOLUTE_DIR($cfg2);
&MKDIR($od);
$od  = &ABSOLUTE_DIR($od);     &MKDIR("$od/work_sh");

$step = $step || 1;

my $startTime = GetTime();
my $user = `whoami`;  chomp $user;
my $workDir = `pwd`;  chomp $workDir;
my $task = "perl $Bin/$Script ".join("\t", @Original_ARGV);

#==================================================================
# bins 
#==================================================================

my $GFFREAD_BIN     = "$config{cufflink}/gffread";     # 2014-12-17 ~ 
my $TOPHAT_BIN      = "$Bin/basic_analysis/v1.2/tophat_cufflinks/tophat/v2.0.13/tophat2";     # 2015-01-28 ~ 

my $HISAT=$config{HISAT};
my $HISAT_BUILD=$config{HISAT_BUILD};

#==================================================================
# load config file 
#==================================================================

my %para;
my %sample;
open (IN,"cat $cfg2 $cfg1|") || die "$!\n";
open(STAT,">$od/work_sh/read_count.sh")||die $!;
`mkdir $od/Read`	unless(-d "$od/Read");
while (<IN>) {
	chomp;
	s/\r$//;s/^\s+//;s/\s+$//;
	next if (/^\#/ || /^$/);	
	my @tmp=split /\s+/,$_;
	if ($tmp[0]=~m/^Sample/) {
		my $fq1=<IN>;  chomp $fq1;
		my $fq2=<IN>;  chomp $fq2;
		my @fq_1=split /\s+/,$fq1;
		$sample{$tmp[1]}{FQ1}=$fq_1[1];
		my @fq_2=split /\s+/,$fq2;
		$sample{$tmp[1]}{FQ2}=$fq_2[1];
		print STAT "less $fq_1[1]|wc -l >$od/Read/$tmp[1].read.xls\n";
	}
	$para{$tmp[0]}=$tmp[1];
}
close IN;	close STAT;
&qsubOrDie("$od/work_sh/read_count.sh","$para{Queue_type1}",80,"$para{Memory}");
my @files=glob("$od/Read/*.read.xls");
open(OUT,">$od/totalRead.stat.xls")||die $!;
foreach my $f(@files){
	my $sam=(split(/\./,basename($f)))[0];
	my $info=`head -n 1 $f`;
	chomp($info);	$info=(split(/\s+/,$info))[0]/2;
	print OUT "$sam\t$info\n";
}
close(OUT);

#######################################
# step 1: check and build bowtie index
######
my $genome;
my $idx_prefix;
my $gtf;
my $gff;
my $index = $para{'Project_key'};

if ($step!=1) {
    $genome = basename($para{Ref_seq});
    $idx_prefix = basename($para{Ref_seq});
    $idx_prefix =~s/.fa$//;
    $gff = basename($para{Ref_ann});
    $gtf = basename($para{Ref_ann});
    $gtf =~s/\.gff3?$//i;
    $gtf.= ".gtf";
}
if ($step==1) {
	print STDOUT "=== check and build bowtie index ===\n";

	&MKDIR("$od/Ref_Genome");
	chdir "$od/Ref_Genome";
	$genome = basename($para{Ref_seq});
	$idx_prefix = basename($para{Ref_seq});
	$idx_prefix =~s/.fa$//;
	if ( !-f "$idx_prefix.hdrs" ) {
        	`grep '>' $para{Ref_seq} > $idx_prefix.hdrs`;
	}
	my $genome_path = dirname($para{Ref_seq});
	if (!-f "$idx_prefix.1.ht2" and !-f "$idx_prefix.1.ht2l" ) {
		if (-e "$genome_path/$idx_prefix.1.ht2") {
			system "ln -s $genome_path/$genome $genome_path/$idx_prefix.*.ht2 ./";
		} elsif (-e "$genome_path/$idx_prefix.1.ht2l") {   # for large genomes
			system "ln -s $genome_path/$genome $genome_path/$idx_prefix.*.ht2l ./";
		}else {
			system "ln -s $genome_path/$genome ./";
			system "$HISAT_BUILD -p 4 $genome $idx_prefix";
		}
	}
	if (!-f "$idx_prefix.fai"){
	        if (-e "$genome_path/$idx_prefix.fai"){
			system "ln -s $genome_path/$idx_prefix.fai ./";
        	}else{
			system "$config{samtools} faidx $idx_prefix.fa";
        	}
	}

	################################## gff2gtf
	$gff = basename($para{Ref_ann});
	$gtf = basename($para{Ref_ann});
	$gtf =~s/\.gff3?$//i;
	$gtf.= ".gtf";
	system "ln -s $para{Ref_ann} ./" unless (-f $gff);
	system "$GFFREAD_BIN $gff -T -o $gtf";
	system "perl $Bin/basic_analysis/v1.2/bin/check_gtf.pl $gtf $gff $gtf.DEU";
	system "ln -s $para{Ref_seq} ./" unless (-f "$idx_prefix.fa");
	chdir "../";

	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
}

#################################### 
# step 2: Align the RNA-seq read to genome using Hisat2
#########
my $strand;
if ($para{'Lib_type'} eq "fr-firststrand" ){
        $strand="rf";
}elsif ($para{'Lib_type'} eq "fr-secondstrand"){
    $strand="fr";
}elsif ($para{'Lib_type'} eq "fr-unstranded") {
    $strand="ff";
}

if ($step==2) {
	print STDOUT "=== Align the RNA-seq read to genome using hisat + stringTie  ===\n";
	print STDOUT "Mapping shell file: $od/work_sh/hisat.sh\n";
	
	open(OUT,">$od/work_sh/Hisat.sh") || die $!;
	open(SAM,">$od/work_sh/samtools.sh")||die $!;	
	&MKDIR("$od/Hisat");

	foreach my $sam (sort keys %sample) {
		&MKDIR("$od/Hisat/$sam");
		print OUT "cd $od/Hisat/$sam && $HISAT  --dta --$strand  -p 6 -x $od/Ref_Genome/$idx_prefix -1 $sample{$sam}{FQ1} -2 $sample{$sam}{FQ2} -S $sam.HISAT_aln.sam ";
		print OUT " --phred64 "		if(exists $para{Qphred} && $para{Qphred} eq "64");
		print OUT "\n";
		print SAM "cd $od/Hisat/$sam && $config{samtools} view -F 4 -Su $sam.HISAT_aln.sam | $config{samtools} sort - $sam.HISAT_aln.sorted && $config{samtools} index $sam.HISAT_aln.sorted.bam $sam.HISAT_aln.sorted.bam.bai \n";
	}
	close OUT;
	close SAM;
	$para{Memory}||="15G";
	$para{Queue_type1}||="medical.q";

	&qsubOrDie("$od/work_sh/Hisat.sh","$para{Queue_type1}",18,"$para{Memory}");
	&qsubCheck("$od/work_sh/Hisat.sh");
        &qsubOrDie("$od/work_sh/samtools.sh","$para{Queue_type1}",18,"$para{Memory}");
        &qsubCheck("$od/work_sh/samtools.sh");
	
	print STDOUT "\n";
}
`rm $od/Hisat/*/*.HISAT_aln.sam`;
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
Program: Hisat_Analysis Procedure
Version: $version
Contact: Meng Fei <mengf\@biomarker.com.cn>
Version: $version
Modifier: Wang Yajing <wangyj\@biomarker.com.cn>

Description:
	This program is a Procedure deal with RNA_Analysis
	Hisat Designed for RNA Analysis with a Reference Genome

	The program will calculate the Junctions

Usage:
    -cfg1             data config, rawdata & refseq path
    -cfg2             detail config, analysis parameters
	-od               output dir            must be given
	-s                step of the program   option,default 1;
                            1  Check the index of Genome & build index for alignment
                            2  run Hisat analysis
	-oneStepOnly      perform one of the pipeline steps above mentioned
USAGE
	print $usage;
	exit;
}
