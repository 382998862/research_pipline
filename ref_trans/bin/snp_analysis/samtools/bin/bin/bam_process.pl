#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my ($od,$bam);
use newPerlBase;
my %config=%{readconf("$Bin/../../../../../config/db_file.cfg")}; 


GetOptions(
			"help|?"=>\&USAGE,
			"od:s"=>\$od,
			"bam:s"=>\$bam,
			)or &USAGE;
&USAGE unless ($bam and $od);
`mkdir $od ` unless (-d $od);

$od=ABSOLUTE_DIR($od);
$bam=ABSOLUTE_DIR($bam);
my $sample=basename(dirname($bam));


my $temp=`$config{samtools} view -h $bam |grep  "[0-9]*N[0-9]*D[0-9]*N" `;  ###GATK　can not recognize
chomp$temp;
if ($temp=~/^$/) {
	my $step=` ln -s $bam $od/$sample.bam`;
}
else {
	$temp=`$config{samtools} view -h $bam |grep -v  "[0-9]*N[0-9]*D[0-9]*N"  >$od/$sample.sam`;
	$temp=`$config{samtools} view -bS $od/$sample.sam > $od/$sample.bam`;
	$temp=`rm $od/$sample.sam `;
}
##################################################

################################################################################
sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:  blast_supply
Version:	$version
Contact:	Wang Yajing <wangyj\@biomarker.com.cn> 
Program Date:   2015/5/27
Usage:
Options:
		-od     output dir
		-bam     bam file
		-h      Help

USAGE
	print $usage;
	exit;
}
