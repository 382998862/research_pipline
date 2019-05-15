#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

my ($i,$o);
GetOptions(
				"h|?" =>\&USAGE,
				"i:s"=>\$i,
				"o:s"=>\$o,
				) or &USAGE;
&USAGE unless ($i and $o);

$i=ABSOLUTE_DIR($i);
open(IN,"$i")||die $!;
open(OU,">$o")||die $!;
my %value;
while(<IN>){
	chomp;
	my @line=split/\s+/,$_;
	my $key="$line[0]\t$line[3]\t$line[4]\t$line[2]\t0\t$line[6]";
	$value{$key}=1;
}
foreach(keys %value){
	print OU "$_\n";
}
close OU;
close IN;



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
ProgramName:
produce intron gff file based on  gff3 file
Contact:	Liuxiaoshuang <liuxs\@biomarker.com.cn> 
Program Date:   2015/8/3
Usage:
  Options:
  -i <file>   input file
  
  
  -o  <file>   output file
  
  -h         Help

USAGE
	print $usage;
	exit;
}