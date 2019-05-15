#!/usr/bin/perl -w
# 
# Copyright (c) BMK 2013
# Writer:         zhangqx <zhangqx@biomarker.com.cn>
# Program Date:   2014/11/7.
# Modifier:       dengdj <dengdj@biomarker.com.cn>
# Last Modified:  2014/11/7.
# ------------------------------------------------------------------
# setting version number
# ------------------------------------------------------------------
my $version="1.0";
use strict;
use Cwd;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);



# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------


my %opts;
GetOptions(\%opts, "in=s", "prefix=s", "outdir=s", "h" ) or USAGE();

# check
USAGE() if(!defined($opts{"in"}) || !defined($opts{"prefix"}) || !defined($opts{"outdir"}) 
	|| defined($opts{"h"}));

# set default value

# ------------------------------------------------------------------
# main pipeline
# ------------------------------------------------------------------
# Time start
my $BEGIN=time();
log_current_time("Start Time");


my @keywords = split(/,/,$opts{"prefix"}) ;
my %keywords;
my $out =  ABSOLUTE_DIR($opts{"outdir"});
foreach  (@keywords) {
	$keywords{$_}= $1;	
}

open (IN, $opts{"in"}) or die $!;

while (<IN>) {
		chomp;
		next if ($_=~/^\s*#/);
		my $line = $_ ;
		my (undef,undef,$feature)=split/\t/,$_;
		
		if (exists $keywords{$feature}) {
			open (A,">>$out/$feature.gff");
			print  A $line,"\n";
		}
}


# Time end
log_current_time("End Time");
run_time($BEGIN);





# ------------------------------------------------------------------
# mkdir dir and check, if failed, then program will die
# ------------------------------------------------------------------
# ------------------------------------------------------------------
# calculate the total run time
# ------------------------------------------------------------------
sub run_time {
	# get parameter
	my ($start_time) = @_;

	# calculate the run time
	my $run_time = time() - $start_time;

	# log
	log_current_time("Total elapsed time: ${run_time}s");
}



# ------------------------------------------------------------------
# log_current_time() and exit(1)
# ------------------------------------------------------------------
sub log_and_exit {
	# log
	log_current_time(shift);
	# exit
	exit(1);
}



# ------------------------------------------------------------------
# print current time with some information
# ------------------------------------------------------------------
sub log_current_time {
	# get parameter
	my ($info) = @_;

	# get current time with string
	my $curr_time = format_datetime(localtime(time()));

	# print info with time
	print "[$curr_time] $info\n";
}



# ------------------------------------------------------------------
# format date time and return
# ------------------------------------------------------------------
sub format_datetime {
	# get parameter
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;


	#$wday = $yday = $isdst = 0;
	my $format_time = sprintf("%4d-%02d-%02d %02d:%02d:%02d", 
		$year+1900, $mon+1, $day, $hour, $min, $sec);

	# return
	return $format_time;
}



# ------------------------------------------------------------------
# get absolute path for file or dir
# ------------------------------------------------------------------
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



# ------------------------------------------------------------------
# print usage information and exit
# ------------------------------------------------------------------
sub USAGE {
	my $usage=<<"USAGE";
Version: 	Version$version
Writer: 	lium <lium\@biomarker.com.cn>
Program Date: 	2013/1/1.
Modifier: 	dengdj <dengdj\@biomarker.com.cn>
Last Modified:  2013/9/16.

Description: 	this program is code template for perl
Usage:
  Example: 	perl gff_prepare.pl -in in.txt -prefix exon,intron -out ./out
  Options:
  Forced parameters:
  -in 		<str> 	the file name of gff file 
  -prefix 	<str> 	the prefix for output
  -out 	<str> 	the dir of output
  -h 		<none> 	Help

USAGE
	print $usage;
	exit(1);
}










