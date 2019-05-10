#!/usr/bin/perl -w
# 
# Copyright (c) BIO_MARK 2014
# Writer:         Dengdj <dengdj@biomarker.com.cn>
# Program Date:   2014
# Modifier:       Dengdj <dengdj@biomarker.com.cn>
# Last Modified:  2014.
my $ver="1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

######################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作

my %opts;
GetOptions(\%opts,"i=s","o=s","m=s","x=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{i}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver

	Usage:

		-i           dis infile                                must be given

		-o           out png file                              must be given

		-m           minimum indel length, default -10         optional

		-x           maximum indel length, default 10          optional

		-h           Help document

	Usage End.

	exit;
}

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
################

# get parameters
my $infile  = $opts{i} ;
my $outfile = $opts{o} ;
my $min_length = defined $opts{m} ? $opts{m} : -10 ;
my $max_length = defined $opts{x} ? $opts{x} : 10 ;

my $disfile = &refine_file($infile);

my $cmd = "$Bin/bar_split.r --infile $disfile --outfile $outfile --main.col 1 --sub.col 2 --value.col 3 --main.lab \"Length of Indel\" --sub.lab \"Region\" --value.lab \"Number of Indel\" --title.lab \"Indel length distribution\" --legend.xpos 0.85 --legend.ypos 0.85" ;
&run_or_die($cmd);

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

###############Subs
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;
	$cur_dir =~ s/\n$//;
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;
		$dir =~ s/\n$// ;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;
		$return =~ s/\n$// ;
	}
	else
	{
		warn "Warning just for file and dir\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}

# &show_log("txt")
sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}

#&run_or_die($cmd);
sub run_or_die()
{
	my ($cmd) = @_ ;
	&show_log($cmd);
	my $flag = system($cmd) ;
	if ($flag != 0){
		&show_log("Error: command fail: $cmd");
		exit(1);
	}
	&show_log("done.");
	return ;
}

## qsub

#my $disfile = &refine_file($infile);
sub refine_file()
{
	my ($infile) = @_ ;
	(my $disfile = $infile) =~ s/$/.cut/ ;
	open (IN, $infile) || die "Can't open $infile, $!\n" ;
	open (OUT, ">$disfile") || die "Can't creat $disfile, $!\n" ;
	while (<IN>){
		chomp ;
		next if (m/^\s*$/) ;
		if (m/^-?\d+/){
			my $length = (split)[0] ;
			if ($length >= $min_length && $length <= $max_length){
				print OUT $_,"\n" ;
			}
		}
		else{
			print OUT $_,"\n" ;
		}
	}
	close(IN);
	close(OUT);

	return ($disfile);
}


