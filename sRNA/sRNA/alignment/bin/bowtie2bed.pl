#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(basename dirname);
my ($i,$o);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$i,
				"o:s"=>\$o,
				) or &USAGE;
&USAGE unless ($i and $o);

$i=ABSOLUTE_DIR($i);

open(IN,"$i")||die $!;
open(OUT,">$o")||die $!;
while(<IN>){
	chomp;
	next if /^\@/;
	my @line=split/\t/,$_;
	next if $line[1]==4;
	if($line[1]==0){
		my $chr=$line[2];
		my $start=$line[3];
		my $len=length($line[9]);
		my $end=$start+$len;
		my $id=$line[0];
		my $score=$line[4];
		print OUT "$chr\t$start\t$end\t$id\t$score\t+\n";
	}
	elsif($line[1]==16){
		my $chr=$line[2];
		my $start=$line[3];
		my $len=length($line[9]);
		my $end=$start+$len;
		my $id=$line[0];
		my $score=$line[4];
		print OUT "$chr\t$start\t$end\t$id\t$score\t-\n";
	}
}
close IN;
close OUT;

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
translast bowtie result(sam) to bed format
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