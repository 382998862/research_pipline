#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.0";

my($i,$o);
GetOptions(
	"help|?"=>\&USAGE,
	"i:s"=>\$i,
	"o:s"=>\$o,
) or &USAGE;
&USAGE unless ($i and $o);

$i=abs_path($i);
$o=abs_path($o);

my %ID;
open(IN,$i) or die $!;
open (OUT,">$o") or die $!;
$/=">";
<IN>;
while(<IN>){
	chomp;
	my ($id,$seq)=(split /\n/,$_,2)[0,1];
	$seq=~s/\s+//g;
	$seq=~tr/a-z/A-Z/;
	$seq=~tr/U/T/;
	if(length $seq>0){
		print OUT ">$id\n$seq\n";	
	}
}
close IN;
close OUT;
$/="\n";




sub USAGE {#
	my $usage=<<"USAGE";
Program: $Script
Version: $version
Usage:
  Options:
	-i	the input file,must be fa format
	-o	the output file,also fa format
	-h	Help

USAGE
	print $usage;
	exit;
}
