#!/usr/bin/perl -w
#Writer           songmm <songmm@biomarker.com.cn>
my $version="1.0.0";
my $BEGIN=time();

use strict;
use warnings;
use experimental 'smartmatch';
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd 'abs_path';
my $programe_dir=basename($0);
my $path=dirname($0);
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($indir,$o);
GetOptions(
	"help|?" =>\&USAGE,
	"indir:s"=>\$indir,
	"o:s"=>\$o
	) or &USAGE;
&USAGE unless ($indir and $o) ;
my %uniq;
my @files = glob("$indir/*.fa");
$/=">";
open OUT ,">$o" or die $!;
foreach my $file(@files){
	open IN,$file or die $!;
	while(<IN>)
	{	
		chomp;
		next if(/^\s*$/);
		my ($head,$seq)=split/\n/,$_,2;
		if(exists $uniq{$head}){next;}
		else
		{
			print OUT ">$_";
			$uniq{$head}=1;
		}
	}
	close IN;
}
close OUT;
sub USAGE
{
    	my $usage=<<"USAGE";
Program: combination all circRNA sequence
Version: $version
Contact: songmm <songmm\@biomarker.com.cn>

Description:

Usage:
  -indir  <dir>  the circRNA fasta sequence dir of all sample 
  -o   <file>   output file  must be given;
USAGE
	print $usage;
	exit;
}
