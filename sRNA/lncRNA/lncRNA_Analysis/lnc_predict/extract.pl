#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($i,$o,$cutoff);

GetOptions(
    "i:s"=>\$i,
	"o:s"=>\$o,
	"cutoff:s"=>\$cutoff,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($i and $o and $cutoff);

open (IN,$i) or die $!;
open (OUT,">$o") or die $!;
print OUT "#ID\tmRNA_size\tORF_size\tFickett_score\tHexamer_score\tcoding_prob\tcoding_noncoding\n";
my $n=1;
<IN>;
while (<IN>){
	chomp;
	my $coding_prob=(split /\s+/,$_)[5];
	if($coding_prob < $cutoff){
		print OUT "$_\tno\n";
	}
	$n++;
}
close IN;
close OUT;	
	
#############################################################################################################
sub USAGE {#
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: Yao Bei<yaob\@biomarker.com.cn> 
      Date: 

     Usage:
            -i	cpat.txt
			-o	cpat_filter.txt
			-cutoff	cutoff

   Example:
            perl $Script -i cpat.txt -o cpat_filter.txt -cutoff cutoff 

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}

