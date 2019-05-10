#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
###############################################################
my %opts;
my ($i,$o);
GetOptions(
    "i:s"   =>\$i,
    "o:s"   =>\$o,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($o and $i);
###############################################################################
$i=abs_path($i);
$o=abs_path($o);
##############################################################################
#my %pfa;
open IN,"$i" || die;
open OUT,">$o" || die;
print OUT "Fusion_Gene_ID\tPfam_acc\tPfam_name\tPfam_description\tBit_score\tE-value\n";
while (<IN>) {
	chomp;
	next if (/^$/);
	next if (/^#/);
	my ($id,$acc,$name,$des,$scor,$val)=(split /\s+/,$_,19)[2,1,0,18,5,4];
	print OUT "$id\t$acc\t$name\t$des\t$scor\t$val\n";
}
close IN;
close OUT;
####################################################################################################3
sub USAGE {
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------------------------------
   Program: $Script
   Date: 2016-09-11
   Contact: MA lm
     Usage:
            --i        <file>   input file
            --o        <file>   output file
            --h                 help documents

   Example:
            perl $Script -i pfam.an  -o pfam.o

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
