#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
#$Script=abs_path($Script);
#---------------------------------------
# GetOptions
#---------------------------------------
my ($abbre,$ilist,$omature,$opre,$olist);
GetOptions(
	"abbre:s"	=>\$abbre,		## the abbreviation in miRbase
	"ilist:s"	=>\$ilist,		## All_miRNA.list
	"omature:s"	=>\$omature,
	"opre:s"	=>\$opre,
	"olist:s"	=>\$olist,
	"help|?"	=>\&USAGE,
) or &USAGE;
&USAGE unless ($abbre and $ilist and $omature and $opre);

open(IN,$ilist) or die $!;
#my $head=<IN>;
open (MATURE,">$omature") or die $!;
my %MAture;
my %Pairs;
open (PRE,">$opre") or die $!;
open (OUT,">$olist") or die $!;
print OUT "#Pre_ID\tMature_ID\n";
while(<IN>){
	chomp;
	next if (/^$/||/^#/||/^Accession/);
	my @lines=(split /\t/,$_);
	my @ids=(split /\-/,$lines[1]);
	if($ids[0] eq $abbre){
		$lines[3]=~tr/atcguU/ATCGTT/;
		$lines[3]=~tr/VDBHWSKMYRvdbhwskmyr/AATAACTACAaataactaca/;
		print PRE ">$lines[1]\n$lines[3]\n";
		$lines[6]=~tr/atcguU/ATCGTT/;
		$lines[6]=~tr/VDBHWSKMYRvdbhwskmyr/AATAACTACAaataactaca/;
		$MAture{$lines[5]}=$lines[6];
		#print MATURE ">$lines[5]\n$lines[6]\n";
		print OUT "$lines[1]\t$lines[5]";
		if(defined $lines[8]){
			$lines[9]=~tr/atcguU/ATCGTT/;
			$lines[9]=~tr/VDBHWSKMYRvdbhwskmyr/AATAACTACAaataactaca/;
			$MAture{$lines[8]}=$lines[9];
			#print MATURE ">$lines[8]\n$lines[9]\n";
			print OUT ",","$lines[8]";
		}
		print OUT "\n";
	}
}
close IN;
close OUT;
close PRE;

foreach my $id (keys %MAture){
	print MATURE ">$id\n$MAture{$id}\n";
}
close MATURE;



##########sub
sub USAGE{
	my $usage=<<"USAGE";
ProgramName: $Script
ProgramData: 2016.08.10
Usage:
	-abbre      the abbreviation of organism in miRbase, such as 'hsa'
	-ilist      All_miRNA.list
	-omature    the output mature fa file of this organism
	-opre       the output pre fa file of this organism
	-olist      the output miRNA list of this organism
	-h          help document

	Example:perl $Script -abbre hsa -ilist All_miRNA.list -omature hsa.mature.fa -opre hairpin_hsa.fa -olist hsa.miRNA.list
		
USAGE
	print $usage;
	exit;
}