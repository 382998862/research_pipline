#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################

# ==============================================================
# Get Options
# ==============================================================
my (@fIn,$fSNP,$fKey,$dOut,$log);

GetOptions(
				"help|?" =>\&USAGE,
				"i:s{,}"=>\@fIn,
				"s:s"=>\$fSNP,
				"d:s"=>\$dOut,
				) or &USAGE;
&USAGE unless (@fIn and $fSNP);

#===============================================================
# Default optional value 
#===============================================================
$dOut||="./";
mkdir $dOut unless (-d $dOut) ;

#===============================================================
# Global value
#===============================================================

my %stat = ();

#===============================================================
# Get Data
#===============================================================

foreach my $infile (@fIn) {
	open (IN,"$infile") or die $!;
	my ($sam) = basename($infile) =~/(\w+)\.stat\.xls$/;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/) ;
		my ($chro, $snpNum, $transition, $tranversion, $heter) = split;

		$stat{$chro}{$sam}{'nsnp'} = $snpNum;
		$stat{$chro}{$sam}{'transition'} = $transition;
		$stat{$chro}{$sam}{'tranversion'} = $tranversion;
		$stat{$chro}{$sam}{'heter'} = $heter;
	}
	close (IN) ;
}

my @sam = ();

open (IN,"$fSNP") or die $!;
while (<IN>) {
	chomp;
	next if (/^$/) ;

	if (/^\#/) {
		my ($chro, $pos, $geneID, $refBase, @snpInfo) = split;
		
		for (my $i=0; $i<@snpInfo; $i+=4) {
			my ($sam, undef) = split /_/, $snpInfo[$i];
			push @sam, $sam;
		}
		next;
	}

	my ($chro, $pos, $geneID, $refBase, @snpInfo) = split;
	for (my $i=0, my $j=0; $i<@snpInfo; $i+=4,$j++) {
		next if ($snpInfo[$i] eq '.') ;
		if ($geneID ne 'intergenic') {
			$stat{$chro}{$sam[$j]}{'gene_region'}++;
			$stat{'Total'}{$sam[$j]}{'gene_region'}++;
		}else{
			$stat{$chro}{$sam[$j]}{'intergenic'}++;
			$stat{'Total'}{$sam[$j]}{'intergenic'}++;
		}

		die "wrong format of merge snp file" if ($j > $#sam) ; 
	}
}
close (IN) ;

#print Dumper %stat; die;

#===============================================================
# output
#===============================================================

foreach my $sample (@sam) {

	open (OUT,">$dOut/$sample.snp.stat.xls") or die $!;

	print OUT join("\t",
		'#Chro',
		(map {("Snp_number",
		"Snp_number_gene",
		"Snp_number_intergenic",
		"Transition",
		"Transversion",
		"Heterozygosity")} ($sample)), 
	),"\n";

	foreach my $chro (sort {$a cmp $b} keys %stat) {
		next if ($chro eq 'Total') ;
		map {
			$stat{$chro}{$_}{'gene_region'}||=0;
			$stat{$chro}{$_}{'intergenic'}||=0;
			$stat{$chro}{$_}{'nsnp'}||=0;
		    $stat{$chro}{$_}{'transition'}||=0;
		    $stat{$chro}{$_}{'tranversion'}||=0;
		    $stat{$chro}{$_}{'heter'}||=0;
		}($sample);

		print OUT join("\t",
			$chro, 
			(
				map {
					($stat{$chro}{$_}{'nsnp'}, $stat{$chro}{$_}{'gene_region'},$stat{$chro}{$_}{'intergenic'},$stat{$chro}{$_}{'transition'},$stat{$chro}{$_}{'tranversion'},$stat{$chro}{$_}{'heter'});
				} ($sample)
			), 
		),"\n";
	}
	map {
		$stat{'Total'}{$_}{'gene_region'}||=0;
		$stat{'Total'}{$_}{'intergenic'}||=0;
		$stat{'Total'}{$_}{'nsnp'}||=0;
		$stat{'Total'}{$_}{'transition'}||=0;
		$stat{'Total'}{$_}{'tranversion'}||=0;
		$stat{'Total'}{$_}{'heter'}||=0;
	}($sample);
	print OUT join("\t",
		'Total', 
		(
			map {
				($stat{'Total'}{$_}{'nsnp'}, $stat{'Total'}{$_}{'gene_region'},$stat{'Total'}{$_}{'intergenic'},$stat{'Total'}{$_}{'transition'},$stat{'Total'}{$_}{'tranversion'},$stat{'Total'}{$_}{'heter'});
			} ($sample)
		), 
	),"\n";

	close (OUT) ;


}

#open (OUT,">$dOut/sam.merge.snp.stat.xls") or die $!;
#
#print OUT join("\t",
#	'#Chro',
#	(map {("Snp_number($_)",
#	"Snp_number_gene($_)",
#	"Snp_number_intergenic($_)",
#	"Transition($_)",
#	"Transversion($_)",
#	"Heterozygosity($_)")} @sam), 
#),"\n";
#
#foreach my $chro (sort {$a cmp $b} keys %stat) {
#	next if ($chro eq 'Total') ;
#	map {
#		$stat{$chro}{$_}{'gene_region'}||=0;
#		$stat{$chro}{$_}{'intergenic'}||=0;
#	}@sam;
#
#	print OUT join("\t",
#		$chro, 
#		(
#			map {
#				($stat{$chro}{$_}{'nsnp'}, $stat{$chro}{$_}{'gene_region'},$stat{$chro}{$_}{'intergenic'},$stat{$chro}{$_}{'transition'},$stat{$chro}{$_}{'tranversion'},$stat{$chro}{$_}{'heter'});
#			} @sam
#		), 
#	),"\n";
#}
#map {
#	$stat{'Total'}{$_}{'gene_region'}||=0;
#	$stat{'Total'}{$_}{'intergenic'}||=0;
#}@sam;
#print OUT join("\t",
#	'Total', 
#	(
#		map {
#			($stat{'Total'}{$_}{'nsnp'}, $stat{'Total'}{$_}{'gene_region'},$stat{'Total'}{$_}{'intergenic'},$stat{'Total'}{$_}{'transition'},$stat{'Total'}{$_}{'tranversion'},$stat{'Total'}{$_}{'heter'});
#		} @sam
#	), 
#),"\n";
#
#close (OUT) ;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ==============================================================
# sub function
# ==============================================================
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Ma chouxian <macx\@biomarker.com.cn> 
Discription:

Usage:
  Options:
  -i	<files>	snp stat file, sep by whitespace, forced
  -s	<files>	merged snp file, forced
  -d	<str>	Directory where output file produced,optional,default [./]
  -h		Help

USAGE
	print $usage;
	exit;
}
