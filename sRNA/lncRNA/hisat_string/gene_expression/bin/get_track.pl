#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use Data::Dumper;
use FindBin        qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd            'abs_path';

my $BEGIN_TIME=time();

# ------------------------------------------------------------------
my ($gtf, $out, $index);
GetOptions(
				"help|?"      =>\&USAGE,
				"gtf:s"       =>\$gtf,
				"out:s"        =>\$out,
				"index:s"     =>\$index,
				) or &USAGE;
&USAGE unless ($gtf and  $out ) ;
################################
my $notename=`hostname`;chomp $notename;

$gtf = abs_path($gtf);
my %newGene_track;
open IN,"$gtf" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp = split/\t+/,$_;
	my @info = split/;/,$tmp[8];
	my ($gene_track,$geneID,$iso_track,$isoID,$class_code);
	if ($tmp[8]=~/class_code \"u\";/) {
		$tmp[8]=~/transcript_id\s\"([^\"]+)\";\sgene_id\s\"([^\"]+)\";\sxloc/;
		$gene_track=$2;
		$iso_track=$1;
		$newGene_track{$gene_track}{$iso_track}=1;
		next;
	}
}
close IN;
open OUT,">$out" || die $!;
my $i=1;
foreach my $gene_track (sort keys %newGene_track) {
	my $gene_name = "$index"."_newGene"."_$i";
	my $j = 1;
	foreach my $iso_track (sort keys %{$newGene_track{$gene_track}}) {
		my $iso_name = "$gene_name".".$j";
		print OUT "$gene_track\t$gene_name\t$iso_track\t$iso_name\n";
		$j++;
	}
	$i++;
}
close OUT;




sub USAGE {#
	my $usage=<<"USAGE";

Usage:
	-gtf              Compare.gtf fil_trackle       must be given
	-out               Human.newGene.track.list.info     must be given  ;
	
USAGE
	print $usage;
	exit;
}
