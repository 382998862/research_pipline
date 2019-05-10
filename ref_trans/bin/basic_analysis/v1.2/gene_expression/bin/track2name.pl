#!/usr/bin/env perl
use strict;
use Data::Dumper;

if (@ARGV != 2) {
	print "\n\tFunction: Extract geneID & isoformID to tracking ID list\n\n";
	print "\tUsage: perl track2name.pl <merge.gtf> <out_index>\n\n";
	exit;
}

my %gene_track2name;
my %iso_track2name;
my %newGene_track;

open IN,"$ARGV[0]" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp = split/\t+/,$_;
	next if ($tmp[2]!~/transcript/);
	my @info = split/;/,$tmp[8];
	my ($gene_track,$geneID,$iso_track,$isoID,$class_code);

	next if ($tmp[8]!~/class_code\s\"j\"/ && $tmp[8]!~/class_code\s\"=\"/);
	next unless /gene_name/;
	$tmp[8]=~/transcript_id\s\"([^\"]+)\";\sgene_id\s\"([^\"]+)\";\sgene_name \"([^\"]+)\";.*cmp_ref \"([^\"]+)\";.*class_code \"([^\"]+)\";/;
	$gene_track=$2;
	$iso_track=$1;
	$geneID=$3;
	$isoID=$4;
	$class_code=$5;

	$gene_track2name{$gene_track}{$geneID}{$iso_track}=1;
	$iso_track2name{$iso_track}{Class}=$class_code;
	$iso_track2name{$iso_track}{oId}=$isoID;
}
close IN;


open OUT1,">$ARGV[1].track2name.list" || die $!;
foreach my $gene_track (sort keys %gene_track2name) {
	foreach my $gene_name (sort keys %{$gene_track2name{$gene_track}}) {
		foreach my $iso_track (sort keys %{$gene_track2name{$gene_track}{$gene_name}}) {
			print OUT1 "$gene_track\t$gene_name\t$iso_track\t$iso_track2name{$iso_track}{oId}\t$iso_track2name{$iso_track}{Class}\n";
		}
	}
}
close OUT1;

