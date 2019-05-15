#!/usr/bin/perl
use strict;

if (@ARGV!=2) {
	print "\n\tFunction: Extract splice gene list from cuffcompare OUT cuffcmp.tracking\n\n";
	print "\tUsage: perl cufftracking2splice_list.pl <cuffcmp.tracking> <altsplice_gene_list>\n\n";
	exit;
}

my %Alt_gene;

open IN,"$ARGV[0]" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp=split/\t+/,$_;
	next if ($tmp[3]!~/j/);
	my ($gene_id,$main_trans,$alt_trans,$alt_gene);
	$tmp[2]=~/(.*)\|(.*)/;
	$gene_id=$1;
	$main_trans=$2;
	$tmp[4]=~s/^q\d://;
	my @alt_info=split/\|/,$tmp[4];
	$alt_gene=$alt_info[0];
	$alt_trans=$alt_info[1];
	$Alt_gene{$gene_id}{$main_trans}=1;
	$Alt_gene{$gene_id}{$alt_trans}=1;
}
close IN;

open OUT,">$ARGV[1]" || die $!;
foreach my $gene (sort keys %Alt_gene) {
	foreach (sort keys %{$Alt_gene{$gene}}) {
		print OUT "$gene\t$_\n";
	}
}
close OUT;