#!/usr/bin/perl -w
use strict;

if (@ARGV != 2)
{
	print "\n";
	print "    Function: Extract splice gene list from cuffcompare OUT cuffcmp.tracking\n\n";
	print "    Usage: perl cufftracking2splice_list.pl <cuffcmp.tracking> <altsplice_gene_list>\n\n";
	exit;
}

my %Alt_gene;

open (IN, $ARGV[0]) || die "Open $ARGV[0] failed!\n";
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	
	my ($cuff_trans_id, $cuff_locus_id, $ref_gene_trans, $class_node, $attribute) = split;
	next unless ($class_node =~ /j/); # j: Potentially novel isoform (fragment): at least one splice junction is shared with a reference transcript

	my ($gene_id, $main_trans) = $ref_gene_trans =~ /(.*)\|(.*)/;

	$attribute =~ s/^q\d://;
	my ($alt_gene, $alt_trans) = split /\|/, $attribute;
	
	$Alt_gene{$gene_id}{$main_trans}++;
	$Alt_gene{$gene_id}{$alt_trans}++;
}
close IN;

open OUT, ">", $ARGV[1] || die "Open $ARGV[1] failed!\n";
foreach my $gene (sort keys %Alt_gene) 
{
	foreach (sort keys %{$Alt_gene{$gene}}) 
	{
		print OUT join("\t", $gene, $_), "\n";
	}
}
close OUT;
