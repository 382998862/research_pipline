#!/usr/nim/perl -w
use strict;
use Data::Dumper;


if (@ARGV != 3) {
	print "\n\tFunction: Compare gene_boundary files create gene structure optimizition Result\n\n";
	print "\tInfiles: AAA.final_tracking.list & AAA.newGene.final.list.info & Cuffdiff_dir & out_dir\n\n";
	print "\tUsage: perl list2expression.pl <final_list> <Cufflinks_dir> <Sample>\n\n";
	exit;
}

my %gene_track2name;
my %iso_track2name;
my %gene_iso;
my %gene_info;
my %iso_info;

open IN,"$ARGV[0]" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp=split /\t+/,$_;
	if ($tmp[1]=~/SERPINB10/) {
        print "$tmp[1]\n";
    }
    
	$gene_track2name{$tmp[1]}=$tmp[2];
	$iso_track2name{$tmp[6]}=$tmp[7];
	$gene_iso{$tmp[1]}{$tmp[6]}=1;
	$gene_info{$tmp[1]}{LOCU}="$tmp[0]:"."$tmp[3]"."-"."$tmp[4]";
	$gene_info{$tmp[1]}{Strand}=$tmp[5];
	$iso_info{$tmp[6]}{LOCU}="$tmp[0]:"."$tmp[8]"."-"."$tmp[9]";
	$iso_info{$tmp[6]}{Strand}=$tmp[5];
	$iso_info{$tmp[6]}{Length}=$tmp[-1];
}
close IN;

my %gene_length;
foreach my $gene (keys %gene_iso) {
	my $len = 0;
	foreach my $iso (keys %{$gene_iso{$gene}}) {
		if ($iso_info{$iso}{Length} > $len) {
			$len = $iso_info{$iso}{Length};
		}
	}
	$gene_length{$gene} = $len;
}

my %gene_fpkm;
open IN,"$ARGV[1]/genes.fpkm_tracking" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /\#/ || /^tracking_id/);
	my @tmp=split /\t+/,$_;
	$gene_fpkm{$tmp[0]}{FPKM}=$tmp[-4];
}
close IN;

my %iso_fpkm;
open IN,"$ARGV[1]/isoforms.fpkm_tracking" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /\#/ || /^tracking_id/);
	my @tmp=split /\t+/,$_;
	$iso_fpkm{$tmp[0]}{FPKM}=$tmp[-4];
	$iso_fpkm{$tmp[0]}{COV}=$tmp[-5];
}
close IN;

open OUT,">$ARGV[2].geneExpression.xls" || die $!;
print OUT "GeneID\tGene_position\tLength\tStrand\tfpkm\n";
foreach my $gene_track (sort keys %gene_info) {
	print OUT "$gene_track2name{$gene_track}\t$gene_info{$gene_track}{LOCU}\t$gene_length{$gene_track}\t$gene_info{$gene_track}{Strand}\t";
	print OUT "$gene_fpkm{$gene_track}{FPKM}\n";
}
close OUT;

open OUT,">$ARGV[2].isoExpression.xls" || die $!;
print OUT "isoformID\tisoform_position\tLength\tStrand\tcoverage\tfpkm\n";
foreach my $iso_track (sort keys %iso_info) {
	print OUT "$iso_track2name{$iso_track}\t$iso_info{$iso_track}{LOCU}\t$iso_info{$iso_track}{Length}\t$iso_info{$iso_track}{Strand}\t";
	print OUT "$iso_fpkm{$iso_track}{COV}\t$iso_fpkm{$iso_track}{FPKM}\n";
}
close OUT;