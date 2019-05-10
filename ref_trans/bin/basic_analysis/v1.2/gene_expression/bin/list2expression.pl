#!/usr/nim/perl -w
use strict;
use Data::Dumper;


if (@ARGV != 3) {
	print "\n\tFunction: Compare gene_boundary files create gene structure optimizition Result\n\n";
	print "\tInfiles: AAA.final_tracking.list & AAA.newGene.final.list.info & Cuffdiff_dir & out_dir\n\n";
	print "\tUsage: perl list2expression.pl <final_list> <Cuff_dir> <geneExpression>\n\n";
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

my %iso_fpkm;
my %iso_external_count;
my %iso_internal_count;
my %iso_read_group;
my %Sam;

open IN,"$ARGV[1]/isoforms.read_group_tracking" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /\#/ || /^tracking_id/);
	my @tmp=split /\t+/,$_;
	$Sam{$tmp[1]}=1;
	$iso_read_group{$tmp[1]}{$tmp[0]}=$tmp[-6];
	$iso_internal_count{$tmp[1]}{$tmp[0]}=$tmp[-5];
	$iso_external_count{$tmp[1]}{$tmp[0]}=$tmp[-4];
	$iso_fpkm{$tmp[1]}{$tmp[0]}=$tmp[-3];
}
close IN;

my %gene_fpkm;
my %gene_external_count;
my %gene_internal_count;
my %gene_read_group;
foreach my $sam (sort keys %Sam) {
	foreach my $gene_track (sort keys %gene_iso) {
		my ($fpkm,$external_count,$internal_count,$read_group);
		foreach my $iso_track (sort keys %{$gene_iso{$gene_track}}) {
			if (!defined $iso_fpkm{$sam}{$iso_track}) {
				print "$sam\t$gene_track\t$iso_track\n";
				delete $iso_info{$iso_track};
				next;
			}
			$fpkm+=$iso_fpkm{$sam}{$iso_track};
			$external_count+=$iso_external_count{$sam}{$iso_track};
			$internal_count+=$iso_internal_count{$sam}{$iso_track};
			$read_group+=$iso_read_group{$sam}{$iso_track};
		}
		next if (!defined $fpkm);
		$gene_fpkm{$sam}{$gene_track}=$fpkm;
		$gene_external_count{$sam}{$gene_track}=$external_count;
		$gene_internal_count{$sam}{$gene_track}=$internal_count;
		if ($read_group=~/\.\d+/) {
			$read_group=int($read_group+1);
		}
		$gene_read_group{$sam}{$gene_track}=$read_group;
	}
}

foreach my $sam (sort keys %Sam) {
	open OUT,">$ARGV[2]/$sam.geneExpression.xls" || die $!;
	print OUT "#GeneID\tLocus\tLength\tStrand\traw_fragment_number\tnormalized_count\tFPKM\n";
	foreach my $gene_track (sort keys %{$gene_fpkm{$sam}}) {
		print OUT "$gene_track2name{$gene_track}\t$gene_info{$gene_track}{LOCU}\t$gene_length{$gene_track}\t$gene_info{$gene_track}{Strand}\t";
		print OUT "$gene_read_group{$sam}{$gene_track}\t$gene_external_count{$sam}{$gene_track}\t$gene_fpkm{$sam}{$gene_track}\n";
	}
	close OUT;

	open OUT,">$ARGV[2]/$sam.isoExpression.xls" || die $!;
	print OUT "#GeneID\tLocus\tLength\tStrand\traw_fragment_number\tnormalized_count\tFPKM\n";
	foreach my $iso_track (sort keys %iso_info) {
		print OUT "$iso_track2name{$iso_track}\t$iso_info{$iso_track}{LOCU}\t$iso_info{$iso_track}{Length}\t$iso_info{$iso_track}{Strand}\t";
		print OUT "$iso_read_group{$sam}{$iso_track}\t$iso_external_count{$sam}{$iso_track}\t$iso_fpkm{$sam}{$iso_track}\n";
	}
	close OUT;
}