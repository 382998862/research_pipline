#!/usr/nim/perl -w
use strict;
use Data::Dumper;


if (@ARGV != 3) {
	print "\n\tFunction: Extract Gene Expression \n\n";
	print "\tInfiles: All.final_tracking.list &  Cuffnorm_dir & out_dir\n\n";
	print "\tUsage: perl list2expression.pl <final_list> <Cuff_dir> <geneExpression>\n\n";
	exit;
}

my %gene_track2name;
my %iso_track2name;
my %gene_iso;
my %gene_info;
my %iso_info;
my %iso_in;
open IN,"$ARGV[0]" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp=split /\t+/,$_;
	$gene_track2name{$tmp[1]}=$tmp[2];  #############
	$iso_track2name{$tmp[6]}=$tmp[7];   #############
	$gene_iso{$tmp[1]}{$tmp[6]}=1;
	$iso_in{$tmp[6]}=$tmp[2];
#	$gene_info{$tmp[1]}{LOCU}="$tmp[0]:"."$tmp[3]"."-"."$tmp[4]";
	$gene_info{$tmp[1]}{Strand}=$tmp[5]; #################
#	$iso_info{$tmp[6]}{LOCU}="$tmp[0]:"."$tmp[8]"."-"."$tmp[9]";
	$iso_info{$tmp[6]}{Strand}=$tmp[5];  ##################
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
my %iso_count;
my %iso_locus;
my %gene_fpkm;
my %gene_count;
my %gene_locus;
my @samples;
#my %iso_fpkm;
#my %iso_external_count;
#my %iso_internal_count;
#my %iso_read_group;
my %Sam;

open IN,"$ARGV[1]/genes.attr_table" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /\#/ || /^tracking_id/);
	my @tmp=split /\t+/,$_;
	$gene_locus{$tmp[0]}=$tmp[6];

}
close IN;

open IN,"$ARGV[1]/isoforms.attr_table" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /\#/ || /^tracking_id/);
	my @tmp=split /\t+/,$_;
	$iso_locus{$tmp[0]}=[$tmp[6],$tmp[7]];
}
close IN;

open IN,"$ARGV[1]/genes.count_table" || die $!;
while (<IN>) {
	chomp;
	if (/^tracking_id/){
		$_=~s/_0//g;
		@samples=split(/\t+/,$_);
		shift@samples;
		next;
	}
	my @tmp=split /\t+/,$_;
	my $gene=shift@tmp;
	for(my $i=0;$i<=$#tmp;$i++){
		$gene_count{$samples[$i]}{$gene}=$tmp[$i];
	}
}
close IN;

open IN,"$ARGV[1]/genes.fpkm_table" || die $!;
open OUT,">$ARGV[2]/All_geneExpression.list" || die $!;
while (<IN>) {
	chomp;
	if (/^tracking_id/){
		$_=~s/_0//g;
		#print OUT $_,"\n";
		@samples=split(/\t+/,$_);
		shift@samples;
		print OUT "GeneID\t",join("\t",@samples),"\n";
		next;
	}
	my @tmp=split /\t+/,$_;
	my $gene=shift@tmp;
	print OUT $gene_track2name{$gene},"\t",join("\t",@tmp),"\n";
	for(my $i=0;$i<=$#tmp;$i++){
		$gene_fpkm{$samples[$i]}{$gene}=$tmp[$i];
	}
}
close IN;
close OUT;

open IN,"$ARGV[1]/isoforms.count_table" || die $!;
while (<IN>) {
	chomp;
	if (/^tracking_id/){
		$_=~s/_0//g;
		@samples=split(/\t+/,$_);
		shift@samples;
		next;
	}
	my @tmp=split /\t+/,$_;
	my $transcript=shift@tmp;
	for(my $i=0;$i<=$#tmp;$i++){
		$iso_count{$samples[$i]}{$transcript}=$tmp[$i];
	}
}
close IN;
open IN,"$ARGV[1]/isoforms.fpkm_table" || die $!;
open OUT,">$ARGV[2]/All_isoExpression.list" || die $!;
while (<IN>) {
	chomp;
	if (/^tracking_id/){
		$_=~s/_0//g;
		#print OUT $_,"\n";
		@samples=split(/\t+/,$_);
		shift@samples;
		print OUT "Transcript\t",join("\t",@samples),"\n";
		next;
	}
	my @tmp=split /\t+/,$_;
	my $transcript=shift@tmp;
	print OUT $iso_track2name{$transcript},"\t",join("\t",@tmp),"\n";
	for(my $i=0;$i<=$#tmp;$i++){
		$iso_fpkm{$samples[$i]}{$transcript}=$tmp[$i];
	}
}
close IN;
close OUT;

foreach my $sam (sort @samples) {
	open OUT1,">$ARGV[2]/$sam.geneExpression.xls" || die $!;
	print OUT1 "#GeneID\tLength\tFPKM\tLocus\tStrand\tNormalized_Count\n";
	open OUT2,">$ARGV[2]/$sam.isoExpression.xls" || die $!;
	print OUT2 "#Transcript_id\tLength\tFPKM\tLocus\tStrand\tGeneID\tNormalized_Count\n";
	foreach my $gene (sort keys $gene_fpkm{$sam}) {
		print OUT1 "$gene_track2name{$gene}\t$gene_length{$gene}\t$gene_fpkm{$sam}{$gene}\t$gene_locus{$gene}\t$gene_info{$gene}{Strand}\t$gene_count{$sam}{$gene}\n";
	}
	foreach my $transcript (sort keys $iso_fpkm{$sam}) {
		print OUT2 "$iso_track2name{$transcript}\t$iso_locus{$transcript}->[1]\t$iso_fpkm{$sam}{$transcript}\t$iso_locus{$transcript}->[0]\t$iso_info{$transcript}{Strand}\t$iso_in{$transcript}\t$iso_count{$sam}{$transcript}\n";
	}
	close OUT1;
	close OUT2;
}
