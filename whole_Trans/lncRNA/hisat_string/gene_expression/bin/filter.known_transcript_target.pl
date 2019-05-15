#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use Data::Dumper;
use FindBin        qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd            'abs_path';

my $BEGIN_TIME=time();

# ------------------------------------------------------------------
my ($gff, $count, $fpkm, $known_trans_fa, $new_trans_fa , $out_fa, $out_count, $out_fpkm);
GetOptions(
				"help|?"      =>\&USAGE,
				"gff:s"	      =>\$gff,
				"count:s"     =>\$count,
				"fpkm:s"      =>\$fpkm,
				"kfa:s"       =>\$known_trans_fa,
				"nfa:s"       =>\$new_trans_fa,
				"out_fa:s"    =>\$out_fa,
				"out_count:s" =>\$out_count,
				"out_fpkm:s"  =>\$out_fpkm,
				) or &USAGE;
&USAGE unless ($gff and $known_trans_fa and $new_trans_fa and $out_fa and $out_count and $out_fpkm) ;
################################

$gff = abs_path($gff);
my %hash;
open IN,"$gff" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp = split/\t+/,$_;
	my @info = split/;/,$tmp[8];
	if($tmp[2] eq "gene"){
		if($info[0] =~ /ID=(.*)/){
			$hash{$1}=1;
		}
	}	
}

my %hc;
open(COUNT,"$count") or die $!;
open(OC,">$out_count") or die $!;
while(<COUNT>){
	chomp;
	if(/^#/){print OC "$_\n";next;}
	my @count = split /\t/,$_;
	my $total_count = 0;
	for (my $i=1;$i<@count;$i++){
		$total_count += $count[$i];
	}
	next if($total_count ==0);
	$hc{$count[0]}=1;
	print OC "$_\n";
}
close(COUNT);
close(OC);

open(FPKM,"$fpkm") or die $!;
open(OF,">$out_fpkm") or die $!;
while(<FPKM>){
	chomp;
	if(/^#/){print OF "$_\n";next;}
	my @fpkm = split /\t/,$_;
	my $total_fpkm =0;
	for (my $j=1;$j<@fpkm;$j++){
		$total_fpkm += $fpkm[$j];
	}
	if($total_fpkm ==0){
		if(exists $hc{$fpkm[0]}){
			print "total count of all samples is not zero and total fpkm is zero,please note!\n ";
		}
	}
	if(exists $hc{$fpkm[0]}){
		print OF "$_\n";
	}
}
close(FPKM);
close(OF);

my %kh;
$/=">";
open(IN1,"$known_trans_fa") or die $!;
<IN1>;
open(OUT,">$out_fa") or die $!;
while(<IN1>){
	chomp;
	my ($name,$fa)=split /\n/,$_,2;
	my ($id,$trans) = split /\s+/,$name,2;
	$kh{$id}=$fa;
	if(exists $hash{$id} && exists $hc{$id}){
		print OUT ">$id $trans\n$kh{$id}";
	}
}
$/="\n";

open(IN2,"$new_trans_fa") or die $!;
while(<IN2>){
	chomp;
	print OUT "$_\n";
}

close OUT;
close IN1;
close IN;



sub USAGE {#
	my $usage=<<"USAGE";

Usage:
	-gff	gene.gff3
	-count	All_gene_count.list
	-fpkm	All_gene_fpkm.list
	-kfa	known longest transcript fa
	-nfa	new longest transcript fa
	-out_fa	all longest transcript fa (expression > 0)
	-out_count	All_gene_count.list (expression > 0)
	-out_fpkm	All_gene_fpkm.list (expression > 0)

Example: perl $0 -gff gene.gff3 -count All_gene_count.list -fpkm All_gene_fpkm.list -kfa Known_longest_transcript.fa -nfa Human_New_longest_transcript.fa -out_count All_gene_count.list1 -out_fpkm All_gene_fpkm.list1

USAGE
	print $usage;
	exit;
}
