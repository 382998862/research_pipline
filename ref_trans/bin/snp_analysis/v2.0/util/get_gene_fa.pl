#!/usr/bin/perl -w
use strict;
use Data::Dumper;

if (@ARGV != 3) {
	print "\n\tFunction: Extract gene to transcript list from gff\n\n";
	print "\tUsage: <genome_fa> <gff_file> <out_index>\n\n";
	print "\tperl $0 genome_fa formated_gff out_index\n\n";
	exit;
}

my %gene_seq;

open IN,"$ARGV[0]" || die $!;
$/='>';
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my ($head,$seq) = split /\n+/,$_,2;
	my $id = (split /\s+/,$head)[0];
	$seq =~ s/\s+//g;
	$gene_seq{$id}=$seq;
}
$/="\n";
close IN;



open IN,"$ARGV[1]" || die $!;
open OUT,">$ARGV[2].fa" ||die $!;
while (<IN>) {
	chomp;
	next if (/^\*$/ || /^\#/);
	$_=~s/\s*$//;
	my @tmp = split /\t+/,$_;
	next unless ($tmp[2] =~ /gene/i);
	my $gene =(split(/ID=|;/,$tmp[8]))[1];
	my $start=($tmp[3]>$tmp[4])?$tmp[4]:$tmp[3];
	my $end=($tmp[3]<$tmp[4])?$tmp[4]:$tmp[3];
	my $len=$end-$start+1;
	my $gene_fa;
	if ($tmp[6] eq '+') {
		$gene_fa=substr($gene_seq{$tmp[0]},$start-1,$len);
	}
	else {
		$gene_fa=substr($gene_seq{$tmp[0]},$start-1,$len);
		$gene_fa =~tr/atgcATCGuU/tacgTAGCAA/;
		$gene_fa =reverse $gene_fa;
	}
	print OUT ">$gene\n$gene_fa\n";
	
}