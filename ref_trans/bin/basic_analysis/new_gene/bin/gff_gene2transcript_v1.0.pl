#!/usr/bin/env perl
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

my %gene2trans;
my %trans_info;

open IN, $ARGV[1] || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp = split /\t+/,$_;
	if ($tmp[2] =~ /mRNA/) {
		$tmp[8] =~ /.*ID=([^;\s\"]+).*Parent=([^;\s\"]+)/;
		$gene2trans{$tmp[0]}{$2}{$1} = 1;
		$trans_info{$1}{Start} = $tmp[3];
		$trans_info{$1}{End} = $tmp[4];
		$trans_info{$1}{Strand} = $tmp[6];
	}
	if ($tmp[2] =~ /exon/) {
		$tmp[8] =~ /.*Parent=([^;\s\"]+)/;
		push @{$trans_info{$1}{Exonarray}},($tmp[3],$tmp[4]);
	}
}
close IN;

my %gene_longest_trans;

open OUT1,">$ARGV[2].transcript.fa" || die $!;

foreach my $chr (sort keys %gene2trans) {
	foreach my $gene (sort keys %{$gene2trans{$chr}}) {
		my ($trans_max_len,$trans_max_id,$trans_max_seq);
		foreach my $iso (sort keys %{$gene2trans{$chr}{$gene}}) {
			my ($trans_len,$trans_seq);
			my @exons = sort {$a<=>$b} @{$trans_info{$iso}{Exonarray}};
			for (my $i=0 ; $i<$#exons ; $i+=2) {
				my $cds_len = ($exons[$i+1] - $exons[$i] + 1);
				$trans_len += $cds_len;
				$trans_seq .= substr ($gene_seq{$chr},$exons[$i] - 1,$cds_len);
			}
			if ($trans_info{$iso}{Strand} eq "-") {
				$trans_seq=~tr/ATCGatcg/TAGCtagc/;
				$trans_seq = reverse $trans_seq;
			}
			print OUT1 ">$iso\n$trans_seq\n";
		}
	}
}
close OUT1;
