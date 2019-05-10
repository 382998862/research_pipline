#!/usr/bin/perl -w
use strict;
use Data::Dumper;

if (@ARGV != 3) {
	print "\n\tFunction: Extract final gff from track list & merged.gtf\n\n";
	print "\tInfiles: AAA.uniq_track2name.list.info & AAA.multi_track2name.list.info & AAA.multi_name2track.list.info\n\n";
	print "\tperl final_track_gff.pl <merged.gtf> <in_index> <out_dir>\n\n";
	exit;
}

my %iso_info;
my %track2chr;
my %chr2track;

open IN,"$ARGV[0]" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp=split /\t+/,$_;
	next if($tmp[2]!~/transcript/);
	next if ($tmp[8]!~/class_code\s\"j\"/ && $tmp[8]!~/class_code\s\"=\"/);
	next unless /gene_name/;
	$tmp[8]=~/transcript_id\s\"([^\"]+)\";\sgene_id\s\"([^\"]+)\";\s/;
	$track2chr{$tmp[0]}{$2}=1;
	$chr2track{$2}=$tmp[0];
	push @{$iso_info{$2}{$1}{Exonarray}},($tmp[3],$tmp[4]);
	$iso_info{$2}{$1}{Strand}=$tmp[6];
}
close IN;

my %gene_track2name;
my %iso_track2name;
my %gene_iso;

open IN,"$ARGV[1].uniq_track2name.list.info" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp=split /\t+/,$_;
	$gene_track2name{$tmp[0]}=$tmp[1];
	$iso_track2name{$tmp[2]}=$tmp[3];
	$gene_iso{$tmp[0]}{$tmp[2]}=1;
	my @exons=sort {$a<=>$b} @{$iso_info{$tmp[0]}{$tmp[2]}{Exonarray}};
	$iso_info{$tmp[0]}{$tmp[2]}{Start}=$exons[0];
	$iso_info{$tmp[0]}{$tmp[2]}{End}=$exons[-1];
}
close IN;

my %gene_model_divided;

open IN,"$ARGV[1].multi_track2name.list.info" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp=split /\t+/,$_;
	my $name;
	if ($tmp[1]=~/(.*)\-\d+/) {
		$name = $1;
	}
	else {
		$name = $tmp[1];
	}
	$gene_model_divided{$name}{$tmp[0]}=$tmp[1];
	$gene_track2name{$tmp[0]}=$tmp[1];
	$iso_track2name{$tmp[2]}=$tmp[3];
	$gene_iso{$tmp[0]}{$tmp[2]}=1;
	my @exons=sort {$a<=>$b} @{$iso_info{$tmp[0]}{$tmp[2]}{Exonarray}};
	$iso_info{$tmp[0]}{$tmp[2]}{Start}=$exons[0];
	$iso_info{$tmp[0]}{$tmp[2]}{End}=$exons[-1];
}
close IN;

my %gene_model_integration;

open IN,"$ARGV[1].multi_name2track.list.info" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp=split /\t+/,$_;
	$gene_model_integration{$tmp[0]}=$tmp[5];
	$gene_track2name{$tmp[0]}=$tmp[1];
	$iso_track2name{$tmp[2]}=$tmp[3];
	$gene_iso{$tmp[0]}{$tmp[2]}=1;
	my @exons=sort {$a<=>$b} @{$iso_info{$tmp[0]}{$tmp[2]}{Exonarray}};
	$iso_info{$tmp[0]}{$tmp[2]}{Start}=$exons[0];
	$iso_info{$tmp[0]}{$tmp[2]}{End}=$exons[-1];
}
close IN;

my %gene_boundary;
foreach my $gene_track (sort keys %iso_info) {
	my ($gene_start,$gene_end,$gene_strand);
	foreach my $iso_track (sort keys %{$iso_info{$gene_track}}) {
		if (!defined $gene_start) {
			$gene_start=$iso_info{$gene_track}{$iso_track}{Start};
			$gene_end=$iso_info{$gene_track}{$iso_track}{End};
			$gene_strand=$iso_info{$gene_track}{$iso_track}{Strand};
		}
		else {
			if (!defined $iso_info{$gene_track}{$iso_track}{Start}) {
				print "$gene_track\n$iso_track\n";
				die;
			}
			if ($iso_info{$gene_track}{$iso_track}{Start} < $gene_start) {
				$gene_start=$iso_info{$gene_track}{$iso_track}{Start};
			}
			if ($iso_info{$gene_track}{$iso_track}{End} > $gene_end) {
				$gene_end=$iso_info{$gene_track}{$iso_track}{End};
			}
		}
	}
	$gene_boundary{$gene_track}{Start}=$gene_start;
	$gene_boundary{$gene_track}{End}=$gene_end;
	$gene_boundary{$gene_track}{Strand}=$gene_strand;
}

open OUT,">$ARGV[2]/$ARGV[1].final_tracking.list" || die $!;
open OUT1,">$ARGV[2]/$ARGV[1].final.gff" || die $!;

foreach my $chr (sort keys %track2chr) {
	foreach my $gene_track (sort keys %{$track2chr{$chr}}) {
		print OUT1 "$chr\tStringTie\tgene\t$gene_boundary{$gene_track}{Start}\t$gene_boundary{$gene_track}{End}\t.\t$gene_boundary{$gene_track}{Strand}\t.\tID=$gene_track2name{$gene_track}\n";
		foreach my $iso_track (sort keys %{$gene_iso{$gene_track}}) {
			print OUT1 "$chr\tStringTie\tmRNA\t$iso_info{$gene_track}{$iso_track}{Start}\t$iso_info{$gene_track}{$iso_track}{End}\t.\t$iso_info{$gene_track}{$iso_track}{Strand}\t.\tID=$iso_track2name{$iso_track}".";Parent=$gene_track2name{$gene_track}\n";
			my @exons = sort {$a<=>$b} @{$iso_info{$gene_track}{$iso_track}{Exonarray}};
			my $length_0;
			for (my $i=0 ; $i < $#exons ; $i+=2) {
				print OUT1 "$chr\tStringTie\tCDS\t$exons[$i]\t$exons[$i+1]\t.\t$iso_info{$gene_track}{$iso_track}{Strand}\t.\tParent=$iso_track2name{$iso_track}\n";
				$length_0 += ($exons[$i+1] - $exons[$i] + 1);
			}
			print OUT "$chr\t$gene_track\t$gene_track2name{$gene_track}\t$gene_boundary{$gene_track}{Start}\t$gene_boundary{$gene_track}{End}\t$gene_boundary{$gene_track}{Strand}\t$iso_track\t$iso_track2name{$iso_track}\t$iso_info{$gene_track}{$iso_track}{Start}\t$iso_info{$gene_track}{$iso_track}{End}\t$length_0\n";
		}
	}
}
close OUT;
close OUT1;

open OUT2,">$ARGV[2]/$ARGV[1].gene_model_divided.gff" || die $!;
foreach my $gene_name (sort keys %gene_model_divided) {
	foreach my $gene_track (sort keys %{$gene_model_divided{$gene_name}}) {
		print OUT2 "$chr2track{$gene_track}\tStringTie\tgene\t$gene_boundary{$gene_track}{Start}\t$gene_boundary{$gene_track}{End}\t.\t$gene_boundary{$gene_track}{Strand}\t.\tID=$gene_track2name{$gene_track}".";RefID=$gene_name\n";
		foreach my $iso_track (sort keys %{$gene_iso{$gene_track}}) {
			print OUT2 "$chr2track{$gene_track}\tStringTie\tmRNA\t$iso_info{$gene_track}{$iso_track}{Start}\t$iso_info{$gene_track}{$iso_track}{End}\t.\t$iso_info{$gene_track}{$iso_track}{Strand}\t.\tID=$iso_track2name{$iso_track}".";Parent=$gene_track2name{$gene_track}\n";
			my @exons = sort {$a<=>$b} @{$iso_info{$gene_track}{$iso_track}{Exonarray}};
			for (my $i=0 ; $i < $#exons ; $i+=2) {
				print OUT2 "$chr2track{$gene_track}\tStringTie\tCDS\t$exons[$i]\t$exons[$i+1]\t.\t$iso_info{$gene_track}{$iso_track}{Strand}\t.\tParent=$iso_track2name{$iso_track}\n";
			}
		}
	}
}
close OUT2;

open OUT3,">$ARGV[2]/$ARGV[1].gene_model_integrated.gff" || die $!;
foreach my $gene_track (sort keys %gene_model_integration) {
	print OUT3 "$chr2track{$gene_track}\tStringTie\tgene\t$gene_boundary{$gene_track}{Start}\t$gene_boundary{$gene_track}{End}\t.\t$gene_boundary{$gene_track}{Strand}\t.\tID=$gene_track2name{$gene_track}".";RefID=$gene_model_integration{$gene_track}\n";
	foreach my $iso_track (sort keys %{$gene_iso{$gene_track}}) {
		print OUT3 "$chr2track{$gene_track}\tStringTie\tmRNA\t$iso_info{$gene_track}{$iso_track}{Start}\t$iso_info{$gene_track}{$iso_track}{End}\t.\t$iso_info{$gene_track}{$iso_track}{Strand}\t.\tID=$iso_track2name{$iso_track}".";Parent=$gene_track2name{$gene_track}\n";
		my @exons = sort {$a<=>$b} @{$iso_info{$gene_track}{$iso_track}{Exonarray}};
		for (my $i=0 ; $i < $#exons ; $i+=2) {
			print OUT3 "$chr2track{$gene_track}\tStringTie\tCDS\t$exons[$i]\t$exons[$i+1]\t.\t$iso_info{$gene_track}{$iso_track}{Strand}\t.\tParent=$iso_track2name{$iso_track}\n";
		}
	}
}
close OUT3;
