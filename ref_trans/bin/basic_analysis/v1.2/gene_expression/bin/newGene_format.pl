#!/usr/bin/perl -w
use strict;
use Data::Dumper;

if (@ARGV != 3) {
	print "\n\tFunction: Extract final gff from track list & merged.gtf\n\n";
	print "\tInfiles: merged.gtf & AAA.newGene.track.list.info\n\n";
	print "\tperl final_track_gff.pl <merged.gtf> <in_index> <out_dir>\n\n";
	exit;
}

my %iso_info;
my %track2chr;
my %known_genename;
=c
open IN,"$ARGV[0]" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp=split /\t+/,$_;

	next if ($tmp[8]!~/class_code\s\"u\"/);
	$tmp[8]=~/gene_id\s\"([^\"]+)\";\stranscript_id\s\"([^\"]+)\";\sexon_number\s\"(\d+)\";/;
	$track2chr{$tmp[0]}{$1}=1;
	push @{$iso_info{$1}{$2}{Exonarray}},($tmp[3],$tmp[4]);
	$iso_info{$1}{$2}{Strand}=$tmp[6];
}
close IN;
=cut     ##############20151207

open IN,"$ARGV[0]" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp=split /\t+/,$_;

	next if ($tmp[8]!~/class_code\s\"j\"/ && $tmp[8]!~/class_code\s\"=\"/);
	next unless /gene_name/;
    $tmp[8]=~/gene_id\s\"([^\"]+)\";\stranscript_id\s\"([^\"]+)\";.*gene_name \"([^\"]+)\";.*oId \"([^\"]+)\";.*class_code \"([^\"]+)\";/ ;
    push@{$known_genename{$1}},$3 ;
}
seek IN,0,0;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp = split/\t+/,$_;
	my @info = split/;/,$tmp[8];

	my ($gene_track,$geneID,$iso_track,$isoID,$class_code);
	if ($tmp[8]=~/class_code \"u\";/) {
		$tmp[8]=~/gene_id\s\"([^\"]+)\";\stranscript_id\s\"([^\"]+)\";\sexon_number.*class_code \"([^\"]+)\";/;
		$gene_track=$1;
		$iso_track=$2;
        $class_code=$3;
        unless (exists $known_genename{$gene_track}) {
            $track2chr{$tmp[0]}{$gene_track}=1;
            push @{$iso_info{$gene_track}{$iso_track}{Exonarray}},($tmp[3],$tmp[4]);
            $iso_info{$gene_track}{$iso_track}{Strand}=$tmp[6];
        }
	}
}
close IN;

my %new_gene_track2name;
my %new_iso_track2name;
my %new_gene_iso;

open IN,"$ARGV[1].newGene.track.list.info" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp=split /\t+/,$_;
	$new_gene_track2name{$tmp[0]}=$tmp[1];
	$new_iso_track2name{$tmp[2]}=$tmp[3];
	$new_gene_iso{$tmp[0]}{$tmp[2]}=1;
	my @exons=sort {$a<=>$b} @{$iso_info{$tmp[0]}{$tmp[2]}{Exonarray}};
	$iso_info{$tmp[0]}{$tmp[2]}{Start}=$exons[0];
	$iso_info{$tmp[0]}{$tmp[2]}{End}=$exons[-1];
}
close IN;

my %gene_boundary;
foreach my $gene_track (sort keys %new_gene_iso) {
	my ($gene_start,$gene_end,$gene_strand);
	foreach my $iso_track (sort keys %{$new_gene_iso{$gene_track}}) {
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

open OUT,">$ARGV[2]/$ARGV[1].newGene_final_tracking.list" || die $!;
open OUT1,">$ARGV[2]/$ARGV[1].newGene_final.gff" || die $!;

foreach my $chr (sort keys %track2chr) {
	foreach my $gene_track (sort keys %{$track2chr{$chr}}) {
		print OUT1 "$chr\tStringTie\tgene\t$gene_boundary{$gene_track}{Start}\t$gene_boundary{$gene_track}{End}\t.\t$gene_boundary{$gene_track}{Strand}\t.\tID=$new_gene_track2name{$gene_track}\n";
		foreach my $iso_track (sort keys %{$new_gene_iso{$gene_track}}) {
			print OUT1 "$chr\tStringTie\tmRNA\t$iso_info{$gene_track}{$iso_track}{Start}\t$iso_info{$gene_track}{$iso_track}{End}\t.\t$iso_info{$gene_track}{$iso_track}{Strand}\t.\tID=$new_iso_track2name{$iso_track}".";Parent=$new_gene_track2name{$gene_track}\n";
			my @exons = sort {$a<=>$b} @{$iso_info{$gene_track}{$iso_track}{Exonarray}};
			my $length;
			for (my $i=0 ; $i < $#exons ; $i+=2) {
				print OUT1 "$chr\tStringTie\tCDS\t$exons[$i]\t$exons[$i+1]\t.\t$iso_info{$gene_track}{$iso_track}{Strand}\t.\tParent=$new_iso_track2name{$iso_track}\n";
				$length += ($exons[$i+1] - $exons[$i] + 1);
			}
			print OUT "$chr\t$gene_track\t$new_gene_track2name{$gene_track}\t$gene_boundary{$gene_track}{Start}\t$gene_boundary{$gene_track}{End}\t$gene_boundary{$gene_track}{Strand}\t$iso_track\t$new_iso_track2name{$iso_track}\t$iso_info{$gene_track}{$iso_track}{Start}\t$iso_info{$gene_track}{$iso_track}{End}\t$length\n";
		}
	}
}
close OUT;
close OUT1;
