#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use Data::Dumper;
use FindBin        qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd            'abs_path';

my $BEGIN_TIME=time();

# ------------------------------------------------------------------
my ($gtf, $out, $index,$track,$out1,$out2);
GetOptions(
				"help|?"      =>\&USAGE,
				"gtf:s"       =>\$gtf,
				"track:s"	=>\$track,
				"out1:s"	=>\$out1,
				"out2:s"	=>\$out2,

				) or &USAGE;
&USAGE unless ($gtf and  $track and $out1 and $out2 ) ;
$gtf = abs_path($gtf);
my %iso_info;
my %track2chr;

open IN,"$gtf" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp=split /\t+/,$_;
	if ($tmp[2]=~/transcript/ and $tmp[8]=~/class_code\s\"u\"/) {
        $tmp[8]=~/transcript_id\s\"([^\"]+)\";\sgene_id\s\"([^\"]+)\";\sxloc/;
		$track2chr{$tmp[0]}{$2}=1;
		$iso_info{$2}{$1}{Strand}=$tmp[6];
    }
	if ($tmp[2] =~/exon/) {
        $tmp[8]=~/transcript_id\s\"([^\"]+)\";\sgene_id\s\"([^\"]+)\";\sexon_number/;
		if (exists $iso_info{$2}{$1}) {
            push @{$iso_info{$2}{$1}{Exonarray}},($tmp[3],$tmp[4]);
        }
        
    }	
}
close IN;

my %new_gene_track2name;
my %new_iso_track2name;
my %new_gene_iso;

#open IN,"$ARGV[1].newGene.track.list.info" || die $!;
open IN,"$track" || die $!;
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

#open OUT,">$ARGV[2]/$ARGV[1].newGene_final_tracking.list" || die $!;
open OUT,">$out1" || die $!;
#open OUT1,">$ARGV[2]/$ARGV[1].newGene_final.gff" || die $!;
open OUT1,">$out2" || die $!;
print OUT1 "#Seq_ID\tSource\tType\tStart\tEnd\tScore\tStrand\tPhase\tAttributes\n";
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
sub USAGE {#
	my $usage=<<"USAGE";

Usage:
	-gtf              Compare.gtf fil_trackle       must be given;
	-track			human.newGene.track.list.info
	-out1			human,newGene_final_tracking.list
	
	-out2           Human.newGene_final.gff     must be given  ;
	
USAGE
	print $usage;
	exit;
}
