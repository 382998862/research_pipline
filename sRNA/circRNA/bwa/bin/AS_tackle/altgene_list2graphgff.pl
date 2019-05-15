#!/usr/bin/perl -w
use strict;

if (@ARGV!=3) {
	print "\n\tFunction: Use splice gene list & cuffcmp.combined.gtf of cuffcompare Altgene graphgff for SpliceGrapher\n\n";
	print "\tUsage: perl cufftracking2splice_list.pl <altsplice_gene_list> <cuffcmp.combined.gtf> <graphgff_od>\n\n";
	exit;
}

my %Alt_gene;

open IN,"$ARGV[0]" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/);
	my ($gene,$trans)=split/\t+/,$_;
	$Alt_gene{$gene}{$trans}=1;
}
close IN;

my %Alt_gene_info;

open IN,"$ARGV[1]" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp=split/\t+/,$_;
	my $chr=$tmp[0];
	next if ($tmp[8]!~/gene_name\s\".*oId\s\"/);
	$tmp[8]=~/gene_name\s\"([^\"]+)(.*)oId\s\"([^\"]+)/;
	my $gene_name=$1;
	my $trans_name=$3;
	if (exists $Alt_gene{$gene_name}{$trans_name}) {
		$Alt_gene_info{$chr}{$gene_name}{$trans_name}{Strand}=$tmp[6];
		push @{$Alt_gene_info{$chr}{$gene_name}{$trans_name}{ExonArray}},($tmp[3],$tmp[4]);
	}
}
close IN;

############### gene & trans boundary define
my %gene_boundary;

foreach my $chr (sort keys %Alt_gene_info) {
	foreach my $gene (sort keys %{$Alt_gene_info{$chr}}) {
		my ($gene_start,$gene_end,$gene_strand);
		foreach my $trans (sort keys %{$Alt_gene_info{$chr}{$gene}}) {
			my @Exons=sort {$a<=>$b} @{$Alt_gene_info{$chr}{$gene}{$trans}{ExonArray}};
			$Alt_gene_info{$chr}{$gene}{$trans}{Start}=$Exons[0];
			$Alt_gene_info{$chr}{$gene}{$trans}{End}=$Exons[-1];
			if (!defined $gene_start) {
				$gene_start=$Exons[0];
				$gene_end=$Exons[-1];
				$gene_strand=$Alt_gene_info{$chr}{$gene}{$trans}{Strand};
			}
			else {
				$gene_start=$Exons[0] if ($Exons[0]<$gene_start);
				$gene_end=$Exons[-1] if ($Exons[-1]>$gene_end);
			}
		}
		$gene_boundary{$chr}{$gene}{Start}=$gene_start;
		$gene_boundary{$chr}{$gene}{End}=$gene_end;
		$gene_boundary{$chr}{$gene}{Strand}=$gene_strand;
	}
}

open OUT,">$ARGV[2]/Alt_gene.graphgff" || die $!;
foreach my $chr (sort keys %Alt_gene_info) {
	system "mkdir $ARGV[2]/$chr";
	foreach my $gene (sort keys %{$Alt_gene_info{$chr}}) {
		open OUT1,">$ARGV[2]/$chr/$gene.gff" || die $!;
		print OUT "$chr\tCufflinks\tgene\t$gene_boundary{$chr}{$gene}{Start}\t$gene_boundary{$chr}{$gene}{End}\t.\t$gene_boundary{$chr}{$gene}{Strand}\t.\tID=$gene;Nmae=$gene\n";
		print OUT1 "$chr\tCufflinks\tgene\t$gene_boundary{$chr}{$gene}{Start}\t$gene_boundary{$chr}{$gene}{End}\t.\t$gene_boundary{$chr}{$gene}{Strand}\t.\tID=$gene;Nmae=$gene\n";
		foreach my $trans (sort keys %{$Alt_gene_info{$chr}{$gene}}) {
			print OUT "$chr\tCufflinks\tmRNA\t$Alt_gene_info{$chr}{$gene}{$trans}{Start}\t$Alt_gene_info{$chr}{$gene}{$trans}{End}\t.\t$Alt_gene_info{$chr}{$gene}{$trans}{Strand}\t.\tID=$trans;Name=$trans;Parent=$gene\n";
			print OUT1 "$chr\tCufflinks\tmRNA\t$Alt_gene_info{$chr}{$gene}{$trans}{Start}\t$Alt_gene_info{$chr}{$gene}{$trans}{End}\t.\t$Alt_gene_info{$chr}{$gene}{$trans}{Strand}\t.\tID=$trans;Name=$trans;Parent=$gene\n";
			my @Exons=sort {$a<=>$b} @{$Alt_gene_info{$chr}{$gene}{$trans}{ExonArray}};
			for (my $i=0;$i<$#Exons;$i+=2) {
				print OUT "$chr\tCufflinks\tCDS\t$Exons[$i]\t$Exons[$i+1]\t.\t$Alt_gene_info{$chr}{$gene}{$trans}{Strand}\t.\tParent=$trans\n";
				print OUT "$chr\tCufflinks\texon\t$Exons[$i]\t$Exons[$i+1]\t.\t$Alt_gene_info{$chr}{$gene}{$trans}{Strand}\t.\tParent=$trans\n";
				print OUT1 "$chr\tCufflinks\tCDS\t$Exons[$i]\t$Exons[$i+1]\t.\t$Alt_gene_info{$chr}{$gene}{$trans}{Strand}\t.\tParent=$trans\n";
				print OUT1 "$chr\tCufflinks\texon\t$Exons[$i]\t$Exons[$i+1]\t.\t$Alt_gene_info{$chr}{$gene}{$trans}{Strand}\t.\tParent=$trans\n";
			}
		}
		close OUT1;
	}
}
close OUT;
