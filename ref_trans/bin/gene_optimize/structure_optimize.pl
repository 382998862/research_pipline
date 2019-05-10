#!/usr/nim/perl -w
use strict;
use Data::Dumper;
use File::Basename qw(basename dirname);

if (@ARGV != 3) {
	print "\n\tFunction: Compare gene_boundary files create gene structure optimizition Result\n\n";
	print "\tInfiles: Old_gene_boundary.list & merged_gene_boundary.list\n\n";
	print "\tUsage: perl gene_sturcture_optimize.pl <merged_gene_gff> <Old_gene_gff> <out_index>\n\n";
	exit;
}

my $od = &dirname($ARGV[2]);
`mkdir $od` unless (-d $od);

my %merged_gene;
&Track_gff ($ARGV[0],\%merged_gene);

#print Dumper %merged_gene;die;

my %old_gene;
&Track_gff ($ARGV[1],\%old_gene);

#print Dumper %old_gene;die;

################# bioinfor_pipeline.log ###################
&show_log2("step_7: Gene_structure_optimize analysis start.");
###########################################################

open OUT,">$ARGV[2].geneStructure.optimize.xls" || die $!;
print OUT "#GeneID\tLocus\tStrand\tSite\tOriginalRegion\tOptimizeRegion\n";
foreach my $chr (sort keys %merged_gene) {
	foreach my $gene (sort keys %{$merged_gene{$chr}}) {
		next if (!exists $old_gene{$chr}{$gene});
		next if (!defined $old_gene{$chr}{$gene}{First_E} || !defined $old_gene{$chr}{$gene}{Last_E}) ;

		my @new_gene_FE = sort {$a<=>$b} @{$merged_gene{$chr}{$gene}{First_E}};
		my @new_gene_LE = sort {$a<=>$b} @{$merged_gene{$chr}{$gene}{Last_E}};
		
		my @old_gene_FE = sort {$a<=>$b} @{$old_gene{$chr}{$gene}{First_E}};
		my @old_gene_LE = sort {$a<=>$b} @{$old_gene{$chr}{$gene}{Last_E}};

		my $gene_LOCU = "$chr:"."$merged_gene{$chr}{$gene}{Start}"."-"."$merged_gene{$chr}{$gene}{End}";
		if ($new_gene_FE[0] < $old_gene_FE[0] && $new_gene_FE[1] > $old_gene_FE[0]) {
			my $original_region="$old_gene_FE[0]"."-"."$old_gene_FE[1]";
			my $optimize_region="$new_gene_FE[0]"."-"."$old_gene_FE[1]";
			if ($merged_gene{$chr}{$gene}{Strand} eq "+") {
				print OUT "$gene\t$gene_LOCU\t$merged_gene{$chr}{$gene}{Strand}\t5'\t$original_region\t$optimize_region\n";
			}
			else {
				print OUT "$gene\t$gene_LOCU\t$merged_gene{$chr}{$gene}{Strand}\t3'\t$original_region\t$optimize_region\n";
			}
		}
		if ($new_gene_LE[1] > $old_gene_LE[1] && $new_gene_LE[0] < $old_gene_LE[1]) {
			my $original_region="$old_gene_LE[0]"."-"."$old_gene_LE[1]";
#			my $optimize_region="$old_gene_LE[0]"."-"."$new_gene_FE[1]";
#---------------------------------------- by Simon Young 2014-12-12 ----------------------------------------------------
			my $optimize_region="$old_gene_LE[0]"."-"."$new_gene_LE[1]";
			if ($merged_gene{$chr}{$gene}{Strand} eq "+") {
				print OUT "$gene\t$gene_LOCU\t$merged_gene{$chr}{$gene}{Strand}\t3'\t$original_region\t$optimize_region\n";
			}
			else {
				print OUT "$gene\t$gene_LOCU\t$merged_gene{$chr}{$gene}{Strand}\t5'\t$original_region\t$optimize_region\n";
			}
		}
	}
}
close OUT;

##############################
&show_log2("step_7: Gene_structure_optimize analysis finished.");
##############################

################# sub
sub Track_gff { #&Track_gff($gff,\%gene_info);
	my $gff = shift @_;
	my $gene_info = shift @_;

	my %gene_iso;
	my %exons;

	open IN,"$gff" || die $!;
	while (<IN>) {
		chomp;
		my $exon_type = 0;
		next if (/^$/ || /\#/);
		my @tmp=split /\t+/,$_;
		if ($tmp[2] eq "gene") {
			$tmp[8]=~/ID=([^;]+)/;
			$gene_info->{$tmp[0]}{$1}{Start}=$tmp[3];
			$gene_info->{$tmp[0]}{$1}{End}=$tmp[4];
			$gene_info->{$tmp[0]}{$1}{Strand}=$tmp[6];
		}
		if ($tmp[2] eq "mRNA") {
			$tmp[8]=~/ID=([^;]+).*Parent=([^;]+)/;
			$gene_iso{$2}{$1}=1;
		}
		if ($tmp[2] eq "exon") {
			$tmp[8]=~/Parent=([^;]+)/;
			push @{$exons{$1}},($tmp[3],$tmp[4]);
			$exon_type = 1;
		}
		if ($tmp[2] eq "CDS" && $exon_type == 0) {
			$tmp[8]=~/Parent=([^;]+)/;
			push @{$exons{$1}},($tmp[3],$tmp[4]);
		}
	}
	close IN;

	foreach my $chr (sort keys %{$gene_info}) {
		foreach my $gene (sort keys %{$gene_info->{$chr}}) {
			my ($gene_FE_start,$gene_FE_end,$gene_LE_start,$gene_LE_end);
			next if (scalar keys %{$gene_iso{$gene}} == 0) ; ## skip geno with no mRNA 

			foreach my $iso (sort keys %{$gene_iso{$gene}}) {
				my @exon = sort {$a<=>$b} @{$exons{$iso}};
				if (!defined $gene_FE_start) {
					$gene_FE_start = $exon[0];
					$gene_FE_end = $exon[1];
					$gene_LE_start = $exon[-2];
					$gene_LE_end = $exon[-1];
					next;
				}
				if ($exon[0] <= $gene_FE_start) {
					if ($exon[1] < $gene_FE_start || $exon[1] > $gene_FE_end) {
						$gene_FE_end = $exon[1];
					}
					$gene_FE_start = $exon[0];
				}
				if ($exon[-1] >= $gene_LE_end) {
					if ($exon[-2] < $gene_LE_start || $exon[-2] > $gene_LE_end) {
						$gene_LE_start = $exon[-2];
					}
					$gene_LE_end = $exon[-1];
				}
			}
			push @{$gene_info->{$chr}{$gene}{First_E}},($gene_FE_start,$gene_FE_end) if (defined $gene_FE_start and defined $gene_FE_end) ;
			push @{$gene_info->{$chr}{$gene}{Last_E}},($gene_LE_start,$gene_LE_end) if (defined $gene_LE_start and defined $gene_LE_end) ;
		}
	}
}
#############################################################################################################
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    $wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

# &show_log1("cmd")
sub show_log1()
{
    open LOG, ">$od/../../bioinfor_pipeline.log";
	my ($txt) = @_ ;
    my $time = time();
    my $Time = &sub_format_datetime(localtime($time));
    print LOG "$Time:\t$txt\n" ;
    return ($time) ;
	close LOG;
}

# &show_log2("cmd")
sub show_log2()
{
    open LOG, ">>$od/../../bioinfor_pipeline.log";
	my ($txt) = @_ ;
    my $time = time();
    my $Time = &sub_format_datetime(localtime($time));
    print LOG "$Time:\t$txt\n" ;
    return ($time) ;
	close LOG;
}
#############################################################################################################