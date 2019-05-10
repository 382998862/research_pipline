#!/usr/local/bin/perl -w
use strict;

use Cwd 'abs_path';
use Getopt::Long;

my $usage = "
##*********************************************************************
  Destription: This script is used to extract result from vcf format.
  Usaage:
      -i    <str>      input file, vcf format..............[required]
      -q    <int>      quality value of snp probability....[default:20]
      -d    <int-int>  covrage depth range.................[default:2-20]
      -r    <int>      snps distance cutoff................[default:5]
      -key  <str>      output file prefix..................[required]
      -o    <str>      output dir..........................[required]
##*********************************************************************

";

my $snp_vcf;
my $snp_qual;
my $dis;
my $depth_range;
my $key;
my $out_dir;

GetOptions(
	"i:s" =>\$snp_vcf,  "q:i" =>\$snp_qual,  "d:s" =>\$depth_range,  "o:s" =>\$out_dir,	
	"key:s"=>\$key,     "r:i" =>\$dis,
);

die $usage unless ($snp_vcf);

$out_dir     ||= ".";      mkdir $out_dir unless (-d $out_dir);
$snp_qual    ||= 20;
$depth_range ||= "5-100";
$dis         ||= 5;
$snp_vcf = abs_path($snp_vcf);
$out_dir = abs_path($out_dir);

my ($min_depth, $max_depth) = $depth_range =~ /(\d+)-(\d+)/;

### Some definition
# $snp_type       homo, hetro
# $variance_type  transition, transversion

my %hetro_type = ('AG'=>'R', 'CT'=>'Y', 'GT'=>'K', 
	              'AC'=>'M', 'CG'=>'S', 'AT'=>'W',
				  'A' =>'A', 'C' =>'C', 'G' =>'G', 'T'=>'T');
my $pre_masked_vcf = "$out_dir/All_Sample.pre_masked.vcf";
my $filtered_vcf = "$out_dir/final.snp.vcf";
my $filtered_indel = "$out_dir/final.indel.vcf";
#my $snp_file     = "$out_dir/$key.snp";
#my $snp_stat     = "$out_dir/$key.stat.xls";
my %snp_type_num;
my %variance_type;
pre_mask();
filter();			  

		
exit;    
## end
###************************************###
## sub functions part
#

sub pre_mask
{## This sub function is used to mask those snp locus too closed.
	open (IN, $snp_vcf) || die "Open $snp_vcf failed!\n";
	open OUT, ">", $pre_masked_vcf;
	my $line_flag = 1;
	my $flag = 0;
	my $last_line;
	while (<IN>)
	{
		chomp;          next if (/^#/ || /^\s*$/);
		if($line_flag == 1)
		{
			$last_line = $_;     $line_flag = 0;    next;
		}
		my ($last_chr, $last_pos) = (split /\t+/, $last_line)[0, 1];
		my ($chr,      $pos     ) = (split /\t+/, $_)[0, 1];
		if($chr eq $last_chr && $pos - $last_pos <= $dis)
		{
			$flag = 1;
		}else
		{
#			print OUT $last_line."\n" unless $flag;     $flag = 0;
#-----------------------------by Simon Young 2015-05-05 -----------------------------
            unless ($flag) {
                print OUT $last_line."\n";
                print OUT $_."\n" if eof(IN);
            }
            $flag = 0;
		}
		$last_line = $_;
	}
	close OUT;
	close IN;
}

sub filter
{
	open (IN, $pre_masked_vcf) || die "Open $pre_masked_vcf failed!\n";
	open (VCF, ">", $filtered_vcf) || die "Open $filtered_vcf failed!\n";
	open (IND, ">", $filtered_indel) || die "Open $filtered_indel failed!\n";
	#open (OUT, ">", $snp_file) || die "Open $snp_file failed!\n";
	#print OUT join ("\t", "#Chr", "Pos", "Ref", "Alt", "Qual", "Depth", "Genotype"), "\n";
	while(<IN>)
	{
		chomp;
		if (/^#/ || /^\s*$/){
			print VCF "$_\n";
			print IND "$_\n";
		}
		my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $flag) = split /\t+/,$_,10;
		next if ($ref =~ /[^ACGT]/ || $alt =~ /[^ACGT]/);
		my @flag=split/\t+/,$flag;
		my $Qual=0;
		#my $len_ref=length($ref);
		#my $len_alt=length($alt);
		foreach(@flag){
			my $Qual1      = (split /:/, $flag)[-1];
			if ($Qual1>$Qual) {
                $Qual=$Qual1;
            }
		}
		
		#my $base_type = (split /\//, $flag)[0];
		#print "$info\n";
		my ($depth)   = $info =~ /DP=(\d+)/;
		
		next unless $depth >= $min_depth && $depth <= $max_depth;
		next unless $Qual  >= $snp_qual;
		if ($info=~/INDEL/) {
            print IND $_."\n";
        }
		else{
			print VCF $_."\n";
		}
        
	

		#$variance_type{$chr}{'snp'}++;
		#my $snp_genotype;
		#if($base_type eq '1')
		#{
		#	$snp_genotype = $alt;
		#}else
		#{
		#	my $alleles   = exists $hetro_type{$ref.$alt} ? ($ref.$alt) : ($alt.$ref);
		#	$snp_genotype = $hetro_type{$alleles};
		#	$variance_type{$chr}{'hetro'}++;
		#}
		#print OUT join ("\t", $chr, $pos, $ref, $alt, $Qual, $depth, $snp_genotype), "\n";

		#my $variance_type = variance_type($ref, $alt);
		#$variance_type{$chr}{$variance_type}++;
	}
	close IND;
	close VCF;
	close IN;
=c
	open (STAT, ">", $snp_stat) || die "Open $snp_stat failed!\n";
	print STAT join ("\t", "#Chr", "Snp_number", "Transition", "Transversion", "Hetrozygosity"), "\n";
	my ($total_snp_num, $total_transversion_num, $total_transition_num, $total_hetro_num);
	foreach my $chr (sort keys %variance_type)
	{
		$variance_type{$chr}{'snp'}         ||=0;
		$variance_type{$chr}{'hetro'}       ||=0;
		$variance_type{$chr}{'transversion'}||=0;
		$variance_type{$chr}{'transition'}  ||=0;
		my $snp_num          = $variance_type{$chr}{'snp'};            next unless $snp_num > 0;
		my $transversion_num = $variance_type{$chr}{'transversion'};
		my $transition_num   = $variance_type{$chr}{'transition'};
		my $hetro_num        = $variance_type{$chr}{'hetro'};
		my $transversion_ratio = sprintf "%.2f", $transversion_num/$snp_num * 100;
		my $transition_ratio   = sprintf "%.2f", $transition_num/$snp_num * 100;
		my $hetro_ratio        = sprintf "%.2f", $hetro_num/$snp_num * 100;
		print STAT join ("\t", $chr, $snp_num, $transition_ratio, $transversion_ratio, $hetro_ratio), "\n";

		$total_snp_num          += $snp_num;
		$total_transversion_num += $transversion_num;
		$total_transition_num   += $transition_num;
		$total_hetro_num        += $hetro_num;
	}
	
	my $total_transversion_ratio = sprintf "%.2f", $total_transversion_num/$total_snp_num * 100;
	my $total_transition_ratio   = sprintf "%.2f", $total_transition_num/$total_snp_num * 100;
	my $total_hetro_ratio        = sprintf "%.2f", $total_hetro_num/$total_snp_num * 100;
	print STAT join ("\t", "Total", $total_snp_num, $total_transition_ratio, $total_transversion_ratio, $total_hetro_ratio), "\n";

	close STAT;
=cut
}





sub variance_type
{
	my ($ref, $variance) = @_;
	my $type;
	my $genotype = $ref.$variance;
	if($genotype eq 'AG' || $genotype eq 'CT' || $genotype eq 'GA' || $genotype eq 'TC')
	{
		$type = 'transition';
	}else
	{
		$type = 'transversion';
	}
	return $type;
}


