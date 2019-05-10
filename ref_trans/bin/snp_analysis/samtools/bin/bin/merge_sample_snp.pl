#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;

my (@fIn, $gff, $fKey, $dOut, $log);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s{,}" =>\@fIn,
				"gff:s"  =>\$gff,
				"d:s"    =>\$dOut,
				) or &USAGE;
&USAGE unless (@fIn and $gff);

$dOut ||= "./";     mkdir $dOut unless (-d $dOut);
$dOut   = abs_path($dOut);

my %snp = ();
my %sampleLabel = ();

#===============================================================
# Get Data
#===============================================================

foreach my $fIn (@fIn) {
	my ($sample) = basename($fIn) =~/(\w+)\.snp$/;
	$sampleLabel{$sample} = 1;
	
	open (IN,"$fIn") or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/) ;
		my ($chro, $pos, $refBase, @snpInfo) = split;

		my $pos_int = pos_unit($pos);
		push @{$snp{$chro}{$pos_int}{'pos_unit'}}, $pos;
		
		#@snpInfo: REF     Sample_Base     Qual    Depth   Genotype 
		$snp{$chro}{$pos}{'sam'}{$sample} = \@snpInfo;
		$snp{$chro}{$pos}{'processed'}    = 0;
		$snp{$chro}{$pos}{'refBase'}      = $refBase;
	}
	close (IN);
}

my @sample_list = sort {$a cmp $b} keys %sampleLabel;

#===============================================================
# Process
#===============================================================
open (OUT,">$dOut/sam.merge.snp") or die $!;
$|=1;
print OUT join("\t",
	'#Chro',
	'Pos',
	'Genic/Intergenic',
	'RefBase',
	(map {
		($_."_Base", $_."_Qual", $_."_Depth", $_."_Genotype");
	} @sample_list),
),"\n";

open (IN, $gff) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/) ;

	my ($chro, undef, $feture, $start, $end, undef, $strand, undef, $attribute) = split;
	next if ($feture ne "gene") ;

	my ($geneID) = $attribute =~/ID=([^;]+)/;
	($start, $end) = sort {$a <=> $b} ($start, $end);
	my @pos = find_snps_in_gene($chro, $start, $end, \%snp);

#	print join("\t",@pos),"\n";
#	<STDIN>;
	
	next unless (@pos) ; ## no snp found in this gene 
	
	foreach my $snpPos (@pos) {
		print OUT join("\t",
			$chro,
			$snpPos,
			$geneID,
			$snp{$chro}{$snpPos}{'refBase'},
			(map {
				if (defined $snp{$chro}{$snpPos}{'sam'}{$_}){
					@{$snp{$chro}{$snpPos}{'sam'}{$_}};
				}else{
					(".", ".", ".", ".");
				}
			} @sample_list),
		),"\n";
		
		$snp{$chro}{$snpPos}{'processed'} = 1;
	}
}

#
# output intergenic snps 
#
foreach my $chro (keys %snp) {
	foreach my $snpPos (keys %{$snp{$chro}}) {
		next unless exists $snp{$chro}{$snpPos}{'processed'};
		next if ($snp{$chro}{$snpPos}{'processed'} == 1) ;
		print OUT join("\t",
			$chro,
			$snpPos,
			'intergenic',
			$snp{$chro}{$snpPos}{'refBase'},
			(map {
				if (defined $snp{$chro}{$snpPos}{'sam'}{$_}){
					@{$snp{$chro}{$snpPos}{'sam'}{$_}};
				}else{
					(".", ".", ".", ".");
				}
			} @sample_list),
		),"\n";
	}
}

close (IN) ;
close (OUT) ;

# ==============================================================
# sub function
# ==============================================================
sub find_snps_in_gene {#
	my ($chro, $start, $end, $ref_snp) = @_;

	my $start_int = pos_unit($start);
	my $end_int   = pos_unit($end);
	my %tmp;
	for(my $i=$start_int; $i<=$end_int; $i++)
	{
		next unless exists $snp{$chro}{$i}{'pos_unit'};
		my @Pos = @{$snp{$chro}{$i}{'pos_unit'}};
		foreach (@Pos)
		{
			$tmp{$_}++;
		}
	}
	my @sortPos = sort {$a<=>$b} keys %tmp;
	
	my @res = ();

	foreach my $snpPos (@sortPos) {
		last if ($snpPos > $end) ;
		next if ($snpPos < $start) ;
		push @res, $snpPos;
	}

	return @res;
}

sub pos_unit
{
	my ($pos) = @_;
	my $unit  = 1e4;
	my $pos_int = int($pos/$unit);
	return $pos_int;
}

sub USAGE {#
	my $usage=<<"USAGE";
Usage:
  -i	<files>	SNP files, sep by whitespace, snp format, forced
  -gff	<file>	GFF file, forced
  -d	<str>	Directory where output file produced,optional,default [./]
  -h		Help

USAGE
	print $usage;
	exit;
}
