#!/usr/bin/env perl
use strict;
use FindBin qw($Bin $Script);
use Cwd 'abs_path';
use Getopt::Long;
use Data::Dumper;
my ($species,$genome,$new_dir,$singleexon);
GetOptions(
				"help|?"      =>\&USAGE,
				"s:s"       =>\$species,
				"g:s"    =>\$genome,
				"newdir:s"  =>\$new_dir,
				"singleexon:s"=>\$singleexon,
				) or &USAGE;
&USAGE unless ($species and  $genome and $new_dir) ;
$genome  = abs_path($genome);
$new_dir = abs_path($new_dir);

my (%filtered_1, %filtered);

chdir($new_dir);
`perl $Bin/gff_gene2transcript_v1.0.pl $genome $new_dir/$species.newGene_final.gff $species.newGene`;
print "perl $Bin/gff_gene2transcript_v1.0.pl $genome $new_dir/$species.newGene_final.gff $species.newGene \n";
my $fa_file = "$new_dir/$species.newGene.transcript.fa";
mkdir("tmp_orf") unless -d 'tmp_orf';
`perl $Bin/Getorf_Extract.pl -fa $fa_file -od $new_dir/tmp_orf`;
open (IN,"$new_dir/tmp_orf/$species.newGene.transcript.orf") || die $!;
my %len;
while(<IN>)
{
	chomp;
	next unless /^>/;
	my $trans_id = (split /\s+/)[0];
	$trans_id =~ />(\S+)_(\d+)$/;
	$trans_id = $1;
	my $seq = <IN>;             chomp $seq;
	my $len = length($seq);
	$len{$trans_id} = $len unless exists $len{$trans_id};
	$len{$trans_id} = ($len > $len{$trans_id}) ? $len : $len{$trans_id};
}
close IN;
foreach my $trans_id(keys %len)
{
	$filtered_1{$trans_id}++ if $len{$trans_id} >= 50;
}

my $gff_file = "$new_dir/$species.newGene_final.gff";
open (IN, $gff_file) || die $!;
my (%info, %exon_num, %exon_parent);
while(<IN>)
{
	chomp;
	my @data = split /\t+/;
	if($data[2] eq 'gene')
	{
		$data[-1] =~ /ID=(\S+)$/;
		my $gene_id = $1;
		$info{$gene_id}{'gene'}{$gene_id} = "$_\n";
	}elsif($data[2] eq 'mRNA')
	{
		$data[-1] =~ /ID=(\S+);Parent=(\S+)$/;
		my ($trans_id, $gene_id) = ($1, $2);
		$exon_parent{$trans_id} = $gene_id;
		$info{$gene_id}{'mRNA'}{$trans_id} = "$_\n";
	}elsif($data[2] eq 'exon')
	{
		$data[-1] =~ /Parent=(\S+)$/;
		my $trans_id = $1;
		my $gene_id = $exon_parent{$trans_id};
		$info{$gene_id}{'mRNA'}{$trans_id} .= "$_\n";
		$exon_num{$trans_id}++;
	}
}
close IN;

my $filter_gff = "$new_dir/$species.newGene_final.filtered.gff";
open OUT, ">$filter_gff";
print OUT "#Seq_ID\tSource\tType\tStart\tEnd\tScore\tStrand\tPhase\tAttributes\n";	
foreach my $gene_id(sort keys %info)
{
	my $lines = $info{$gene_id}{'gene'}{$gene_id};
	my $flag;
	foreach my $trans_id(sort keys %{$info{$gene_id}{'mRNA'}})
	{
		next unless exists $filtered_1{$trans_id};
		my $exon_num = $exon_num{$trans_id};
		next if ($exon_num == 1 and !defined $singleexon);
		$filtered{$trans_id}++;
		$lines .= $info{$gene_id}{'mRNA'}{$trans_id};
		$flag++;
	}
	print OUT $lines if $flag;
}
close OUT;
open IN, "$new_dir/$species.newGene_final_tracking.list";
open OUT, ">$new_dir/$species.newGene_filtered_final_tracking.list";
while(<IN>)
{
	chomp;
	my @data = split /\t+/;
	next unless exists $filtered{$data[7]};
	print OUT "$_\n";
}
close IN;
close OUT;

system "rm -r $new_dir/tmp_orf";



sub USAGE {#
	my $usage=<<"USAGE";

Usage:
	-s	       index of species,
	-g 		genome,
	-newdir		final_track,
	-singleexon	singleexon,

	
USAGE
	print $usage;
	exit;
}

