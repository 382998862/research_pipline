#!/usr/bin/perl -w
use Getopt::Std;
use Getopt::Long;
use List::Util qw(sum);
use Cwd qw(abs_path);

my ($known_count,$known_fpkm,$known_gff,$known_gtf);
GetOptions(
	"known_lncrna_count:s" =>\$known_count,
	"known_lncrna_fpkm:s" =>\$known_fpkm,
	"known_lncrna_gff:s" =>\$known_gff,
	"known_lncrna_gtf:s" =>\$known_gtf,
	"help|h" =>\&USAGE,
) or &USAGE;

&USAGE unless (defined $known_count && defined $known_fpkm && defined $known_gff && defined $known_gtf);
$known_count=abs_path($known_count);
$known_fpkm=abs_path($known_fpkm);
$known_gff=abs_path($known_gff);
$known_gtf=abs_path($known_gtf);

`mv $known_count $known_count.tmp`;
print("mv $known_count $known_count.tmp\n");
`mv $known_fpkm $known_fpkm.tmp`;
print("mv $known_fpkm $known_fpkm.tmp\n");
`mv $known_gff $known_gff.tmp`;
print("mv $known_gff $known_gff.tmp\n");
`mv $known_gtf $known_gtf.tmp`;
print("mv $known_gtf $known_gtf.tmp\n");

my (@zero,%h);
open(IN,"$known_count.tmp") or die $!;
open(OUT,">$known_count");
while(<IN>){
	chomp;
	if(/^#/){print OUT "$_\n";next;}
	my ($id,$count)=split /\t/,$_,2;
	my @A = split /\t/,$count;
	my $sum = sum @A;
	if($sum ==0){
		push @zero,$id;
		next;
	}else{
		$h{$id}=1;
		print OUT "$_\n";
	}
}
close(IN);
close(OUT);
print ("creat filtered known_lncRNA.count file finish! \n");

open(IN,"$known_fpkm.tmp");
open(OUT,">$known_fpkm");
while(<IN>){
	chomp;
	if(/^#/){print OUT "$_\n";next;}
	my($id,$fpkm)=split /\t/,$_,2;
	if(exists $h{$id}){
		print OUT "$_\n";
	}
}
close(IN);
close(OUT);
print ("creat filtered known_lncRNA.fpkm file finish! \n");
my(%h1,%h2,@gene);
open(IN,"$known_gff.tmp");
open(OUT,">$known_gff");
while(<IN>){
	chomp;
	next if(/^#|^\s$/);
	my($chr,$from,$type,$start,$end,$dot,$strand,$score,$att)=split;
	if($type=~/lncRNA_gene|processed_transcript|lincRNA_gene|gene/){
		if($att=~/ID=(.*?);/){
			$h1{$1}{gene}=$_;
		}
	}
	if($type=~/lncRNA|lincRNA|transcript/){
		if($att=~/ID=(.*?);Parent=(.*?);/){
			if(exists $h{$1}){
				push @gene,$2;
				$h2{$2}=1;
				$h1{$2}{rna}{$1}=$_;
			}
		}
	}
	if($type eq "exon"){
		if($att=~/ID=(.*?);Parent=(.*?);/){
			$h1{$2}{exon}{$1}=$_;
		}
	}
}

foreach my $g (keys %h1){
	if(exists $h2{$g}){
		print OUT "$h1{$g}{gene}\n";
		foreach my $rna(keys %{$h1{$g}{rna}}){
				print OUT "$h1{$g}{rna}{$rna}\n";
				foreach my $exon(keys %{$h1{$rna}{exon}}){
					print OUT "$h1{$rna}{exon}{$exon}\n";
				}
		}
	}
}
close(IN);
close(OUT);

open(IN,"$known_gtf.tmp");
open(OUT,">$known_gtf");
while(<IN>){
	chomp;
	my($chr,$from,$type,$start,$end,$dot,$strand,$score,$att)=split /\t/,$_;
	if($att=~/transcript_id "(.*?)"; gene_id "(.*?)";/){
		if(exists $h{$1}){
			print OUT "$_\n";
		}
	}
}
close(IN);
close(OUT);

############################################################################
sub USAGE {
        my $usage=<<"USAGE";
--------------------------------
Program:
Version:
Usage:
	-known_lncrna_count
	-known_lncrna_fpkm
	-known_lncrna_gff
	-known_lncrna_gtf
	-h help documents

Example: perl filter_fpkm.pl -known_lncrna_count lncRNA_count.list -known_lncrna_fpkm lncRNA_fpkm.list -known_lncrna_gff Homo_sapiens.GRCh38.lncRNA.85.gff3 -known_lncrna_gtf Known_lncRNA.gtf
--------------------------------
USAGE
        print $usage;
        exit;
}

