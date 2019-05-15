#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fa,$gtf,$od);

GetOptions(
    "fa:s" =>\$fa,
    "gtf:s" =>\$gtf,
    "o:s"   =>\$od,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($fa and $gtf and $od);
#######################
my %gene_seq;
open IN,"$fa" || die $!;
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
my %trans_info;
#my %gene2trans;
open IN,"$gtf" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^FusionID/ || /^\#/);
	my @tmp = split /\t+/,$_;
	#$gene2trans{$tmp[0]}{$1} = 1;
        $trans_info{$tmp[0]}{gene1}{Start} = $tmp[7];
         $trans_info{$tmp[0]}{gene1}{End} = $tmp[8];
         $trans_info{$tmp[0]}{gene1}{Strand} = $tmp[5];
         $trans_info{$tmp[0]}{gene1}{Chr} = $tmp[6];
         $trans_info{$tmp[0]}{gene2}{Start} = $tmp[12];
         $trans_info{$tmp[0]}{gene2}{End} = $tmp[13];
         $trans_info{$tmp[0]}{gene2}{Strand} = $tmp[10];
         $trans_info{$tmp[0]}{gene2}{Chr} = $tmp[11];
	}

close IN;
open OUT,">$od" || die $!;
foreach my $gene (sort keys %trans_info) {
			my ($trans_len,$trans_seq);
			my $gen1_len = ($trans_info{$gene}{gene1}{End} - $trans_info{$gene}{gene1}{Start} + 1);
			my $gen2_len = ($trans_info{$gene}{gene2}{End} - $trans_info{$gene}{gene2}{Start} + 1);	
			my $trans_seq1 = substr ($gene_seq{$trans_info{$gene}{gene1}{Chr}}, $trans_info{$gene}{gene1}{Start} - 1,$gen1_len);
			if ($trans_info{$gene}{gene1}{Strand} eq "-") {
				$trans_seq1 =~tr/atgcATCGuU/tacgTAGCAA/;
				$trans_seq1 = reverse $trans_seq1;
			}
			my $trans_seq2 = substr ($gene_seq{$trans_info{$gene}{gene2}{Chr}}, $trans_info{$gene}{gene2}{Start} - 1,$gen2_len);
			if ($trans_info{$gene}{gene2}{Strand} eq "-") {
				$trans_seq2 =~tr/atgcATCGuU/tacgTAGCAA/;
				$trans_seq2 = reverse $trans_seq2;
			}
			$trans_seq=$trans_seq1.$trans_seq2;
			
			print OUT ">$gene\n$trans_seq\n";
}
close OUT;
###########
sub USAGE {
	my $usage=<<"USAGE";
 ProgramName:
     Contact:	Ma Liming <malm\@biomarker.com.cn> 
Program Date:	2016.08.29
      Modify:	
 Description:	This program is used to ......
       Usage:
		Options:
		
		-fa     <file>    the genome file 
	        -gtf    <file>    the gff file 
	        -o     <file>     output file of gene fasta 
		-h		help
Example:
            perl $Script --fa  genemo.fa --gtf genemo.gtf --o gene.fa

USAGE
	print $usage;
	exit;
}
