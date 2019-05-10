#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename;
use Cwd 'abs_path';

my ($infile, $outfile);
GetOptions(
				"help|?" =>\&USAGE,
                "i:s"=>\$infile,
                "o:s"=>\$outfile,
				) or &USAGE;
&USAGE unless ($infile and $outfile);



my @file=split/,/,$infile;
my %hash;
my %region;
my %all_sample;
foreach my $anno(@file){
    my $basename=basename $anno;
    my ($sample)=$basename=~/^(.+?)\./;
    $all_sample{$sample}=1;
    open(IN,"$anno")||die $!;
    while (<IN>) {
        chomp;
        next if /^#/;
        my @line=split/\t+/,$_;
        if ($line[7]=~/SNPEFF_EFFECT\=(.+?)\;/) {
            $hash{$1}{$sample}++;
            $region{$line[0]}{$line[1]}=$1;
        }
        else{
            $hash{other}{$sample}++;
            $region{$line[0]}{$line[1]}="other";
        }
    }
    close IN;
}

my @sample=sort keys %all_sample;
my $sample=join("\t",@sample);
my %final_type;
foreach my $chr(keys %region){
    foreach my $pos(keys %{$region{$chr}}){
        my $type=$region{$chr}{$pos};
        $final_type{$type}++;
    }
}
open(OUT,">$outfile")||die $!;

print OUT "#Region\tType\t$sample\tall\n";

foreach my $type(sort keys %hash){
    print OUT "--\t$type\t";
    #foreach my $sample (sort keys %{$hash{$type}}){
    foreach my $sample(@sample){
        $hash{$type}{$sample}||=0;
        print OUT "$hash{$type}{$sample}\t";
    }
    print OUT "$final_type{$type}\n";
}
close OUT;


sub USAGE {#
	my $usage=<<"USAGE";
Discription:
		Produce snp/indel anno stat file
        Author:Liuxs(liuxs\@biomarker.com.cn)
        Data:2016/09/06
Usage:
  -i			<file>	input file ,split by comma ,eg T01.snp.anno.vcf,T02.snp.anno.vcf
  
  -o			<file>	output file

  -h		Help

USAGE
	print $usage;
	exit;
}
