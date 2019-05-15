#!/usr/bin/perl -w
#Writer           songmm <songmm@biomarker.com.cn>
my $version="1.0.0";
my $BEGIN=time();

use strict;
use warnings;
use experimental 'smartmatch';
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd 'abs_path';
my $programe_dir=basename($0);
my $path=dirname($0);
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($circexplorer,$outfile);
GetOptions(
				"help|?" =>\&USAGE,
                "c:s"=>\$circexplorer,
				"o:s"=>\$outfile,
				) or &USAGE;
&USAGE unless ($circexplorer and $outfile) ;
open(IN, $circexplorer) or die $!;
open(OUT, ">$outfile") or die $!;
my %hash;
while (<IN>) {
    chomp;
    if(/^#/){print OUT "$_\n";next;}
    my @info = split /\t/,$_;
    my($chr,$start,$end)=@info[0..2];
    my $junction_reads = $info[12];
    if (exists $hash{"$chr;$start;$end"}) {
        my $current_junction_reads = (split /\t/,$hash{"$chr;$start;$end"},2)[0];
        if ($current_junction_reads < $junction_reads) {
            $hash{"$chr;$start;$end"}=$junction_reads."\t".$_;
        }
    }
    else
    {
        $hash{"$chr;$start;$end"}=$junction_reads."\t".$_;
    }
}
close IN;


foreach my $key(keys %hash)
{
    my $line = (split /\t/,$hash{$key},2)[1];
    print OUT "$line\n";
}
close OUT;



sub USAGE
{
    	my $usage=<<"USAGE";
Program: filter the same circRNA by junctions reads;
Version: $version
Contact: songmm <songmm\@biomarker.com.cn>

Description:

Usage:
  -c    <file>  circRNA file must be given;
  -o   <dir>    output file must be given;
USAGE
	print $usage;
	exit;
}
