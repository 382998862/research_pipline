#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($in, $out,$sample);

GetOptions(
    "i:s" =>\$in,
    "o:s" =>\$out,
    "s:s"=>\$sample,
    "help|h" =>\&USAGE
    ) or &USAGE;
&USAGE unless ($in and $out and $sample);
my @files = split /,/,$in;
my @samples = split/,/,$sample;
my %hash;
for(my $i = 0;$i<@files;$i++)
{
	open IN,"$files[$i]" or die $!;
	<IN>;
	while(<IN>)
	{
    		chomp;
    		print OUT "$_\n",next if/#/;
    		my @line = split/\t/,$_;
		$hash{$line[4]}{$samples[$i]}++;
	}
	close IN;
}
open OUT,">$out" or die $!;
foreach my $type(sort keys%hash)
{
    foreach my $s(sort keys% {$hash{$type}})
    {
        print OUT "$s\t$type\t$hash{$type}{$s}\n";
    }
}
 close OUT;


#############################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------------------------------
   Program: $Script
   Contact: songmm <songmm\@biomarker.com.cn>
      Date: 2016-07-20

     Usage:
            --i     <FILE>      input file split by ","
            --o        <DIR>   output file
            --s         <STR>    T01,T02
            --h                 help documents

   Example:
            perl $Script --i T01.circ.intersect,T02.circ.intersect --o circRNA.type.stat --s T01,T02

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}