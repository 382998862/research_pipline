#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my($dir, $out,$type);

GetOptions (
	"help|?"=>\&USAGE,
	"dir:s" =>\$dir,
	"o:s"=>\$out,
	"type:s"=>\$type
)or &USAGE;
&USAGE unless($dir and $out);
my $fastq_dump = "/share/nas34/songmm/cirRNAtest/tool/sratoolkit.2.5.2-centos_linux64/bin/fastq-dump";
my @file = glob("$dir/*.sra");
$type ||= "fq";
foreach(@file)
{
	if($type =~/fa|fasta/)
	{
		`$fastq_dump --fasta 0 --split-files $_  -O $out && touch $out/done.ok`;
	}
	else
	{
		`$fastq_dump --split-files $_  -O $out && touch $out/done.ok`;
	}
}


sub USAGE
{
	my $usage =<<"USAGE";
	ProgramName:
     Contact:   Simon Young <yangxh\@biomarker.com.cn> 
Program Date:   2012.07.02
      Modify:   
 Description:   This program is used to ......
       Usage:
        Options:
        -dir <file>   input dir,xxx format,forced
        -o <file>   output dir,forced
        -type <file>   file type,fq or fa,forced
        -h      help
USAGE
    print $usage;
    exit;

}