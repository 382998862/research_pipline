#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="2.0.0";
use newPerlBase;
my %config=%{readconf("$Bin/../../db_file.cfg")};

#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut);
GetOptions(
				"help|?"  => \&USAGE,
				"i:s" => \$fIn,
				"o:s" => \$fOut,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);

#my $ConvertSVGtoPNG="/share/nas1/litc/tools/svg2xxx_release/svg2xxx";
my $ConvertSVGtoPNG=$config{svg2xxx1}; # 2014-12-11 ~ 

# ------------------------------------------------------------------
# Class Anno
# ------------------------------------------------------------------
my %Class=(
	"J" => [1,"Translation, ribosomal structure and biogenesis"],
	"A" => [2,"RNA processing and modification"],
	"K" => [3,"Transcription"],
	"L" => [4,"Replication, recombination and repair"],
	"B" => [5,"Chromatin structure and dynamics"],
	"D" => [6,"Cell cycle control, cell division, chromosome partitioning"],
	"Y" => [7,"Nuclear structure"],
	"V" => [8,"Defense mechanisms"],
	"T" => [9,"Signal transduction mechanisms"],
	"M" => [10,"Cell wall/membrane/envelope biogenesis"],
	"N" => [11,"Cell motility"],
	"Z" => [12,"Cytoskeleton"],
	"W" => [13,"Extracellular structures"],
	"U" => [14,"Intracellular trafficking, secretion, and vesicular transport"],
	"O" => [15,"Posttranslational modification, protein turnover, chaperones"],
	"C" => [16,"Energy production and conversion"],
	"G" => [17,"Carbohydrate transport and metabolism"],
	"E" => [18,"Amino acid transport and metabolism"],
	"F" => [19,"Nucleotide transport and metabolism"],
	"H" => [20,"Coenzyme transport and metabolism"],
	"I" => [21,"Lipid transport and metabolism"],
	"P" => [22,"Inorganic ion transport and metabolism"],
	"Q" => [23,"Secondary metabolites biosynthesis, transport and catabolism"],
	"R" => [24,"General function prediction only"],
	"S" => [25,"Function unknown"],
);
my $TILE="eggNOG Function Classification of Consensus Sequence";
my $Xtile="Function Class";
my $Ytile="Frequency";

# ------------------------------------------------------------------
# Get Data
# ------------------------------------------------------------------
my %Data;
&LoadData($fIn,\%Data);
open (OUT,">$fOut.stat") || die "can't open file  $fOut.stat\n ";
&OutHash(\%Data);
close(OUT);
`$config{Rscript}  $Bin/eggNOG_anno_plot.r  $fOut.stat  $fOut.png >eggNOG.classfy.plot.log 2>&1`;


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub LoadData {#
	my ($fIn,$Data)=@_;
    my %tmp;
	open (IN,"<",$fIn) or die $!;
	<IN>;
	while (<IN>) {
		chomp;
        next if /^\#/;
        my ($gi,undef,undef,undef,$fun)=split(/\t/);
        next if exists $tmp{$gi};
        $tmp{$gi}=$fun;
		foreach my $elem (split "",$fun) {
			$$Data{$elem}++;
		}
	}
	close (IN) ;
}
sub OutHash
{
	my ($hash)=@_;
	print OUT "#ID\tClass_Name\tNumbers\n";
	foreach my $sty (sort {$Class{$a}[0] <=> $Class{$b}[0]} keys %Class)
	{
		print OUT "$sty\t"."$Class{$sty}[1]\t".(${$hash}{$sty}||0)."\n";
	}
}
sub USAGE {#
	my $usage=<<"USAGE";
Program: eggNOGFunClassDrawer.pl
Version: $version
Contact: Li Tiancheng <litc\@biomarker.com.cn> <ltc_gs\@qq.com>
Revised: Yangsh <2011-07-20>
Revised: Wangyj <2016-03-08>
Description:
	Make eggNOG Annotation Classify Graph

Usage:
  Options:
  -i <file>  input eggNOG annotation result file, forced
  -o <file>  png output file.


USAGE
	print $usage;
	exit;
}
