#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
#my %config=%{readconf("$Bin/../../../../../config/db_file.cfg")}; 
my %config=%{readconf("$Bin/../../../../Config/lncRNA_pip.cfg")};
my $BEGIN_TIME=time();
my $version="2.0.0";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($Index,$bam,$od,$total_num);
GetOptions(
				"help|?"      =>\&USAGE,
				"i:s"         =>\$Index,
				"totalRead:s" =>\$total_num,
				"bam:s"       =>\$bam,
				"od:s"        =>\$od,
				) or &USAGE;
&USAGE unless ($Index and $bam and $od and $total_num) ;

$bam = &ABSOLUTE_DIR($bam);       &MKDIR($od);
$od  = &ABSOLUTE_DIR($od);
print STDERR "bam file: $bam\n";
print STDERR "Out file: $Index.mapped.stat.xls && $Index.insertSize.png\n";
###############################
my (%left_mapped_reads, %left_uniq_mapped_reads);
my (%right_mapped_reads, %right_uniq_mapped_reads);
my %perfect_mapped_reads;
my %pair_mapped_reads;
my %insertsize;

system "$config{samtools} view $bam > $bam.txt";
open IN, "$bam.txt" or die;
 #bedtools bamtobed -i 
while (<IN>) {
    chomp;
    my ($qid, $flag, $ref_id, $pos, $mapQ, $CIGAR, $rnext, $pnext, $tlen, $seq, $qual, $more) = (split /\t/,$_,12);
    next if ($ref_id eq '*');

    if (64 & $flag) {
        $left_mapped_reads{$qid} = 1;
        $left_uniq_mapped_reads{$qid} = 1 if ($mapQ >=40);
    } else {
        $right_mapped_reads{$qid} = 1;
        $right_uniq_mapped_reads{$qid} = 1 if ($mapQ >=40);
    }

    if (1 & $flag) {
        unless (12 & $flag) {
            $pair_mapped_reads{$qid} = 1;
            $insertsize{abs($tlen)}++ if ($rnext eq '=');
        }
    }
}

close IN;


my %map_strand;
`$config{bedtools} bamtobed -i $bam > $bam.bed`;
open IN, "$bam.bed" or die;
 #bedtools bamtobed -i 
while (<IN>) {
    chomp;
    my ($ref_id, $start,$end,$qid,undef, $str) = split (/\t/,$_);
	$map_strand{'plus'}{$qid}=1  if ($str eq '+');
	$map_strand{'minus'}{$qid}=1 if ($str eq '-');
}
close IN;

foreach my $keys (keys %{$map_strand{'plus'}}) {
	if (exists $map_strand{'plus'}{$keys} and exists $map_strand{'minus'}{$keys}) {
		delete $map_strand{'plus'}{$keys};
		delete $map_strand{'minus'}{$keys};
	}
}

my $left_mapped_num = keys %left_mapped_reads;
my $right_mapped_num = keys %right_mapped_reads;
my $left_uniq_mapped_num = keys %left_uniq_mapped_reads;
my $right_uniq_mapped_num = keys %right_uniq_mapped_reads;
my $plus_strand_num=keys %{$map_strand{'plus'}};
my $minus_strand_num=keys %{$map_strand{'minus'}};


my $mapped_num = $left_mapped_num + $right_mapped_num;
my $mapped_per = sprintf("%.2f", $mapped_num/$total_num*100).'%';

my $uniq_mapped_num = $left_uniq_mapped_num + $right_uniq_mapped_num;
my $uniq_mapped_per = sprintf("%.2f", $uniq_mapped_num/$total_num*100).'%';
my $mult_mapped_num = $mapped_num - $uniq_mapped_num;
my $mult_mapped_per = sprintf("%.2f", $mult_mapped_num/$total_num*100).'%';

my $pair_mapped_num = (keys %pair_mapped_reads) * 2;
my $pair_mapped_per = sprintf("%.2f", $pair_mapped_num/$total_num*100).'%';
my $sngl_mapped_num = $mapped_num - $pair_mapped_num;
my $sngl_mapped_per = sprintf("%.2f", $sngl_mapped_num/$total_num*100).'%';

my $plus_strand_per=sprintf("%.2f", $plus_strand_num/$total_num*100).'%';
my $minus_strand_per=sprintf("%.2f", $minus_strand_num/$total_num*100).'%';

my $print = <<"_STAT_";
Total Reads\t$total_num\t100\%
mapped Reads\t$mapped_num\t$mapped_per

Uniq Map\t$uniq_mapped_num\t$uniq_mapped_per
Multiple Map\t$mult_mapped_num\t$mult_mapped_per

Pair Map\t$pair_mapped_num\t$pair_mapped_per
Single Map\t$sngl_mapped_num\t$sngl_mapped_per

Only Map Plus Strand\t$plus_strand_num\t$plus_strand_per
Only Map Minus Strand\t$minus_strand_num\t$minus_strand_per


_STAT_

my $mappedStat_result = "$od/$Index.mappedStat.xls";
open (OUT, ">", $mappedStat_result) || die;
print OUT $print;
close OUT;

#print Dumper %multiple_dis;

my $max = (sort {$b <=> $a}values %insertsize)[0];
my $insertSize_list = "$od/$Index.insertSize.list";
open (OUT , ">", $insertSize_list) or die $!;
#print OUT "Type:Line","\n";
#print OUT "Width:600","\n";
#print OUT "Height:400","\n";
#print OUT "WholeScale:0.9","\n";
#print OUT "FontSize:25","\n";
#print OUT "X:InsertSize","\n";
#print OUT "Y:Number of Reads","\n";
#print OUT "XStep:100","\n";
#print OUT "YStep:",int($max/10),"\n";
#print OUT "XStart:0","\n";
#print OUT "YStart:0","\n";
#print OUT "XEnd:800","\n";
#print OUT "YEnd:",$max+200,"\n\n";
#print OUT "Color:red","\n";

foreach my $insert (sort {$a <=> $b} keys %insertsize) {
	next if ($insert<=100);
	print OUT $insert,":",$insertsize{$insert},"\n";
	last if ($insert>=800);
}
close (OUT);

chdir($od);
#system ("/share/nas2/genome/biosoft/distributing_svg_4.74/distributing_svg.pl $Index.insertSize.list $Index.insertSize.svg");
#system ("/share/nas2/genome/biosoft/distributing_svg_4.74/svg2xxx_release/svg2xxx $Index.insertSize.svg");

my $Rscript = $config{Rscript};
system "grep \"^[0-9]\" $od/$Index.insertSize.list |sed 's/:/\t/' > $od/$Index.insertSize.r.list";
system "$Rscript $Bin/pointOrLine.r --infile $od/$Index.insertSize.r.list --outfile $od/$Index.insertSize.png --x.col 1 --y.col 2 --x.lab \"Insert Size (bp)\" --y.lab \"Number of Reads\" --is.line --line.color 1 ";
#system "$config{Rscript} $Bin/plot_insertsize.R infile=$od/$Index.insertSize.r.list outfile=$od/$Index.insertSize.r.png bg=F";
system "rm $bam.bed";
#########################subs
sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;chomp $return;
	}
	else
	{
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}
sub USAGE {#
	my $usage=<<"USAGE";
Program: Stat map_ratio of PairEND Read and draw insertsize png (.bam file format);
Version: $version
Contact: Simon Young <yangxh\@biomarker.com.cn>
   Data: 2015-01-30

Usage:
  -i               <str>    Index of inFiles and outFiles. Index.mappedStat.xls, Index.insertSize.png,
  -totalRead       <int>    Total Read Number;
  -bam             <str>    Input file, accepted_hits.bam;
  -od              <str>    OUT file DIR;

USAGE
	print $usage;
	exit;
}
