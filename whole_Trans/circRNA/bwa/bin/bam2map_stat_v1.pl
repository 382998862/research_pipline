#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.1";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($Index,$bam,$od,$totalRead);
GetOptions(
				"help|?"      =>\&USAGE,
				"i:s"         =>\$Index,
				"totalRead:s" =>\$totalRead,
				"bam:s"       =>\$bam,
				"od:s"        =>\$od,
				) or &USAGE;
&USAGE unless ($Index and $bam and $od and $totalRead) ;

$bam = &ABSOLUTE_DIR($bam);       &MKDIR($od);
$od  = &ABSOLUTE_DIR($od);
print STDERR "bam file: $bam\n";
print STDERR "Out file: $Index.mapped.stat.xls && $Index.insertSize.png\n";
###############################
my %map_stat;
my %insertsize;
my %multiple_ID;
my %multiple_dis;

open IN, "samtools view $bam |" || die;
while (<IN>) {
	chomp;
	my @tmp = split /\s+/,$_;
	my ($read, $indel, $type, $insert) = @tmp[0, 5, 6, 8];
	my ($mismatch, $hit_num);
	for (9..$#tmp) {
		if ($tmp[$_] =~ /NM:i:(\d+)/) {
			$mismatch = $1;
		}
		if ($tmp[$_] =~ /NH:i:(\d+)/) {
			$hit_num = $1;
		}
	}
	if ($hit_num == 1) {
		$map_stat{uniq}++;
		if ($indel !~ /N/ && $mismatch == 0) {
			$map_stat{perfect}++;
			$map_stat{mis}{$mismatch}++;
		}
		if ($indel =~ /N/ && $mismatch == 0) {
			$map_stat{Indel}++;
			$map_stat{mis}{$mismatch}++;
		}
		if ($indel!~/N/ && $mismatch!=0) {
			if ($mismatch<=2) {
				$map_stat{mis}{$mismatch}++;
			}
			else {
				$map_stat{mis}{3}++;
			}
		}
		if ($indel=~/N/ && $mismatch!=0) {
			$map_stat{misIndel}++;
		}
		if ($type=~/=/) {
			$map_stat{pair_end}++;
			$insertsize{abs($insert)}++;
		}
		if ($type=~/\*/) {
			$map_stat{single_end}++;
		}
	}
	if ($hit_num!=1 && !defined $multiple_ID{$read}) {
		$multiple_ID{$read} = $hit_num;
		$map_stat{multi}++;
		if ($hit_num<=10) {
			$multiple_dis{low}++;
		}
		if ($hit_num>=10) {
			$multiple_dis{high}++;
		}
		if ($indel!~/N/ && $mismatch==0) {
			$map_stat{'perfect'}++;
			$map_stat{'mis'}{$mismatch}++;
		}
		if ($indel=~/N/ && $mismatch==0) {
			$map_stat{'Indel'}++;
			$map_stat{'mis'}{$mismatch}++;
		}
		if ($indel!~/N/ && $mismatch!=0) {
			if ($mismatch<=2) {
				$map_stat{'mis'}{$mismatch}++;
			}
			else {
				$map_stat{'mis'}{3}++;
			}
		}
		if ($indel=~/N/ && $mismatch!=0) {
			$map_stat{'misIndel'}++;
		}
		if ($type=~/=/) {
			$map_stat{'pair_end'}++;
			$insertsize{abs($insert)}++;
		}
		if ($type=~/\*/) {
			$map_stat{'single_end'}++;
		}
	}
}
my $total_map = $map_stat{'uniq'} + $map_stat{'multi'};
my $map_ratio = ($total_map/$totalRead) * 100;
my $perfect_ratio = $map_stat{'perfect'}/$total_map * 100;


my $mappedStat_result = "$od/$Index.mappedStat.xls";
open (OUT, ">", $mappedStat_result) || die;
print  OUT "Total Reads\t$totalRead\t100%\n";
printf OUT "%s\t%s\t%.2f", "mapped Reads", $total_map, $map_ratio;
print  OUT "%\n";
printf OUT "%s\t%s\t%.2f", "Perfect Map", $map_stat{perfect}, $perfect_ratio;
print  OUT "%\n\n";
print  OUT "mismatch\n";
foreach (sort {$a<=>$b} keys %{$map_stat{mis}}) {
	my $mis_ratio;
	if ($_ <= 2) {
		print OUT "$_\t";
		print OUT "$map_stat{mis}{$_}\t";
		$mis_ratio = $map_stat{mis}{$_}/$total_map*100;
		printf OUT "%.2f",$mis_ratio;
		print OUT "%\n";
	}
	else {
		print OUT ">=3\t";
		print OUT "$map_stat{mis}{$_}\t";
		$mis_ratio = $map_stat{mis}{$_}/$total_map*100;
		printf OUT "%.2f",$mis_ratio;
		print OUT "%\n";
	}
}
print OUT "Indel\n";
my $Indel_ratio = $map_stat{Indel}/$total_map*100;
printf OUT "%s\t%s\t%.2f","1", $map_stat{Indel}, $Indel_ratio;
print  OUT "%\n";

my $mis_Indel_ratio = $map_stat{misIndel}/$total_map*100;
printf OUT  "%s\t%s\t%.2f", "mismatch+Indel", $map_stat{misIndel}, $mis_Indel_ratio;
print  OUT "%\n\n";

my $uniq_ratio = $map_stat{uniq}/$total_map*100;
printf OUT "%s\t%s\t%.2f","Uniq Map", $map_stat{uniq}, $uniq_ratio;
print  OUT "%\n";
my $multi_ratio = $map_stat{multi}/$total_map*100;
printf OUT  "%s\t%s\t%.2f", "Multiple Map", $map_stat{multi}, $multi_ratio;
print  OUT "%\n\n";

my $pair_map_ratio = $map_stat{pair_end}/$total_map*100;
printf OUT "%s\t%s\t%.2f", "Pair Map", $map_stat{pair_end}, $pair_map_ratio;
print  OUT "%\n";
my $single_map_ratio = $map_stat{single_end}/$total_map*100;
printf OUT "%s\t%s\t%.2f", "Single Map", $map_stat{single_end}, $single_map_ratio;
print  OUT "%\n";

close OUT;

#print Dumper %multiple_dis;

my $max = (sort {$b <=> $a}values %insertsize)[0];
my $insertSize_list = "$od/$Index.insertSize.list";
open (OUT , ">", $insertSize_list) or die $!;
print OUT "Type:Line","\n";
print OUT "Width:600","\n";
print OUT "Height:400","\n";
print OUT "WholeScale:0.9","\n";
print OUT "FontSize:25","\n";
print OUT "X:InsertSize","\n";
print OUT "Y:Number of Reads","\n";
print OUT "XStep:100","\n";
print OUT "YStep:",int($max/10),"\n";
print OUT "XStart:0","\n";
print OUT "YStart:0","\n";
print OUT "XEnd:800","\n";
print OUT "YEnd:",$max+200,"\n\n";
print OUT "Color:red","\n";

foreach my $insert (sort {$a <=> $b} keys %insertsize) {
	next if ($insert<=100);
	print OUT $insert,":",$insertsize{$insert},"\n";
	last if ($insert>=800);
}
close (OUT);

chdir($od);
system ("/share/nas2/genome/biosoft/distributing_svg_4.74/distributing_svg.pl $Index.insertSize.list $Index.insertSize.svg");
system ("/share/nas2/genome/biosoft/distributing_svg_4.74/svg2xxx_release/svg2xxx $Index.insertSize.svg");

my $Rscript = "/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript";
system "grep \"^[0-9]\" $od/$Index.insertSize.list |sed 's/:/\t/' > $od/$Index.insertSize.r.list";
system "$Rscript $Bin/pointOrLine.r --infile $od/$Index.insertSize.r.list --outfile $od/$Index.insertSize.r.png --x.col 1 --y.col 2 --x.lab \"Insert Size (bp)\" --y.lab \"Number of Reads\" --is.line --line.color 1 ";


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
Contact: mengf <mengf\@biomarker.com.cn>

Usage:
  -i               <str>    Index of inFiles and outFiles. Index.mappedStat.xls, Index.insertSize.png,
  -totalRead       <int>    Total Read Number;
  -bam             <str>    Infile;
  -od              <str>    OUT file DIR;

USAGE
	print $usage;
	exit;
}
