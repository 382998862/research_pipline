#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($enrich_file,$key,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"enrich_file:s"=>\$enrich_file,
				"key:s"=>\$key,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($enrich_file and $od and $key) ;
################################
$enrich_file=abs_path($enrich_file);
&MKDIR($od);
$od=abs_path($od);

####################################### Make info for the program
my %enrich_info;
&info_stat ($enrich_file,\%enrich_info,$key);

my $scatter = "$od/$key.KEGG.Phase.png";
&Phase_plot ("$od/$key.KEGG.list",$scatter);

######################################### subs
sub info_stat
{
	my $xls = shift;
	my $enrich_info = shift;
	my $key = shift;

	open IN,"$xls" || die $!;
	open OUT,">$od/$key.KEGG.list" || die $!;
	print OUT "#pathway\tko\tenrichment_factor\tcorrect_p\n";
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		my $info = $_;
		my @tmp = split /\t+/,$_;
		my ($m,$n,$M,$N);
		$tmp[2] =~ /(\d+)\sout\sof\s(\d+)/;
		$m = $1;
		$n = $2;
		$tmp[3] =~ /(\d+)\sout\sof\s(\d+)/;
		$M = $1;
		$N = $2;

#		$enrich_info->{$tmp[0]}->{$tmp[1]}{"Rich_f"} = $tmp[-1];
#		$enrich_info->{$tmp[0]}->{$tmp[1]}{"P_value"} = $tmp[-2];
		print OUT "$tmp[0]\t$tmp[1]\t";
		my $rich_factor=($m/$n)/($M/$N);
		print OUT sprintf("%.2f",$rich_factor),"\t$tmp[-1]\n";
#		printf OUT "%.2f",($n/$N)/($m/$M);
#		print OUT "\n"
#		print OUT "\t$tmp[-1]\n";
	}
	close IN;
	close OUT;
}

sub Phase_plot
{
	my $file = shift;
	my $png = shift;
    my $Rscript = "/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript";

	`$Rscript  $Bin/kegg_pathway_enrich_plot.r $file $png `;
}

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
Program: 
Version: $version
Contact: Meng Fei <mengf\@biomarker.com.cn>

Description:
	
Usage:
	-enrich_file      pathway enrichment file              must be given
	-od               output dir                           must be given
	-key              out file prefix                      must be given
USAGE
	print $usage;
	exit;
}
