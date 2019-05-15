#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";

#######################################################################################
# ------------------------------------------------------------------
# GetOptions   
# ------------------------------------------------------------------
my ($align,$key,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"align:s"=>\$align,
				"key:s"=>\$key,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($align and $key);
################################
$align=&ABSOLUTE_DIR($align);
$od=$od || "./";
` mkdir $od ` if (!-d "$od");
$od=&ABSOLUTE_DIR($od);

####################################### align process
open IN,"$align" || die $!;
open OUT,">$od/$key.bowtie_silver.list" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/);
	my @tmp=split /\s+/,$_;
	my $index="rRNA";
	print OUT "$tmp[0]\t$index\t$tmp[1]\t$tmp[2]\n";
}
close OUT;


######################################### subs
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
	-align               Uniq fasta file       must be given
	-key                 sample name           must be given
	-od                  output dir            choice default "./";
USAGE
	print $usage;
	exit;
}
