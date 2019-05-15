#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $version="1.0.0";

#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($quality,$out);
GetOptions(
				"help|?" =>\&USAGE,
				"quality:s"=>\$quality,
				"prefix:s"=>\$out,
				) or &USAGE;
&USAGE unless ($quality  and $out);

my %CFG=%{readconf("$Bin/../../CFG")};
my $Rscript = $CFG{Rscript};

#$quality=&ABSOLUTE_DIR($quality);
#$out=&ABSOLUTE_DIR($out);

open IN, "<$quality";
open OUT,">$out.error.stat";
print OUT "#Sites\tError_ration\n";
while(<IN>){
	chomp;
	next if(/^$/ || /^\#/);
#	print $_;
	my @line = split /\t/,$_;
	my $site = $line[0];
	my $q = 0;
	for (my $i=1;$i<@line ;$i++) {
#		next if ($i==0);
		my $p = $line[$i]*($i-1)/100;
#		print "$line[$i]\n";
		$q += $p;
	}
	my $error = 10**(-$q/10)*100;
	print OUT "$site\t$error\n";
}

close IN;
close OUT;


`$Rscript $Bin/simpleBar.r --infile "$out.error.stat" --outfile "$out.error.png" --x.col 1 --y.col 2 --x.lab Position --y.lab Error_rate  --axis.size 8  --title.lab "$out Base Calling Error Rate per Cycle"`;

sub ABSOLUTE_DIR{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}



sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Usage:
  Options:
  -quality    <dir>    quality file from Data Assess;
  -prefix  <num>   Erorr Ration File prefix;
  -h         Help

USAGE
	print $usage;
	exit;
}
