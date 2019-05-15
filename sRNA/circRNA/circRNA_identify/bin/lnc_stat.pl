#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my %config=%{readconf("$Bin/../../project.cfg")};
my ($lengthFile,$odir,$key);
GetOptions(
	"help|?"=>\&USAGE,
	"i:s"   =>\$lengthFile,
	"k:s"=>\$key,
	"o:s"   =>\$odir
	)or &USAGE;
&USAGE unless ($lengthFile and $odir and $key);
my %len;
open IN ,"$lengthFile";
while (<IN>){
	next if($.==1);
	my @tmp= (split /\s+/,$_);
	if ($tmp[1]>=0 && $tmp[1]<200){$len{$tmp[0]}{200}++;}	
	elsif ($tmp[1]>=200 && $tmp[1]<400){$len{$tmp[0]}{400}++;}
	elsif ($tmp[1]>=400 && $tmp[1]<600){$len{$tmp[0]}{600}++;}	
	elsif ($tmp[1]>=600 && $tmp[1]<800){$len{$tmp[0]}{800}++;}
	elsif ($tmp[1]>=800 && $tmp[1]<1000){$len{$tmp[0]}{1000}++;}	
	elsif ($tmp[1]>=1000 && $tmp[1]<1200){$len{$tmp[0]}{1200}++;}
	elsif ($tmp[1]>=1200 && $tmp[1]<1400){$len{$tmp[0]}{1400}++;}	
	elsif ($tmp[1]>=1400 && $tmp[1]<1600){$len{$tmp[0]}{1600}++;}
	elsif ($tmp[1]>=1600 && $tmp[1]<1800){$len{$tmp[0]}{1800}++;}	
	elsif ($tmp[1]>=1800 && $tmp[1]<2000){$len{$tmp[0]}{2000}++;}
	elsif ($tmp[1]>=2000 && $tmp[1]<2200){$len{$tmp[0]}{2200}++;}	
	elsif ($tmp[1]>=2200 && $tmp[1]<2400){$len{$tmp[0]}{2400}++;}
	elsif ($tmp[1]>=2400 && $tmp[1]<2600){$len{$tmp[0]}{2600}++;}	
	elsif ($tmp[1]>=2600 && $tmp[1]<2800){$len{$tmp[0]}{2800}++;}
	elsif ($tmp[1]>=2800 && $tmp[1]<3000){$len{$tmp[0]}{3000}++;}	
	elsif ($tmp[1]>=3000){$len{$tmp[0]}{'3000+'}++;}

}
my $all_len;
close IN;
open OUT ,">$odir/$key.stat.len";
foreach my $type(sort {$a cmp $b} keys%len)
{
	for(my $i=200;$i<=3000;$i+=200)
	{
		if(exists $len{$type}{$i})
		{
			print OUT "$type\t$i\t$len{$type}{$i}\n";
		}
	}
	if(exists $len{$type}{'3000+'})
	{
		print OUT "$type\t3000+\t$len{$type}{'3000+'}\n";
	}

}
close OUT;
#plot
`$config{Rscript} $Bin/R/lengthPlot.r --infile $odir/$key.stat.len --outfile $odir/$key.circRNA.length.distribution  --x.col 2 --y.col 3 --skip 0 --x.lab \"Length Range\" --y.lab \"The number of circRNA\" --title.lab \"CircRNA Length Distribution of $key\" --axis.size 5`;
#####################################################################################################################################################################################################################
sub max{
	my $a=shift;
	my $b=shift;
	my $tmp;
	if ($a >= $b){
		$tmp=$a;
	}else{
		$tmp=$b;
	}
	return $tmp;
}


sub USAGE {
	        my $usage=<<"USAGE";
----------------------------------------------------------------------------------------
   Program: 	$Script
   Discription: Extract len and exon_number from filter_final.gtf & merged.gtf\n\n
   USAGE:
   		-i	input file
		-o	output dir directory
		-k 	sample id
		-h	help documents
----------------------------------------------------------------------------------------
USAGE
			print $usage;
			exit ;
		}
