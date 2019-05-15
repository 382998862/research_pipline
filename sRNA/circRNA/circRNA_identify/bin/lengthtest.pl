#!/usr/bin/perl
use strict;
use warnings;
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME = time();
my $version    = "1.0.0";
my ( $fIn, $out);
## input is CIRI output file;
GetOptions(
	"help|?"   => \&USAGE,
	"o:s"    => \$out,
	"i:s"      => \$fIn
) or &USAGE;
&USAGE unless ( $fIn and $out);
open( IN, $fIn ) or die $!;
open( OUT,  ">$out" )       or die $!;
print OUT "range\tcount\n";

my %hash;
while (<IN>) {
	next if (/\#/);
	my @line        = split /\t/, $_;
	my $length      = $line[3] - $line[2] + 1;
	my $type = $line[8];
	for(my $i=0;;$i+=200)
	{
		if($i<=$length && $length<$i+200)
		{
			my $end = $i+200;
			my $range = "$i"."-".$end;
			#print OUT "$type\t$range\t$length\n";
			$hash{$range}++;
			last;
		}
	}
}
close IN;
foreach my $range (keys %hash) {
	print OUT "$range\t$hash{$range}\n";
}

close OUT;

my $dir = dirname($out);
`Rscript $Bin/R/barplot.r --infile $out --outfile $dir/circRNA.length.distribution --x.col 2 --y.col 1 --y.lab "Frequency" --x.lab "Length Range" --skip 1 --title.lab "CircRNA Length Distribution"`;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ", time() - $BEGIN_TIME, "s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

sub ABSOLUTE_DIR {    #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir = `pwd`;
	chomp($cur_dir);
	my ($in) = @_;
	my $return = "";
	if ( -f $in ) {
		my $dir  = dirname($in);
		my $file = basename($in);
		chdir $dir;
		$dir = `pwd`;
		chomp $dir;
		$return = "$dir/$file";
	}
	elsif ( -d $in ) {
		chdir $in;
		$return = `pwd`;
		chomp $return;
	}
	else {
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

################################################################################################################

sub max {    #&max(lists or arry);
	         #求列表中的最大值
	my $max = shift;
	my $temp;
	while (@_) {
		$temp = shift;
		$max = $max > $temp ? $max : $temp;
	}
	return $max;
}

################################################################################################################

sub min {    #&min(lists or arry);
	         #求列表中的最小值
	my $min = shift;
	my $temp;
	while (@_) {
		$temp = shift;
		$min = $min < $temp ? $min : $temp;
	}
	return $min;
}

################################################################################################################

sub revcom() {    #&revcom($ref_seq);
	 #获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
	my $seq = shift;
	$seq =~ tr/ATCGatcg/TAGCtagc/;
	$seq = reverse $seq;
	return uc $seq;
}

################################################################################################################

sub GetTime {
	my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) =
	  localtime( time() );
	return sprintf(
		"%4d-%02d-%02d %02d:%02d:%02d",
		$year + 1900,
		$mon + 1, $day, $hour, $min, $sec
	);
}

################################################################################################################
sub USAGE {
	my $usage = <<"USAGE";
 ProgramName:
     Version:   $version
     Contact:   songmm <songmm\@biomarker.com.cn> 
Program Date:   2016.03.16
      Modify:   
 Description:   This program is used to ......
       Usage:
        Options:
        -i <file>   input file,xxx format,forced
        -o <file>   output dir,forced
        -h      help

USAGE
	print $usage;
	exit;
}
