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
my ( $fIn, $dir);
## input is CIRI output file;
GetOptions(
	"help|?"   => \&USAGE,
	"od:s"    => \$dir,
	"i:s"      => \$fIn
) or &USAGE;
&USAGE unless ( $fIn and $dir);
my @files = split /,/,$fIn;
my %hash;
foreach my $file(@files)
{
	open( IN, $file ) or die $!;
	<IN>;
	my $prefix = basename($file);
	$prefix=~/(^\w+)\./;
	open( OUT,  ">$dir/$1_length.stat" )or die $!;
	print OUT "circRNA_type\trange\tgeneLength\n";
	while (<IN>) {
		next if (/\#/);
		my @line        = split /\t/, $_;
		my $length      = $line[2] - $line[1] + 1;
		my $type = $line[4];
		for(my $i=0;;$i+=250)
		{
			if($i<=$length && $length<$i+250)
			{
				my $end = $i+250;
				my $range = "$i"."-".$end;
				print OUT "$type\t$range\t$length\n";
				#$hash{$type}{$range}++;
				last;
			}
		}
	}
	close IN;
	close OUT;
	`/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript $Bin/R/barplot.r --infile $dir/$1_length.stat --outfile $dir/$1.circRNA.length.distribution --x.col 3 --y.col 1 --y.lab "The number of circRNA" --x.lab "Length Range" --skip 1 --title.lab "CircRNA Length Distribution"`;
}

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
        -od <file>   output file dir ,forced
        -h      help

USAGE
	print $usage;
	exit;
}
