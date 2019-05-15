#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd qw(realpath);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Times = localtime();
my $year=$Times[5]+1900;
my $month=$Times[4]+1;
my $day=$Times[3];

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fin,$fout,$fkey);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fin,
				"o:s"=>\$fout,
				"k:s"=>\$fkey,
				) or &USAGE;
&USAGE unless ($fin and $fkey);
$fin = Cwd::realpath($fin);
$fout ||= './';
$fout = Cwd::realpath($fout);
mkdir($fout) if (!-d $fout);

my %evalue;
my @evalues = ('1E-5to1E-50','1E-50to1E-100','1E-100to1E-150','1E-150to0','0');
my %identity;
my @identities = ('23%-40%','40%-60%','60%-80%','80%-100%','100%');
my $lines_num = 0;
open IN,$fin or die "cannot open file $fin, $!\n";
open OUT1,">$fout/$fkey.evalue.stat" or die "cannot open file $fout/$fkey.evalue.stat, $!\n";
open OUT2,">$fout/$fkey.identity.stat" or die "cannot open file $fout/$fkey.evalue.stat,$!\n";

while (<IN>) {
	chomp;
	next if (/^$/||/\#/);
	$lines_num ++;
	my @elements = split /\t+/,$_;

	if ($elements[2] > 1e-50) {
		$evalue{$evalues[0]}++;
	}elsif($elements[2] > 1e-100){
		$evalue{$evalues[1]}++;
	}elsif($elements[2] > 1e-150){
		$evalue{$evalues[2]}++;
	}elsif($elements[2] > 0){
		$evalue{$evalues[3]}++;
	}else{
		$evalue{$evalues[4]}++;
	}

	my ($ele) = $elements[3] =~ /\(([\d\.]+)\)$/; 
	$ele = $elements[3] if (!defined $ele && $elements[3] <= 100) ;
	if ($ele < 40 ) {
		$identity{$identities[0]}++;
	}elsif($ele < 60 ) {
		$identity{$identities[1]}++;
	}elsif($ele < 80){
		$identity{$identities[2]}++;
	}elsif($ele < 100) {
		$identity{$identities[3]}++;
	}else{
		$identity{$identities[4]}++;
	}
}
die "please check you file ,what is no line\n" if ($lines_num == 0);
print OUT1 "\#Range\tNum\tPercent\n";
foreach  (@evalues) {
	$evalue{$_} ||= 0;
	my $percent = sprintf("%.2f", 100 * $evalue{$_}/$lines_num);
	print OUT1 "$_\t$evalue{$_}\t$percent\%\n";
}

print OUT2 "\#Range\tNum\tPercent\n";
foreach  (@identities) {
	$identity{$_} ||= 0;
	my $percent = sprintf("%.2f", 100 * $identity{$_}/$lines_num);
	print OUT2 "$_\t$identity{$_}\t$percent\%\n";
}
close IN;
close OUT1;
close OUT2;

`perl $Bin/Just_Pie.pl -i $fout/$fkey.evalue.stat -o $fout/$fkey.evalue.svg -note "Pie for E-value" `;
`perl $Bin/Just_Pie.pl -i $fout/$fkey.identity.stat -o $fout/$fkey.identity.svg -note "Pie for Identity" `;

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub AbsolutePath {   #get file realpath
	my $input = shift;
	my $return = Cwd::realpath("$input");
	return $return;
}

sub USAGE {
	my $usage=<<"USAGE";
	Program:$Script
	Version:$version	[$month:$day:$year]
	Contact:Sun Huaiyu <sunhy\@biomarker.com.cn>
	Options:
		-i	<file>	input file                      [required]
		-k	<str>	the prefix                      [required]
		-o	<file>	outdir  default [./]            [optional]
		-h	            This help

USAGE
	print $usage;
	exit;
}
