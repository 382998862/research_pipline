#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $program_name=basename($0);

if (@ARGV<2) {
	print "Usage: perl $0 qufile gcfile outdir\n";
	exit;
}
###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $program_name Time :[$Time_Start]\n\n";
############################################

my @units;
my $i = 0;
my $j = 0;
my %info;
my $tmp;
my $qu = $ARGV[0];
my $gc = $ARGV[1];
my $out= $ARGV[2];
mkdir("$out")unless(-d "$out");
print "$qu\n$gc\n";

open(IN,"$qu")||die"$!";
$tmp = <IN>;
while (<IN>) {
	$i++;
	@units = split;
	for ($j=1;$j<@units ;$j++) {
		$info{$i}{$j} = $units[$j];
	}
}
close IN;
my @tmp=keys %info;
my $len=600/@tmp;


my $graph_qu = basename($qu).".svg";
print "$graph_qu\n";
my $key1;
my $key2;

my @YMark=qw(0 5 10 15 20 25 30 35 40);


open(O,">$out/$graph_qu")||die"$!";
print O <<_FLAG_;
<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg width="800px" height="640px" version="1.1"	xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 640">
  <desc>Solexa Data GC Distribution</desc>
  <desc>lpan\@bmk.com</desc>
  <g font-family="Arial" font-size="30"  font-weight="bold" stroke-width="2" fill="none">
	<g id="top" transform="translate(100,0)">
		<rect width="800" height="80" stroke-fill="black"/>
		<text x="300" y="60" fill="black" font-size="45" style="text-anchor:middle">Quality Distribution</text>
	</g>
	<g id="content" transform="translate(100,80)">
		<rect x="0" y="0" width="608" height="480" stroke-width="4" stroke="black"/>
			<g id="grid">
				<text x="-40" y="10" fill="black">$YMark[8]</text>
				<text x="-40" y="130" fill="black">$YMark[6]</text>
				<text x="-40" y="250" fill="black">$YMark[4]</text>
				<text x="-40" y="370" fill="black">$YMark[2]</text>
				<text x="-40" y="490" fill="black">$YMark[0]</text>

				<line x1="0" y1="120" x2="10" y2="120" stroke="black" stroke-width="4"/>
				<line x1="0" y1="240" x2="10" y2="240" stroke="black" stroke-width="4"/>
				<line x1="0" y1="360" x2="10" y2="360" stroke="black" stroke-width="4"/>
			</g>
			<g id="story" fill="#0000BB">
_FLAG_

my $dis;
if (@tmp%6==0) {
	$dis=@tmp/6;
}
elsif (@tmp%7==0) {
	$dis=@tmp/7;
}
elsif (@tmp%5==0) {
	$dis=@tmp/5;
}

if (@tmp%7!=0 && @tmp%5!=0 && @tmp%6!=0) {
	if (@tmp>=180) {
		$dis=30;
	}
	elsif (@tmp>=150 && @tmp<180) {
		$dis=25;
	}
	elsif(@tmp>90 && @tmp<150){
		$dis=20;
	}
	elsif(@tmp<=90){
		$dis=15;
	}
}
my $step=int(@tmp/$dis);

my @XMark;
for ($i=0; $i<=$step; $i++) {
	$XMark[$i] = $dis * $i;
}
for (my $n=0;$n<=$step;$n++) {
	my $xtxt=-15+$n*608*$dis/@tmp;
	print O "<text x=\"$xtxt\" y=\"508\" fill=\"black\">$XMark[$n]</text>\n";
	my $xline;
	if ($n!=$step) {
		$xline=$n*608*$dis/@tmp+608*$dis/@tmp;
		print O "<line x1=\"$xline\" y1=\"470\" x2=\"$xline\" y2=\"480\" stroke=\"black\" stroke-width=\"4\"/>\n";
	}
}

my $var1;
my $var2;
my $var3;
my $var4;
my @num;
my $num;
my $space;

foreach  $key1 (sort {$a<=>$b} keys %info) {
	foreach $key2 (sort {$a<=>$b} keys %{$info{$key1}}) {
		$var1 = $key1 * $len;
		$var2 = (41 - $key2) * 12;
		$var3 = sqrt($info{$key1}{$key2} / 100);
		$var4=$len-1;
		print O "<rect x=\"$var1\" y=\"$var2\" width=\"$var4\" height=\"11\" opacity=\"$var3\" /> \n";

	}
}


print O <<_FLAG_;
		</g>
	</g>
	<g id="bottom" transform="translate(100,560)">
		<text x="300" y="70" fill="black" font-size="40" style="text-anchor:middle">Cycle Number</text>
	</g>
	<g id="left" transform="translate(0,0)">
		<text x="-120" y="0" style="text-anchor:middle" fill="black" transform="rotate(-90 120,80) " font-size="40">Quality Score</text>
	</g>
  </g>
</svg>
_FLAG_
close O;

##########################################################
#Drawing The Reads GC distribution
################Define main variables here:
my @names;
%info = ();
my %color = (
	"A(%)"=>"green",
	"T(%)"=>"red",
	"C(%)"=>"blue",
	"G(%)"=>"black",
	"N(%)"=>"gray"
);
################
open(IN,"$gc")||die"Can't open $gc\n";
@names = split(/\s+/,<IN>);
my $reads_len;
while (<IN>) {
	@units = split;
	$reads_len++;
	$info{$names[1]}{$units[0]} = $units[1];
	$info{$names[2]}{$units[0]} = $units[2];
	$info{$names[3]}{$units[0]} = $units[3];
	$info{$names[4]}{$units[0]} = $units[4];
	$info{$names[5]}{$units[0]} = $units[5];
}
close IN;

my $graph_gc = basename($gc).".svg";
print"$graph_gc\n";

@YMark=qw(0 10 20 30 40 50);

open(O,">$out/$graph_gc")||die"$!";
print O <<_FLAG_;
<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg width="800px" height="640px" version="1.1"	xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 640">
  <desc>Solexa Data GC Distribution</desc>
  <desc>lpan\@bmk.com</desc>
  <g font-family="Arial" font-size="30"  font-weight="bold" stroke-width="2" fill="none">
	<g id="top" transform="translate(100,0)">
		<rect width="800" height="80" stroke-fill="black"/>
		<text x="300" y="60" fill="black" font-size="45" style="text-anchor:middle">Base Distribution</text>
	</g>
	<g id="content" transform="translate(100,80)">
		<rect x="0" y="0" width="600" height="480" stroke-width="4" stroke="black"/>
			<g id="grid">
				<text x="-40" y="10" fill="black">$YMark[5]</text>
				<text x="-40" y="106" fill="black">$YMark[4]</text>
				<text x="-40" y="202" fill="black">$YMark[3]</text>
				<text x="-40" y="298" fill="black">$YMark[2]</text>
				<text x="-40" y="394" fill="black">$YMark[1]</text>
				<text x="-40" y="490" fill="black">$YMark[0]</text>

				<line x1="0" y1="96" x2="10" y2="96" stroke="black" stroke-width="4"/>
				<line x1="0" y1="192" x2="10" y2="192" stroke="black" stroke-width="4"/>
				<line x1="0" y1="288" x2="10" y2="288" stroke="black" stroke-width="4"/>
				<line x1="0" y1="384" x2="10" y2="384" stroke="black" stroke-width="4"/>

			</g>
			<g id="right" transform="translate(480,0)">
				<rect x="20" y="10" width="90" height="170" stroke="black" stroke-width="2"/>
				<g transform="translate(30,45)" font-family="Arial" font-weight="bold" font-size="30" fill="none">
					<text x="0" y="0" fill="green">A</text>
					<line x1="25" y1="-15" x2="60" y2="-15" stroke="green" stroke-width="4" />
					<text x="0" y="30" fill="red">T</text>
					<line x1="25" y1="15" x2="60" y2="15" stroke="red" stroke-width="4"/>
					<text x="0" y="60" fill="blue">C</text>
					<line x1="25" y1="45" x2="60" y2="45" stroke="blue" stroke-width="4"/>
					<text x="0" y="90" fill="black">G</text>
					<line x1="25" y1="75" x2="60" y2="75" stroke="black" stroke-width="4"/>
					<text x="0" y="120" fill="gray">N</text>
					<line x1="25" y1="105" x2="60" y2="105" stroke="gray" stroke-width="4"/>
				</g>
			</g>
			<g id="story">
_FLAG_
if ($reads_len%6==0) {
	$dis=$reads_len/6;
}
elsif ($reads_len%7==0) {
	$dis=$reads_len/7;
}
elsif ($reads_len%5==0) {
	$dis=$reads_len/5;
}

if ($reads_len%7!=0 && $reads_len%5!=0 && $reads_len%6!=0) {
	if ($reads_len>=180) {
		$dis=30;
	}
	elsif ($reads_len>=150 && $reads_len<180) {
		$dis=25;
	}
	elsif($reads_len>90 && $reads_len<150){
		$dis=20;
	}
	elsif($reads_len<=90){
		$dis=15;
	}
}
$step=int(@tmp/$dis);
@XMark=qw(0 1 2 3 4 5 6 7 8 9 10);
for ($i=0; $i<=$step; $i++) {
	$XMark[$i] = $i*$dis;
}
for (my $n=0;$n<=$step;$n++) {
	my $xtxt=-15+$n*600*$dis/@tmp;
	print O "<text x=\"$xtxt\" y=\"508\" fill=\"black\">$XMark[$n]</text>\n";
	my $xline;
	if ($n!=$step) {
		$xline=$n*600*$dis/@tmp+600*$dis/@tmp;
		print O "<line x1=\"$xline\" y1=\"470\" x2=\"$xline\" y2=\"480\" stroke=\"black\" stroke-width=\"4\"/>\n";
	}
}


foreach  $key1 (sort keys %info) {
	print O "<polyline points=\"";
	foreach $key2 (sort {$a<=>$b} keys %{$info{$key1}}) {
		$var1 = $key2 * $len;
		$var2 = (50 - $info{$key1}{$key2}) * 9.6;
		print O "$var1,$var2\t";
	}
	print O "\" stroke-width=\"2\" stroke=\"$color{$key1}\"\/>\n";
}

print O <<_FLAG_;
		</g>
	</g>
	<g id="bottom" transform="translate(100,560)">
		<text x="300" y="70" fill="black" font-size="40" style="text-anchor:middle">Cycle Number</text>
	</g>
	<g id="left" transform="translate(0,0)">
		<text x="-120" y="0" style="text-anchor:middle" fill="black" transform="rotate(-90 120,80) " font-size="40">Percentage</text>
	</g>
  </g>
</svg>
_FLAG_
close O;
###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $program_name Time :[$Time_End]\n\n";
##################################
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


###################################################################
#Convert Svg To Png 
#my $svg2xxx = dirname($0)."/svg_kit/svg2xxx.pl";
#`perl $svg2xxx ./ graph`;