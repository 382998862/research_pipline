#!/usr/bin/perl -w
use	strict;
use	Getopt::Long;
use	Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;

my $program_name=basename($0);
my %config=%{readconf("$Bin/../../../config/db_file.cfg")}; 
my $ver="1.0";
############################################
my %opts;
GetOptions(\%opts,"i=s","o=s");
if (!defined($opts{i})||!defined $opts{o}) {
	
	print << "	Usage End.";
	Description:
	Contact: He hua <heh\@biomarker.com.cn>
	version:$ver
	Usage:

		-i          infile                           must be given;
		-o         outfile key (out.cycleQ.*)        must be given;
		
	Usage End.
		exit;
}
###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $program_name Time :[$Time_Start]\n\n";
############################################
my $in=$opts{i};
my $out=$opts{o};
my $cycle_Q20=0;
open (IN,"$in")||die "$!";
my %cycle;
while (<IN>) {
	chomp;
	next if(/^\#/||/^$/);
	my @line=split/\s+/,$_;
	my $cycle_Q=0;
	for (my $i=1;$i<@line ;$i++) {
		$cycle_Q+=(($i-1)*$line[$i]/100);
	}
	$cycle{$line[0]}=$cycle_Q;
	$cycle_Q20++ if ($cycle_Q>=20);
}
close IN;
my $cycle_num=keys %cycle;
my $cycle_Q20_percent=100*$cycle_Q20/$cycle_num;
my $x=$cycle_num+1;

my $xstep;
if ($cycle_num%6==0) {
	$xstep=$cycle_num/6;
}
elsif ($cycle_num%7==0) {
	$xstep=$cycle_num/7;
}
elsif ($cycle_num%5==0) {
	$xstep=$cycle_num/5;
}

if ($cycle_num%7!=0 && $cycle_num%5!=0 && $cycle_num%6!=0) {
	if ($cycle_num>=180) {
		$xstep=50;
	}
	elsif ($cycle_num>=150 && $cycle_num<180) {
		$xstep=30;
	}
	elsif($cycle_num>90 && $cycle_num<150){
		$xstep=20;
	}
	elsif($cycle_num<=90){
		$xstep=15;
	}
}

open OUT,">$out.cycleQ.stat"||die;
print OUT "$cycle_Q20_percent\n";
close OUT;
open (OUT1,">$out.cycleQ.list")||die "$!";
print OUT1 << "Usage End.";
Type:Rect
Width:800#图宽
Height:600 #图高
WholeScale:0.9
XStep:$xstep
YStep:10
XScalePos:0 #x坐标标尺的起始位置
XStart:0
YStart:0
XEnd:$x
FontSize:36
YEnd:46
Note:Cycle Average Phred Score
X:Cycle Number
Y:Phred Score
XUnit:0.8 #柱形宽度(无重叠时)。
Scale:
0
Usage End.
for (my $i=1;$i<=$x/$xstep ;$i++) {
	print OUT1 $i*$xstep,"\n";
}
print OUT1 ":End\nColor:\#9153F6\n";

for (my $i=1;$i<=$cycle_num;$i++) {
	print OUT1 $i,":$cycle{$i}\n";
}
print OUT1 ":End\n";
close OUT1;

`perl $config{distributing_svg1} $out.cycleQ.list $out.cycleQ.svg`;

################Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $program_name Time :[$Time_End]\n\n";

###############Subs
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
