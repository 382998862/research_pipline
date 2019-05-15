#!/usr/bin/perl -w
use strict;
use Cwd;
use SVG;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0";
my @Times = localtime();
#######################################################################################

my %CFG=%{readconf("$Bin/../../CFG")};

my $Time_Start = sub_format_datetime(localtime(time()));
print STDOUT "$Script start at:[$Time_Start]\n";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($funiq1,$funiq2,$foutdir,$fsample1,$fsample2);
GetOptions(
				"help|?" =>\&USAGE,
				"1:s"=>\$funiq1,
				"2:s"=>\$funiq2,
				"s1:s"=>\$fsample1,
				"s2:s"=>\$fsample2,
				"od:s"=>\$foutdir,
				) or &USAGE;
&USAGE unless ($funiq1 and $funiq2 and $foutdir and $fsample1 and $fsample2);
$foutdir = Cwd::realpath($foutdir);
mkdir $foutdir if (!-d $foutdir);
my %tmp;
my %common;
my %uniq_s1;
my %uniq_s2;
my $total_s1 = 0;
my $total_s2 = 0;
open IN,$funiq1 or die "cannot open file $funiq1\n";
while (<IN>) {
	chomp;
	next if (/^\s*$/);$_ =~ s/>//;
	my $line =$_;
	chomp(my $seq = <IN>);
	$line =~ /(.*)_(\d+)_x([\d]+)$/;
	$total_s1 += $3;
	$tmp{$seq} = $line;
	$uniq_s1{$seq} = $line;
}
close IN;
open UNI,$funiq2 or die "cannot open file $funiq2\n";
while (<UNI>) {
	chomp;
	next if (/^\s*$/);$_ =~ s/>//;
	my $line = $_;
	chomp(my $seq = <UNI>);
	$line =~ /(.*)_(\d+)_x([\d]+)$/;
	$total_s2 += $3;
	if (defined $tmp{$seq}) {
		my $key = "$tmp{$seq} $line";
		$common{$seq} = $key;
		delete $uniq_s1{$seq};
		next;
	}
	$uniq_s2{$seq} = $line;
}
close UNI;
my $total_specific_s1 = 0;
my $total_specific_s2 = 0;
open OUTC,">$foutdir/$fsample1\_$fsample2.common.fa" or die;
foreach my $key (keys %common) {
	print OUTC ">$common{$key}\n$key\n";
}
close OUTC;
open OUTU1,">$foutdir/$fsample1.specific.fa" or die;
foreach my $key (keys %uniq_s1) {
	print OUTU1 ">$uniq_s1{$key}\n$key\n";
	$uniq_s1{$key} =~ /(.*)_(\d+)_x([\d]+)$/;
	$total_specific_s1 += $3;
}
close OUTU1;
open OUTU2,">$foutdir/$fsample2.specific.fa" or die;
foreach my $key (keys %uniq_s2) {
	print OUTU2 ">$uniq_s2{$key}\n$key\n";
	$uniq_s2{$key} =~ /(.*)_(\d+)_x([\d]+)$/;
	$total_specific_s2 += $3;
}
close OUTU2;

my $total_total_srna = $total_s1+$total_s2;
my $common_total_srna = $total_s1+$total_s2-$total_specific_s1-$total_specific_s2;
my $uniq_specific_s1 = (keys %uniq_s1);
my $uniq_specific_s2 = (keys %uniq_s2);
my $common_uniq_srna = (keys %common);
my $total_uniq_srna = $uniq_specific_s1+$uniq_specific_s2+$common_uniq_srna;

open STAT,">$foutdir/common_specific.stat" or die;
print STAT "Types\tUnique reads\tPercentage\tTotal reads\tPercentage\n";
print STAT "Total_sRNA\t$total_uniq_srna\t100.00\%\t$total_total_srna\t100.00\%\n";

print STAT "$fsample1\&$fsample2\t$common_uniq_srna\t",(sprintf "%.2f", 100*$common_uniq_srna/$total_uniq_srna),"\%\t$common_total_srna\t",(sprintf "%.2f",100*$common_total_srna/$total_total_srna),"\%\n";

print STAT "$fsample1\_specific\t$uniq_specific_s1\t",(sprintf "%.2f", 100*$uniq_specific_s1/$total_uniq_srna),"\%\t$total_specific_s1\t",(sprintf "%.2f",100*$total_specific_s1/$total_total_srna),"\%\n";

print STAT "$fsample2\_specific\t$uniq_specific_s2\t",(sprintf "%.2f", 100*$uniq_specific_s2/$total_uniq_srna),"\%\t$total_specific_s2\t",(sprintf "%.2f",100*$total_specific_s2/$total_total_srna),"\%\n";

close STAT;

&venn_pairwise($total_uniq_srna,$common_uniq_srna,$uniq_specific_s1,$uniq_specific_s2,"Uniq_sRNA","$fsample1\_vs\_$fsample2.Uniq_sRNA");
&venn_pairwise($total_total_srna,$common_total_srna,$total_specific_s1,$total_specific_s2,"Total_sRNA","$fsample1\_vs\_$fsample2.Total_sRNA");

#######################################################################################
my $Time_End = sub_format_datetime(localtime(time()));
print STDOUT "\n$Script Done at: [$Time_End]\t\tTotal elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub venn_pairwise {
        my ($total,$common,$num1,$num2,$name,$file_name) = @_;
		my $num1_per = sprintf "%.2f",100*$num1/$total;
		my $num2_per = sprintf "%.2f",100*$num2/$total;
		my $common_per = sprintf "%.2f",100*$common/$total;
        open OUT,">$foutdir/$file_name.svg" or die "cannot open $foutdir/$file_name.svg , $!\n";
        my $filename = "$file_name.svg";
        my $svg = SVG->new(width=>800,height=>500);
        $svg->circle(cx=>250,cy=>275,r=>150,style=>{'fill'=>'hotpink','stroke'=>'black','stroke-width'=>0,'fill-opacity'=>0.6});
        $svg->circle(cx=>400,cy=>275,r=>150,style=>{'fill'=>'deepskyblue','stroke'=>'black','stroke-width'=>0,'fill-opacity'=>0.6});
        $name = ucfirst($name);
        $svg->text(x=>400,y=>50,-cdata=>"$name",style=>{'font-size'=>30,'font-weight'=>'bold','text-anchor'=>'middle'});
        $svg->text(x=>120,y=>275,-cdata=>"$num1_per\%",style=>{'font-size'=>20,'font-weight'=>'bold'});
        $svg->text(x=>300,y=>275,-cdata=>"$common_per\%",style=>{'font-size'=>20,'font-weight'=>'bold'});
        $svg->text(x=>450,y=>275,-cdata=>"$num2_per\%",style=>{'font-size'=>20,'font-weight'=>'bold'});

		$svg->rect(x=>560,y=>100,width=>15,height=>15,style=>{'fill'=>'hotpink','fill-opacity'=>0.6});
		$svg->rect(x=>560,y=>130,width=>15,height=>15,style=>{'fill'=>'hotpink','fill-opacity'=>0.6});
		$svg->rect(x=>560,y=>130,width=>15,height=>15,style=>{'fill'=>'deepskyblue','fill-opacity'=>0.6});
		$svg->rect(x=>560,y=>160,width=>15,height=>15,style=>{'fill'=>'deepskyblue','fill-opacity'=>0.6});

        $svg->text(x=>580,y=>111,-cdata=>"$fsample1 specific\($num1\)",style=>{'font-size'=>15});
        $svg->text(x=>580,y=>141,-cdata=>"$fsample1\&$fsample2\($common\)",style=>{'font-size'=>15});
		$svg->text(x=>580,y=>171,-cdata=>"$fsample2 specific\($num2\)",style=>{'font-size'=>15});
        print OUT $svg->xmlify(
                -pubid => "-//W3C//DTD SVG 1.0//EN",
                -inline   => 1
        );
        close OUT;
        chdir $foutdir;
        print "$CFG{svg2xxx}  $filename -t png \n";
        system("$CFG{svg2xxx}  $filename -t png");
}


sub sub_format_datetime {   #Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub USAGE {
	my $usage=<<"USAGE";
	Program:$Script
	Version:$version	[2013/4/16]
	Contact:Sun Huaiyu <sunhy\@biomarker.com.cn>
	Purpose:Get the common and specific sRNA ,make statistics
	Options:
		-1	<file>	input sample1 uniq sRNA file
		-2	<file>	input sample2 uniq sRNA file
		-s1	<str>	the sample1 marker
		-s2	<str>	the sample2 marker
		-od	<file>	output dir
		-h	Help

USAGE
	print $usage;
	exit;
}
