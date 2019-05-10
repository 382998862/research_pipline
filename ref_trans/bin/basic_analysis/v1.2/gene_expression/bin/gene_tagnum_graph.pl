#!/usr/local/bin/perl -w
# 
# Copyright (c) BMK 2011
# Writer:          xiayh <xiayh@biomarker.com.cn>
# Program Date:   2011.
# Modifier:        xiayh <xiayh@biomarker.com.cn>
# Last Modified:  2011.
my $ver="1.0";

use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my %config=%{readconf("$Bin/../../../../../config/db_file.cfg")}; 

######################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作

my %opts;
GetOptions(\%opts,"i1=s","i2=s","value=s","o=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{i1}) || !defined($opts{i2}) ||!defined($opts{value}) ||!defined($opts{o}) || defined($opts{h}))
{
	&help;
}

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
################
my $programe_dir=basename($0);

my $path=dirname($0);


my  $in1=$opts{i1};
my $in2=$opts{i2} ;#"/share/nas5/RNA_seq/Project/RNA_seq/Cabbage/Analysis/Map_evaluate_20110715/Cabbage_A-WT/Cabbage_A-WT.geneExpression.xls";##这是RPKM值会发生变化
my $threshold= $opts{value};#1;
my $out=$opts{o}.".express";#"/share/nas1/xiayh/Bin/result/new_expre_analysis1.txt";

my $line_num=`wc -l $in1`;
my $read_num=(split /\s+/,$line_num)[0];

my $total_read=$read_num/4; ###统计reads的数目
my $aa;
my @read_total;
my $read_change_num;
for(my $i=0;$i<=80;$i++)
{
	$aa=1-$i/80;
	$read_change_num=$total_read*($aa);
	push @read_total,$read_change_num;

}
@read_total=sort {$b<=>$a} @read_total;##由大到小的排列
my ($x_end,$x_step)=&define_axis ($read_total[0]/1000000);
$x_end=$x_end;

open FI,"$in2";
my @b;

my $RPKM;
my @gene_RPKM;
while(<FI>)
{
	chomp;
	next if(/^\#Gene/);
	@b=split/\s+/,$_;
	$RPKM=$b[2];
	next if($RPKM=~/^FPKM/);
	next if($RPKM==0);
	push	@gene_RPKM	,$RPKM;
}
close FI;

my $total=0;
for(my $kk=0;$kk<=$#gene_RPKM;$kk++)
	{
		if($gene_RPKM[$kk]>=$threshold)
		{
			$total++;
		}
	}
$total=$total/1000;
my ($y_end,$y_step)=&define_axis ($total);
$y_end=$y_end+4;
open FO,">$out.gene_tag.list";

print  FO <<"	Usage End.";
Type:Point
NoLine:0
PointSize:1
Width:700
Height:500
WholeScale:0.85
XStart:0
YStart:0
XEnd:$x_end
YEnd:$y_end
XStep:$x_step
YStep:$y_step
Fontsize:30
Color:#00BB00
x:Total  Tag Number(X1M)
y: Total Gene  Number(X1K)
Note:FPKM>=$threshold

Color#00BB00
Mark:Gene mapped by reads
	Usage End.

my $gene_num=0;

my $gene_RPKM;
my $RP;
my $yend;
for(my $j=0;$j<=80;$j++)
{
	$gene_num=0;
	for(my $k=0;$k<=$#gene_RPKM;$k++)
	{
		$RP=$gene_RPKM[$k]*(1-$j/80);
		if($RP>=$threshold)
		{
			$gene_num++;
		}
	}
	
	#$read_gene{$read_total[$j]}=$gene_num;
	print FO $read_total[$j]/1000000,":",$gene_num/1000,"\n";
	 
}
close FO;

my $od=dirname $out;
my $out_name=basename $out;
chdir $od;
`perl $config{distributing_svg2} $out_name.gene_tag.list  $out_name.gene_tag.svg` ;
`$config{svg2xxx2}  $out_name.gene_tag.svg `;

##返回坐标轴的最大值及相应的步长信息
sub define_axis () {
	my $i=0;
	my ($max)=@_;
	my $time=1;
	my @ret=();
	my @limit=(1,2,3,4,5,6,8,10,12,15,16,20,24,25,30,40,50,60,80,100,120);
	my @unlim=(0,1,2,3,4,5,6,8 ,10,11,14,15,18.5,21,23,29,37,47,56,76 ,92 ,110);
	while ($max >$unlim[21]) {
		$max=$max/10;
		$time=$time*10;
	}
	for ($i=0;$i<=20 ;$i++) {
		if ($max>$unlim[$i] && $max<=$unlim[$i+1]) {
			$ret[0]=$limit[$i]*$time;
			if ($i==2 || $i==5 || $i==9 || $i==14) {
				$ret[1]=$ret[0]/3;
			}
			elsif ($i==4 || $i==7 || $i==13 || $i==16){
				$ret[1]=$ret[0]/5;
			}
			else {
				$ret[1]=$ret[0]/4;
			}
		}
	}
	return @ret;
}



###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

###############Subs
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
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

sub help
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver

	Usage:

		-i1          infile     must be given 

		-i2          infile     must be given

		-value                 the lower limit of RPKM

		-o          outfile    must be given 

		-h    Help document

	Usage End.

	exit;
}

