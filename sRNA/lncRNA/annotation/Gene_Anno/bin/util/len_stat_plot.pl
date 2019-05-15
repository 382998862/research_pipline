#!/usr/bin/perl -w
# 
#Copyright (c) BMK 2011 
#Writer           Mengf <mengf@biomarker.com.cn>
#Program Date   2011 
#Modifier         Mengf <mengf@biomarker.com.cn>
#Last modified  2011 
my $version="1.2.0";
my $BEGIN=time();

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;

my %config=%{readconf("$Bin/../../db_file.cfg")} ;

my $programe_dir=basename($0);
my $path=dirname($0);
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($Index,$fa,$region,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$Index,
				"s=s"=>\$region,
				"fa=s"=>\$fa,
				"od=s"=>\$od,
				) or &USAGE;
&USAGE unless ($Index and $fa and $od) ;

$fa=&ABSOLUTE_DIR($fa);
&MKDIR($od);
$od=&ABSOLUTE_DIR($od);

$region = $region || "300,500,1000,2000";

my @range=split/,/,$region;

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $programe_dir Time :[$Time_Start]\n\n";

################################ Program
my %length_region_dis;
my $count=0;

open IN,"$fa" || die $!;
$/='>';
while (<IN>) {
	chomp;
	s/^\s+//;s/\s+$//;s/\r+$//;
	next if (/^$/ || /^\#/);
	my ($head,$seq)=split/\n+/,$_,2;
	$count++;
	my $id=(split/\s+/,$head)[0];
	$seq=~s/\s+//g;
	my $len=length($seq);

	my $region_site=int($len/100);
	if ($region_site>=30) {
		$length_region_dis{30}++;
	}
	else {
		$length_region_dis{$region_site}++;
	}
}

close IN;

#####################
#Stat output
################################

open OUT2,">$od/$Index.stat.xls" || die $!;
print OUT2 "#Length_span\tNumbers\tPercent\n";
foreach (0..30) {
	my $site=$_;
	my $min=$site*100;
	my $max=($site+1)*100;
	my $len_span="$min"."~"."$max";
	my $percent;
	if ($site==30 && defined $length_region_dis{30}) {
		$len_span=">3000";
		$percent=($length_region_dis{30}/$count)*100;
		printf OUT2 "%s\t%s\t\t%.2f\n",$len_span,$length_region_dis{30},$percent;
	}
	elsif (!defined $length_region_dis{$site}) {
		print OUT2 "$len_span\t0\t0\n";
	}
	else {
		$percent=($length_region_dis{$site}/$count)*100;
		printf OUT2 "%s\t%s\t%.2f\n",$len_span,$length_region_dis{$site},$percent;
	}
}
close OUT2;


#####################
#Length Distribution Graph
###############################

open OUT,">$od/$Index.distribution.svg" || die $!;
my $svg=&Length_Distribution_svg(\%length_region_dis,$Index);
print OUT "$svg";
close OUT;
chdir $od;
system "perl $config{svg2xxx1} $Index.distribution.svg";  #svg-->png

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $programe_dir Time :[$Time_End]\n\n";
&Runtime($BEGIN);


####################subs
sub Length_Distribution_svg {#
	my ($len_dis,$lable) = @_;
	my $svg;
	my %length_dis=%$len_dis;
	my $max_value=(sort {$a <=> $b} values %length_dis)[-1];
	my $y_max=int(&log10($max_value))+1;
	my $width=720;
	my $height=480;
	my $left_pos=100;my $right_pos=$width-35;
	my $low_pos=$height-100;my $high_pos=50;
	my $yx_coordinate=$left_pos-5;my $yy_coordinate=$low_pos;
	my $xx_coordinate=$left_pos;my $xy_coordinate=$low_pos+5;
	my $x_iner=3;my $x_span=9;
	my $rect_color = "rgb(0,165,210)";

	###############画坐标轴 坐标轴刻度 网格线
	$svg.=&svg_paper($width,$height);
	$svg.=&svg_line($left_pos,$low_pos,$left_pos,$high_pos,"#000000");
	$svg.=&svg_txt_middle($left_pos/2,$height/2,16,"#000000","$lable Number",3);
	$svg.=&svg_line($left_pos,$low_pos,$right_pos,$low_pos,"#000000");
	$svg.=&svg_txt_middle($width/2,$height-40,16,"#000000","Length (nt)");
	$svg.=&svg_txt_middle($width/2,$high_pos/2+10,18,"#000000","$lable Length Distribution");

	my $y_step=($low_pos-$high_pos)/$y_max;
	my $x_step=($right_pos-$left_pos)/31;
	for (0..$y_max) {
		if ($_==0) {
			$svg.=&svg_line($yx_coordinate,$low_pos-($y_step*$_),$xx_coordinate,$low_pos-($y_step*$_),"#000000");
			$svg.=&svg_txt_end($xx_coordinate-10,$low_pos-($y_step*$_)+5,10,"#000000",0);
		}
        elsif ($_<=5) {
			$svg.=&svg_line($yx_coordinate,$low_pos-($y_step*$_),$xx_coordinate,$low_pos-($y_step*$_),"#000000");
			$svg.=&svg_txt_end($xx_coordinate-10,$low_pos-($y_step*$_)+5,10,"#000000",10**$_);
		}
		else {
			$svg.=&svg_line($yx_coordinate,$low_pos-($y_step*$_),$xx_coordinate,$low_pos-($y_step*$_),"#000000");
			$svg.=&svg_txt_end($xx_coordinate-10,$low_pos-($y_step*$_)+5,10,"#000000","1e+$_");
		}
	}
	for (1..30) {
		$svg.=&svg_line($left_pos+($x_step*$_),$yy_coordinate,$left_pos+($x_step*$_),$xy_coordinate,"#000000");
        my $lab = 100*($_-1)."-".100*$_;
		$svg.=&svg_txt_end($left_pos+($x_step*$_)+5-0.5*$x_step,$xy_coordinate+13-0.5*$x_step,8,"#000000",$lab,4);
#		$svg.=&svg_txt_middle($left_pos+($x_step*$_)+5,$xy_coordinate+13,10,"#000000",$lab,4);
	}
	$svg.=&svg_txt_end($right_pos-0.4*$x_step,$xy_coordinate+15-0.5*$x_step,8,"#000000",">3000",4);

	##############画长度分布柱子
	for my $site (0..30) {
		next if (!defined $length_dis{$site}) ;
		my $x=$left_pos+($site)*$x_step+$x_iner;
		my $log=&log10($length_dis{$site});
		my $y=($log/$y_max)*($low_pos-$high_pos);
		$svg.=&svg_rect_event($x,$low_pos-$y,$x_span,$y,$rect_color,$rect_color,1);
	}
	$svg.=&svg_end();
	return ($svg);
}

sub svg_txt (){
	#&svg_txt(x,y,size,color,text,[vertical,0/1/2/3]);						# 写字符文档（x坐标，y坐标，字符大小，颜色，字符文档，，）
	#$anchor 以xy确定的点为准，0：text以点起始；1：以点为中心；2：以点结束。
	#$svg_x[6] 旋转，0:0度;1:90度;2:180度;3:270度;
	my @svg_x=@_;
	if (!defined $svg_x[6]) {
		$svg_x[6]=0;
	}
	my $svg_matrix='';
	if ($svg_x[6]==0) {
		$svg_matrix="1 0 0 1";
	}
	if ($svg_x[6]==1) {
		$svg_matrix="0 1 -1 0";
	}
	if ($svg_x[6]==2) {
		$svg_matrix="-1 0 0 -1";
	}
	if ($svg_x[6]==3) {
		$svg_matrix="0 -1 1 0";
	}
	if (!defined $svg_x[5] || $svg_x[5] == 0) {
		my $line="<text fill=\"$svg_x[3]\"  transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" font-family=\"Arial\" font-size=\"$svg_x[2]\">$svg_x[4]</text>\n";
		return $line;
	}else{
		my $anchor="";
		if ($svg_x[5]==1) {
			$anchor="middle";
		}
		if ($svg_x[5]==2) {
			$anchor="end";
		}
		my $line="<text fill=\"$svg_x[3]\" text-anchor=\"$anchor\" transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" font-family=\"Arial\" font-size=\"$svg_x[2]\">$svg_x[4]</text>\n";
		return $line;
	}
}

sub svg_txt_middle (){#&svg_txt(x,y,size,color,text,[vertical,0/1/2/3]);
	my @svg_x=@_;
	if (!defined $svg_x[5]) {
		$svg_x[5]=0;
	}
	my $svg_matrix='';
	if ($svg_x[5]==0) {
		$svg_matrix="1 0 0 1";
	}
	if ($svg_x[5]==1) {
		$svg_matrix="0 1 -1 0";
	}
	if ($svg_x[5]==2) {
		$svg_matrix="-1 0 0 -1";
	}
	if ($svg_x[5]==3) {
		$svg_matrix="0 -1 1 0";
	}
	my $line="<text fill=\"$svg_x[3]\" transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" font-family=\"ArialNarrow-Bold\" font-size=\"$svg_x[2]\" text-anchor=\"middle\">$svg_x[4]</text>\n";
	return $line;
}

sub svg_txt_end (){#&svg_txt(x,y,size,color,text,[vertical,0/1/2/3]);
	my @svg_x=@_;
	if (!defined $svg_x[5]) {
		$svg_x[5]=0;
	}
	my $svg_matrix='';
	if ($svg_x[5]==0) {
		$svg_matrix="1 0 0 1";
	}
	if ($svg_x[5]==1) {
		$svg_matrix="0 1 -1 0";
	}
	if ($svg_x[5]==2) {
		$svg_matrix="-1 0 0 -1";
	}
	if ($svg_x[5]==3) {
		$svg_matrix="0 -1 1 0";
	}
    if ($svg_x[5]==4) {
		$svg_matrix="1 -0.75 0.75 1";
    }
	my $line="<text fill=\"$svg_x[3]\" transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" font-family=\"ArialNarrow-Bold\" font-size=\"$svg_x[2]\" text-anchor=\"end\">$svg_x[4]</text>\n";
	return $line;
}

sub svg_line (){#&svg_line(x1,y1,x2,y2,color,[width])
	my @svg_x=@_;
	my $line="<line fill=\"$svg_x[4]\" stroke=\"$svg_x[4]\" x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\"/>\n";
	if (defined $svg_x[5]) {
		$line="<line fill=\"$svg_x[4]\" stroke=\"$svg_x[4]\" stroke-width=\"$svg_x[5]\" x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\"/>\n";
	}
	return $line;
}

sub svg_paper (){#&svg_paper(width,height,[color])
	my @svg_x=@_;
	my $line="";
	$line.="<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n";
	$line.="<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 20001102//EN\" \"http://www.w3.org/TR/2000/CR-SVG-20001102/DTD/svg-20001102.dtd\">\n\n";
	$line.="<svg width=\"$svg_x[0]\" height=\"$svg_x[1]\" viewBox=\"0 0 $svg_x[0] $svg_x[1]\">\n";
	$line.="<Date>".(localtime())."</Date>\n";
	if (defined $svg_x[2]) {
		$line.="<rect x=\"0\" y=\"0\" width=\"$svg_x[0]\" height=\"$svg_x[1]\" fill=\"$svg_x[2]\"/>\n";
	}
	return $line;
}

sub svg_rect_event () {#&svg_rect(x,y,width,height,fill_color,stroke_color,strock_width);
	my @svg_x=@_;
	my $line ;
	$line="<rect x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" style=\"fill:$svg_x[4];stroke:$svg_x[5];stroke-width:$svg_x[6];\" />\n";
	if (defined $svg_x[7]){
		$line="<rect x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" style=\"fill:$svg_x[4];stroke:$svg_x[5];stroke-width:$svg_x[6];fill-opacity:$svg_x[7]\" />\n";
	}
	if (defined $svg_x[8]){
		$line="<rect x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" style=\"fill:$svg_x[4];stroke:$svg_x[5];stroke-width:$svg_x[6];fill-opacity:$svg_x[7]\" onclick=\"alert('$svg_x[8]')\" onmousemove=\"window.status='$svg_x[8]'\" />\n";
	}
	if (defined $svg_x[9]){
		$line="<rect x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" style=\"fill:$svg_x[4];stroke:$svg_x[5];stroke-width:$svg_x[6];fill-opacity:$svg_x[7]\" onmouseover=\"changeText(evt,'$svg_x[8]','$svg_x[9]')\" onmouseout=\"changeTextNotOver(evt,'$svg_x[8]')\" onclick=\"changeClick(evt)\" onmousemove=\"window.stsatus='$svg_x[8]'\"  attributeName=\"fill-opacity\" from=\"0\" to=\"0.4\" begin=\"mouseover\" end=\"mouseout\"/> \n";
	}
	return $line;
}

sub svg_end (){#
	return "</svg>\n";
}

sub log10 {#
	my ($n) = @_;
	return log($n)/log(10);
}

sub Runtime
{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total elapsed time: ${t}s\n";
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub GetTMR {#
	#Get Total Mapped Reads from file line which contain string "Mapped Reads\s"
	my $fStat=shift;
	open (IN,"<",$fStat) or die $!;
	while (<IN>) {
		if (/Mapped Reads\s(\d+)/) {
			close (IN) ;
			return $1;
		}
	}
	close (IN) ;
	die "Error Reads Stat file.\n";
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
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
sub USAGE {#
    my $usage=<<"USAGE";
     Program: $0
     Version: $version
     Contact: Meng Fei <mengf\@biomarker.com.cn>
Program Date: 2012-07-02
      Modify: Simon Young <simonyoung8824\@gmail.com> 
 Modify Date: 2014-08-08
 Description: 
        Usage:
          -i       <str>    Index of inFiles and outFiles       must be given;
          -fa      <str>    Infile                              must be given;
          -od      <str>    OUT file DIR                        must be given;

USAGE
	print $usage;
	exit;
}
