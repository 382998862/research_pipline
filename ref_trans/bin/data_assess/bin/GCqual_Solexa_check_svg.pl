#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my %config=%{readconf("$Bin/../../../config/db_file.cfg")}; 
my $program_name=basename($0);

if (@ARGV<3) {
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
mkdir $out unless (-d "$out");
mkdir "$out/PNG" unless (-d "$out/PNG");
$out =~s/\/$//;
print "$qu\n$gc\n";


my $name=basename$qu;
my $Rscript=$config{Rscript};
my $cmd="$Rscript $Bin/quality_bar.R infile=$qu outfile=$out/PNG/$name.png";
print $cmd,"\n";
system$cmd;

$name=basename$gc;
$cmd="$Rscript $Bin/plot_acgtn.R infile=$gc outfile=$out/PNG/$name.png";
print $cmd,"\n";
system$cmd;

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

sub svg_paper (){#&svg_paper(width,height,[color])
	#my $svg_drawer = getlogin()."@".(`hostname`);
	my $svg_drawer = "chenx"."@"."biomarker\.com\.cn";
	chomp $svg_drawer;
	my @svg_x=@_;
	my $line="";
	$line.="<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n";
	$line.="<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 20001102//EN\" \"http://www.w3.org/TR/2000/CR-SVG-20001102/DTD/svg-20001102.dtd\">\n\n";
	$line.="<svg width=\"$svg_x[0]\" height=\"$svg_x[1]\">\n";
	$line.="<Drawer>$svg_drawer</Drawer>\n";
	$line.="<Date>".(localtime())."</Date>\n";
	if (defined $svg_x[2]) {
		$line.="<rect x=\"0\" y=\"0\" width=\"$svg_x[0]\" height=\"$svg_x[1]\" fill=\"$svg_x[2]\"/>\n";
	}
	return $line;
}

sub max 
{
	my ($x1,$x2)=@_;
	my $max;
	if ($x1 > $x2) {
		$max=$x1;
	}
	else {
		$max=$x2;
	}
	return $max;
}

sub min 
{
	my ($x1,$x2)=@_;
	my $min;
	if ($x1 < $x2) {
		$min=$x1;
	}
	else {
		$min=$x2;
	}
	return $min;
}

sub color_gradient  #datanow,data1,data2,color1,color2
{
	my @svg_x=@_;
	my $out_color;
	if ($svg_x[0] >=$svg_x[2]) {
		$out_color=$svg_x[4];
	}
	elsif ($svg_x[0] <=$svg_x[1]) {
		$out_color=$svg_x[3];
	}
	else {
		my $tmp_red1=&hex2ten(substr($svg_x[3],1,2));
		my $tmp_gre1=&hex2ten(substr($svg_x[3],3,2));
		my $tmp_blu1=&hex2ten(substr($svg_x[3],5,2));
		my $tmp_red2=&hex2ten(substr($svg_x[4],1,2));
		my $tmp_gre2=&hex2ten(substr($svg_x[4],3,2));
		my $tmp_blu2=&hex2ten(substr($svg_x[4],5,2));
		my $new_red=int(($svg_x[0]-$svg_x[1])/($svg_x[2]-$svg_x[1])*($tmp_red2-$tmp_red1)+$tmp_red1);
		my $new_gre=int(($svg_x[0]-$svg_x[1])/($svg_x[2]-$svg_x[1])*($tmp_gre2-$tmp_gre1)+$tmp_gre1);
		my $new_blu=int(($svg_x[0]-$svg_x[1])/($svg_x[2]-$svg_x[1])*($tmp_blu2-$tmp_blu1)+$tmp_blu1);
		$new_red=&ten2hex($new_red);$new_red="0$new_red" if(length($new_red)==1);
		$new_gre=&ten2hex($new_gre);$new_gre="0$new_gre" if(length($new_gre)==1);
		$new_blu=&ten2hex($new_blu);$new_blu="0$new_blu" if(length($new_blu)==1);
		$out_color="#$new_red$new_gre$new_blu";
	}
	return $out_color;
}

sub ten2hex  #ten
{
	my $tmp_ten=$_[0];   
	my $hex_value=uc(sprintf("%lx",$tmp_ten));
	return $hex_value;
}

sub hex2ten  #hex
{
	my $tmp_hex=$_[0];
	my $ten_value=0;
	my $tmp_i=0;
	my $tmp_j=0;
	my @tmp_x=split(//,$tmp_hex);
	my %hash_hex=(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,'A',10,'B',11,'C',12,'D',13,'E',14,'F',15);
	for ($tmp_i=@tmp_x-1;$tmp_i>=0 ;$tmp_i--) {
		$ten_value+=$hash_hex{$tmp_x[$tmp_i]}*(16**$tmp_j);
		$tmp_j++;
	}
	return $ten_value;
}

sub svg_polygon  #colorfill,colorstroke,coloropacity,point1,point2,...
{
	my @svg_x=@_;
	my $svg_color=shift(@svg_x);
	my $svg_color2=shift(@svg_x);
	my $svg_trans=shift(@svg_x);
	my $svg_points=join(" ",@svg_x);
	my $line="<polygon fill=\"$svg_color\" stroke=\"$svg_color2\" opacity=\"$svg_trans\" points=\"$svg_points\"/>\n";
	return $line;
}

sub svg_circle  #&svg_circle(x,y,r,color,[info])
{
	my @svg_x=@_;
	my $line="<circle r=\"$svg_x[2]\" cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" fill=\"$svg_x[3]\" />\n";
	if (defined $svg_x[4]) {
		$line="<circle r=\"$svg_x[2]\" cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" fill=\"$svg_x[3]\" onclick=\"alert('$svg_x[4]')\" onmousemove=\"window.status='$svg_x[4]'\" />\n";
	}
	return $line;
}

sub svg_txt  #&svg_txt(x,y,size,color,text,[vertical,0/1/2/3]);
{
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
	my $line="<text fill=\"$svg_x[3]\" transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" font-family=\"ArialNarrow-Bold\" font-size=\"$svg_x[2]\">$svg_x[4]</text>\n";
	return $line;
}

sub svg_mid_txt #&svg_mid_txt(x,y,size,color,text,[vertical,0/1/2/3]);
{
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
	my $line="<text fill=\"$svg_x[3]\" transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" text-anchor=\"middle\" font-family=\"ArialNarrow-Bold\" font-size=\"$svg_x[2]\">$svg_x[4]</text>\n";
	return $line;
}

sub svg_dashed  #&svg_line(x1,y1,x2,y2,color,"10 5",[width])
{
	my @svg_x=@_;
	my $line="<line x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\" style=\"stroke-dasharray:$svg_x[5];fill:none;stroke:$svg_x[4]\"/>\n";
	if (defined $svg_x[6]) {
		$line="<line x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\" style=\"stroke-dasharray:$svg_x[5];fill:none;stroke:$svg_x[4];stroke-width:$svg_x[6]\"/>\n";
	}
	return $line;
}
sub svg_line  #&svg_line(x1,y1,x2,y2,color,[width])
{
	my @svg_x=@_;
	my $line="<line fill=\"$svg_x[4]\" stroke=\"$svg_x[4]\" x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\"/>\n";
	if (defined $svg_x[5]) {
		$line="<line fill=\"$svg_x[4]\" stroke=\"$svg_x[4]\" stroke-width=\"$svg_x[5]\" x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\"/>\n";
	}
	return $line;
}

sub svg_rect  #&svg_rest(x,y,width,height,color,[opacity])
{
	my @svg_x=@_;
	if (!defined $svg_x[5]) {
		$svg_x[5]=1;
	}
	my $line="<rect x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" fill=\"$svg_x[4]\" opacity=\"$svg_x[5]\"/>\n";
	return $line;
}

sub svg_end  #end
{
	return "</svg>\n";
}
