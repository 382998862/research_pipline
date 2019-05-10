#!/usr/bin/perl -w
# 
# Copyright (c) BIO_MARK 2012
# Writer:         Zhangyh <zhangyh@biomarker.com.cn>
# Program Date:   2012/9/3 9:25
# Modifier:       Zhangyh <zhangyh@biomarker.com.cn>
# Last Modified:  2012/9/3

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $ver="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
#GetOptions
# ------------------------------------------------------------------
my ($rawdata,$outdir,$config);
GetOptions(
				"help|?" =>\&USAGE,
				"r:s"=>\$rawdata,
				"c:s"=>\$config,
				"o:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($rawdata and $outdir and $config);
######################################################################################
# ------------------------------------------------------------------
# Main Body
# ------------------------------------------------------------------

mkdir $outdir if (!-d $outdir);
$rawdata = ABSOLUTE_DIR($rawdata);
$outdir = ABSOLUTE_DIR($outdir);
$config = ABSOLUTE_DIR($config);

# ------------------------------------------------------------------
# Get sample names
# ------------------------------------------------------------------

my @samples;
my @posi =  ('left','center','right');
my $flag = 0;

open CONFIG, "<", "$config" or die "Can't open $config:$!\n";
while (<CONFIG>) {
	next unless (/^Sample/);
	s/\s+$|\n$//;
	push @samples,[split/\s+/,$_,-1];
}
close CONFIG;

open RAWDATA, ">", "$outdir/rawdata_genome.html" or die "Can't creat $outdir/rawdata_genome.html:$!\n";

&PART_1();

my $a=0;
open STAT, "<", "$rawdata/AllSample_GC_Q.stat" or die "Can't open $rawdata/AllSample_GC_Q.stat:$!\n";
while (<STAT>) {
	chomp;
	if ($.==1) {
		my @A=split/\s+/,$_;
		$a=@A;
		print RAWDATA "                        </tr>\n" if $a==7;
		#print RAWDATA "                            <th>Q30%</th>\n                        </tr>\n" if $a==8;
		print RAWDATA "                            <th>Q30%</th>\n                        </tr>\n" if $a==8;
		next;
	}
	my @stat = (split/\s+/,$_,-1);
	&Q20stat(\@stat) if $a==7;
	&Q30stat(\@stat) if $a==8;
}

print RAWDATA "                    </table>\n";
print RAWDATA "                    <h3> CycleQ20分布图:</h3>\n";


for my $sample (@samples) {
	print RAWDATA "                    <div id=\"Table3\">\n" if $flag == 0;
	&Q20PNG($posi[$flag],$sample->[1]);
	print RAWDATA "                    </div>\n" if $flag == 2;
	print RAWDATA "                    </div>\n" if ($sample eq $samples[-1] and $flag != 2);
	$flag = $flag==2 ? 0 : ++$flag;
}


&PART_2();

close RAWDATA;


# ------------------------------------------------------------------
# Encode GB2312 to UTF-8
# ------------------------------------------------------------------
`mv $outdir/rawdata_genome.html $outdir/rawdata_genome.html.tmp`;
`perl -MEncode -ne 'print encode("utf8",decode("gb2312",\$_))' $outdir/rawdata_genome.html.tmp > $outdir/rawdata_genome.html`;
`rm $outdir/rawdata_genome.html.tmp`;


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# subroutine for HTML
# ------------------------------------------------------------------

sub PART_1{

	print RAWDATA <<HTML;
<!DOCTYPE HTML>
<html>
<head>
    <meta http-equiv="Content-Type" content="text/html" charset="utf-8" />
    <title>测序数据统计</title>

    <!-- Linking styles -->
    <link rel="stylesheet" href="css/style.css" type="text/css" media="screen">
    <link rel="stylesheet" href="css/menu.css" type="text/css" media="screen">
</head>
<body>

    <!-- Defining the header section of the page -->
    <header id="header">
        <div class="head_container">
            <div class="logo">
                <a href="../index.html" title="返回首页"><img src="./images/logo.gif" height=57 /></a>
            </div>
            <nav>
                <!-- Defining the navigation menu -->
                <ul class="topnav">

                    <li>
                        <a href="rawdata_genome.html">测序数据统计</a>
                    </li>

                    <li>
                        <a href="geneExpression.html">比对效率及基因表达量</a>
                    </li>

                    <li>
                        <a href="Alt_splice.html">可变剪接及新基因分析</a>
                    </li>

                    <li>
                        <a href="deg.html">差异基因分析</a>
                    </li>


                    <li>
                        <a href="../index.html">HOME</a>
                    </li>
                </ul>
            </nav>
        </div>
    </header>
    <!-- Defining the main content section of the page -->
    <div class="primary_content container">

        <!-- column #1 -->
        <div id="column1" class="column">

            <div class="block">

                <div class="bcon">

                    <h3>测序数据统计</h3>

                    <p>采用illumina NGS 测序平台对每个样品进行双端测序，数据量以及数据质量，高质量以及充足的数据量是保证我们后续结果准确性的依据。</p>

                    <table id="customers">
                        <tr>
                            <th>SampleID</th>
                            <th>Reads number</th>
                            <th>basenumber</th>
                            <th>GC%</th>
                            <th>N%</th>
                            <th>Q20%</th>
                            <th>CycleQ20%</th>

HTML
}

sub Q20stat{
	my $stat=shift;
	print RAWDATA <<HTML;
                        <tr>
                            <td>$stat->[0]</td>
                            <td>$stat->[1]</td>
                            <td>$stat->[2]</td>
                            <td>$stat->[3]</td>
                            <td>$stat->[4]</td>
                            <td>$stat->[5]</td>
                            <td>$stat->[6]</td>
                        </tr>

HTML
}
sub Q30stat{
	my $stat=shift;
	print RAWDATA <<HTML;
                        <tr>
                            <td>$stat->[0]</td>
                            <td>$stat->[1]</td>
                            <td>$stat->[2]</td>
                            <td>$stat->[3]</td>
                            <td>$stat->[4]</td>
                            <td>$stat->[5]</td>
                            <td>$stat->[6]</td>
                            <td>$stat->[7]</td>
                        </tr>
HTML
}


sub Q20PNG{
	my ($posi,$sample_name) = @_;
	print RAWDATA <<HTML;
                        <div class="$posi">
                            <div class="bcon">
                                <p>$sample_name：</p>
                            </div>
                            <a target=_blank href="../rawdata/PNG/$sample_name.cycleQ.png"><IMG src="../rawdata/PNG/$sample_name.cycleQ.png" border="0" width="310" height="257"></a>
                        </div>

HTML
}

sub PART_2{
	print RAWDATA <<HTML;
                </div>

            </div>
            <div class="block">
                <a href="../index.html" class="button"><span><span>HOME</span></span></a>
                <a href="#" class="button"><span><span>TOP</span></span></a>
            </div>

        </div>

    </div>


   <footer><!-- Defining the footer section of the page -->


    <div class="copyright">
        <div class="container">百迈客生物科技有限公司 &copy; 2012 | <a href="#">BMK</a></div>
    </div>
    </footer>
</body>
</html>

HTML
}
# ------------------------------------------------------------------
# subroutine
# ------------------------------------------------------------------

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {
	print <<"	Usage End.";
	
	Description:
		The progarm is used to creat rawdata_stat.htm
	Version: $ver
	
	perl $0 -r rawdata -o output
	
	-r	rawdata dir
	-c	config
	-o	output dir
	-h	Help

	Usage End.
	exit;
}
