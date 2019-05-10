

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
my ($alt,$outdir,$config);
GetOptions(
				"help|?" =>\&USAGE,
				"d:s"=>\$alt,
				"c:s"=>\$config,
				"o:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($alt and $outdir and $config);

######################################################################################
# ------------------------------------------------------------------
# Main Body
# ------------------------------------------------------------------

mkdir $outdir if (!-d $outdir);
$alt = ABSOLUTE_DIR($alt);
$outdir = ABSOLUTE_DIR($outdir);
$config = ABSOLUTE_DIR($config);

# ------------------------------------------------------------------
# Get sample names
# ------------------------------------------------------------------

my (@samples,$project);

open CONFIG, "<", "$config" or die "Can't open $config:$!\n";
while (<CONFIG>) {
	if (/^Sample/){
		s/\s+$|\n$//;
		push @samples,[split/\s+/,$_,-1];
	}
	if (/^Project/) {
		chomp;
		$project = (split/\s+/,$_,-1)[1];
	}
}
close CONFIG;

my $sam_num=$#samples+1;

open ALT, ">", "$outdir/Alt_splice.html" or die "Can't creat $outdir/Alt_splice.html:$!\n";

&PART_1();
my $i=0;
for my $ref (glob "$alt/*") {
	&href(basename($ref));
	$i++;
}

if ($i!=$sam_num) {
	print "error: Some samples no Altsplice Results, please Check!\n";die;
}

&newGene();

&PART_2();

close ALT;

# ------------------------------------------------------------------
# Encode GB2312 to UTF-8
# ------------------------------------------------------------------
`mv $outdir/Alt_splice.html $outdir/Alt_splice.html.tmp`;
`perl -MEncode -ne 'print encode("utf8",decode("gb2312",\$_))' $outdir/Alt_splice.html.tmp > $outdir/Alt_splice.html`;
`rm $outdir/Alt_splice.html.tmp`;

#######################################################################################
# ------------------------------------------------------------------
# subroutine for HTML
# ------------------------------------------------------------------
sub PART_1{
print ALT <<HTML;
<!DOCTYPE HTML>
<html>
<head>
    <meta http-equiv="Content-Type" content="text/html" charset="utf-8" />
    <title>可变剪接分析</title>

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

                <h3>可变剪接类型</h3>

                <div class="bcon">
                    <p align="letf">&nbsp;&nbsp;&nbsp;&nbsp;可变剪接是mRNA在由前体加工为成熟mRNA过程中内含子剪切位点的不同形成不同转录本的生物学现象。不同转录本编码不同蛋白，则可能导致生物体表达不同性状。可变剪接一般分为7种类型，分别是：</p>
                    <p>&nbsp;&nbsp;&nbsp;&nbsp;<font size="3" color="#1430E5">Exon Skipping：</font>剪接过程中外显子被跳过的剪接类型；<ul><li><a target="_blank" href = "images/exon_skipping.png">结果示例</a></li></ul></p>
                    <p>&nbsp;&nbsp;&nbsp;&nbsp;<font size="3" color="#1430E5">Intron retention：</font>剪接过程中内含子被保留剪接类型；<ul><li><a target="_blank" href = "images/intron_retention.png">结果示例</a></li></ul></P>
                    <p>&nbsp;&nbsp;&nbsp;&nbsp;<font size="3" color="#1430E5">Alt3' splice：</font>某个外显子的3'剪接位点发生了改变；<ul><li><a target="_blank" href = "images/alt_3_site.png">结果示例</a></li></ul></p>
                    <p>&nbsp;&nbsp;&nbsp;&nbsp;<font size="3" color="#1430E5">Alt5' splice：</font>某个外显子的5'剪接位点发生了改变；<ul><li><a target="_blank" href = "images/alt_5_site.png">结果示例</a></li></ul></p>
                    <p>&nbsp;&nbsp;&nbsp;&nbsp;<font size="3" color="#1430E5">Alt First exon：</font>转录本的5'端外显子剪接位置发生了变化；<ul><li><a target="_blank" href = "images/alt_first_exon.png">结果示例</a></li></ul></p>
                    <p>&nbsp;&nbsp;&nbsp;&nbsp;<font size="3" color="#1430E5">Alt Last exon：</font>转录本的3'端外显子剪接位置发生了变化；<ul><li><a target="_blank" href = "images/alt_last_exon.png">结果示例</a></li></ul></p>
                    <p>&nbsp;&nbsp;&nbsp;&nbsp;<font size="3" color="#1430E5">Mutually exclusive exon：</font>两个外显子分别剪接形成两条不同的转录本；<ul><li><a target="_blank" href = "images/mutually_exclusive_exon.png">结果示例</a></li></ul></p>
                </div>
            </div>

            <div class="block">

                <h3>可变剪接分析及结果</h3>

                <div class="bcon">

                    <p>&nbsp;&nbsp;&nbsp;&nbsp;每个测序样品的Reads与参考基因组比对，通过比对情况确定剪接位置，比较剪接位置的变化从而确定不同的可变剪接类型，具体采用tophat<a target="_blank" href="http://tophat.cbcb.umd.edu/index.html">[1]</a>软件进行序列比对及剪接分析，通过splicegrapher<a target="_blank" href="http://splicegrapher.sourceforge.net/">[2]</a>软件统计可变剪接类型及作图分析。</p>
                    <p>结果表格示例：<ul><li><a target="_blank" href = "Alt_splice_info.html">结果示例</a></li></ul></p>
                    <p>作图结果：<ul><li><a target="_blank" href = "images/altsplice.png">作图结果</a></li></ul></p>

                    <p>项目可变剪接结果目录：Alt_stat为可变剪接结果统计表格目录，pdf为可变剪接做图目录；</p>
HTML
}

sub href{
	my $ref = shift;
	print ALT <<HTML;
                    <p>$ref</p>
                    <ul>
                        <li>
                            <a target="_blank" href="../Alt_splice/$ref">$ref 可变剪接结果目录；</a>
                        </li>
                    </ul>
HTML
}

sub newGene {
    print ALT <<HTML;
                </div>

            </div>
            <div class="block">
                <h3>新基因发掘</h3>

                <div class="bcon">

                    <p>&nbsp;&nbsp;&nbsp;&nbsp;每个测序样品的Reads与参考基因组比对，通过与基因组已有注释基因序列区域的比较，发掘本次测序特异的一些编码区域，作为本次测序发现的新基因信息，这里提供新基因的序列信息，以及格式化的gff注释文件</p>
                    <p>结果目录：<ul><li><a target="_blank" href = "../NewGene">结果目录</a></li></ul></p>
HTML
}

sub PART_2{
    print ALT <<HTML;
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
		The progarm is used to creat DEG.htm
	Version: $ver
	
	perl $0 -d DEG_Analysis -c config -o output
	
	-d	Alt splice dir
	-c	config
	-o	output dir
	-h	Help

	Usage End.
	exit;
}
