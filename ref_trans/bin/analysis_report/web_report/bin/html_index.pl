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
my ($config,$outdir);
GetOptions(
				"help|?" =>\&USAGE,
				"c:s"=>\$config,
				"o:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($config and $outdir);
######################################################################################
# ------------------------------------------------------------------
# Main Body
# ------------------------------------------------------------------

mkdir $outdir if (!-d $outdir);
$config = ABSOLUTE_DIR($config);
$outdir = ABSOLUTE_DIR($outdir);

my ($Project_name,$Project_info,$Contract_NO,$Customer_info);

open CONFIG, "<", "$config" or die "Can't open $config:$!\n";
	while (<CONFIG>) {
		next if (/#/ || /^\s+$/);
		chomp;
		$Project_name = (split/\s+/,$_)[1] if /Project_name/;
		$Project_info = (split/\s+/,$_)[1] if /Project_info/;
		$Contract_NO = (split/\s+/,$_)[1] if /Contract_NO/;
		$Customer_info = (split/\s+/,$_)[1] if /Customer_info/;
	}

close CONFIG;

open INDEX, ">", "$outdir/index.html" or die "Can't creat $outdir/index.html:$!\n";


print INDEX <<HTML;
<!DOCTYPE HTML>
<html lang="en-US">
<head>
    <meta http-equiv="Content-Type" content="text/html" charset="utf-8" />
    <title>BMK technology 项目报告</title>

    <!-- Linking styles -->
    <link rel="stylesheet" href="html/css/index_style.css" type="text/css" media="screen">
    <link rel="stylesheet" href="hyml/css/index_menu.css" type="text/css" media="screen">

    <!-- Linking scripts -->
    <!--[if lt IE 9]>
        <script src="http://html5shim.googlecode.com/svn/trunk/html5.js"></script>
    <![endif]-->
    <script src="html/js/jquery.js" type="text/javascript"></script>
    <script src="html/js/superfish.js" type="text/javascript"></script>
    <script src="html/js/supersubs.js" type="text/javascript"></script>
    <script src="html/js/script.js" type="text/javascript"></script>
    <script src="html/js/jquery.nivo.slider.js" type="text/javascript"></script>
</head>
<body>

    <!-- Defining the header section of the page -->
    <header id="header">
        <div class="container">
            <div class="logo">
                <a href="../index.html" title="返回首页"><img src="html/images/logo.gif" height=57 /></a>
            </div>

        </div>
    </header>



    <!-- Defining the main content section of the page -->
    <div class="primary_content container">

        <!-- column #1 -->
        <div id="column1" class="column">

            <div class="block">
                <h3>项目报告</h3>
                <div class="bcon">
                    <p>项目名称：$Project_name</p>
					<br />
                    <p>项目详细：$Project_info</p>
                    <br />
					<p>合同编号：$Contract_NO</p>
                    <br />
					<p>客户信息：$Customer_info</p>
                </div>

            </div>



        </div>

        <!-- column #2 -->
        <div id="column2" class="column">

            <div class="block2">
                <div class="bg">
                    <h3>信息分析结果</h3>
                    <div class="indexbcon">
                        <ul>
                            <li>
                                <a target="_blank" href="html/rawdata_genome.html" class="date"><span><span>测序数据统计</span></span></a>
                                <p>目录:<a target="_blank" href="rawdata">rawdata</a></p>
                            </li>
                            <li>
                                <a target="_blank" href="html/geneExpression.html" class="date"><span><span>比对效率及基因表达量</span></span></a>
                                <p>目录:<a target="_blank" href="geneExpression">geneExpression</a></p>
                            </li>
                            <li>
                                <a target="_blank" href="html/Alt_splice.html" class="date"><span><span>可变剪接及新基因分析</span></span></a>
                                <p>可变剪接目录:<a target="_blank" href="Alt_splice">Alt_Splice</a></p>
                                <p>新基因目录:<a target="_blank" href="NewGene">NewGene</a></p>
                            </li>
                            <li>
                                <a target="_blank" href="html/deg.html" class="date"><span><span>差异表达分析</span></span></a>
                                <p>目录:<a target="_blank" href="DEG_Analysis"> DEG_Analysis</a></p>
                            </li>
                        </ul>
                    </div>
                    <div>
                        <a href="#" class="button"><span><span>HOME</span></span></a>
                    </div>
                </div>
                <div class="bot"></div>
            </div>

        </div>

    </div>

    <footer><!-- Defining the footer section of the page -->


    <div class="copyright">
        <div class="container">百迈客生物科技有限公司 &copy; 2012 | <a href="#">BMK</a></div>
    </div>
    </footer>
</body></html>
HTML
close INDEX;

# ------------------------------------------------------------------
# Encode GB2312 to UTF-8
# ------------------------------------------------------------------
`mv $outdir/index.html $outdir/index.html.tmp`;
`perl -MEncode -ne 'print encode("utf8",decode("gb2312",\$_))' $outdir/index.html.tmp > $outdir/index.html`;
`rm $outdir/index.html.tmp`;

# ------------------------------------------------------------------
# Creat css and js folder
# ------------------------------------------------------------------

if (!-d "$outdir/html/css"){
    `mkdir -p $outdir/html/css`;
}
if (!-d "$outdir/html/js"){
    `mkdir -p $outdir/html/js`;
}
if (!-d "$outdir/html/images"){
    `mkdir -p $outdir/html/images`;
}

`cp -r $Bin/bin/css $outdir/html/`;
`cp -r $Bin/bin/js $outdir/html`;
`cp -r $Bin/bin/images $outdir/html`;


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

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
		The progarm is used to creat index.htm
	Version: $ver
	
	perl $0 -c Config.txt -o output
	
	-c	config file
	-o	output dir
	-h	Help

	Usage End.
	exit;
}
