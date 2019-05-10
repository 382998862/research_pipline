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
my ($geneExpression,$outdir,$config);
GetOptions(
				"help|?" =>\&USAGE,
				"g:s"=>\$geneExpression,
				"c:s"=>\$config,
				"o:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($geneExpression and $outdir and $config);

######################################################################################
# ------------------------------------------------------------------
# Main Body
# ------------------------------------------------------------------

mkdir $outdir if (!-d $outdir);
$geneExpression = ABSOLUTE_DIR($geneExpression);
$outdir = ABSOLUTE_DIR($outdir);
$config = ABSOLUTE_DIR($config);

# ------------------------------------------------------------------
# Get sample names
# ------------------------------------------------------------------

my (@samples,$project);

open CONFIG, "<", "$config" or die "Can't open $config:$!\n";
while (<CONFIG>) {
	chomp;
	next if (/^$/ || /^\#/);
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

my $flag=0;
my @statMapped;
my @posi =  ('left','center','right');

open EXP, ">", "$outdir/geneExpression.html" or die "Can't creat $outdir/geneExpression.html:$!\n";

&PART_1();

for my $sample (@samples) {
	open IN, "<", "$geneExpression/$sample->[1].mappedStat.xls" or die "Can't open $geneExpression/$sample->[1].mappedStat.xls:$!\n";
	my @statMapped;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		push @statMapped,[split/\t/,$_,-1];
	}
	close IN;
	&Table3(\@statMapped,$posi[$flag],$sample);
	$flag = $flag==2 ? 0 : ++$flag;
	
}
#print EXP "                </div>\n" if (@samples % 3 != 0);
print EXP "                </div>\n" if ($flag != 0);

&PART_2();

my @xls;

open XLS, "<", "$geneExpression/$samples[0][1].geneExpression.xls" or die "Can't open $geneExpression/$samples[0][1].geneExpression.xls:$!\n";
while (<XLS>) {
	last if $. == 5;
	chomp;
	s/^\#//;
	push @xls,[split/\t/,$_,-1];
}
close XLS;

&Example(\@xls);

&PART_3();

close EXP;

# ------------------------------------------------------------------
# Encode GB2312 to UTF-8
# ------------------------------------------------------------------
`mv $outdir/geneExpression.html $outdir/geneExpression.html.tmp`;
`perl -MEncode -ne 'print encode("utf8",decode("gb2312",\$_))' $outdir/geneExpression.html.tmp > $outdir/geneExpression.html`;
`rm $outdir/geneExpression.html.tmp`;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
# ------------------------------------------------------------------
# subroutine for HTML
# ------------------------------------------------------------------
sub PART_1{
print EXP <<HTML;
<!DOCTYPE HTML>
<html>
<head>
    <meta http-equiv="Content-Type" content="text/html" charset="utf-8" />
    <title>测序数据及参考基因组统计</title>

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

                <h3>比对效率统计</h3>

                <div class="bcon">
                    <p>不同样品的测序数据与参考基因组进行比对，统计个样品的比对效率、插入片段及打断随机性，对各染色体的覆盖及Reads的分布作图。采用tophat软件进行比对及可变剪接分析，Cufflinks软件进行基因表达丰度的分析，对样品的比对效率统计如下表：</p>
                </div>
HTML
}

sub Table3{
	my ($s,$p,$sam) = @_;
	print EXP "                <div id=\"Table3\">\n" if ($p eq 'left');
	print EXP <<HTML;
                    <div class="$p">

                        <div class="bcon">
                            <p>$sam->[1]：</p>
                        </div>

                        <table id="customers">
                            <tr>
                                <td>$s->[0][0]</td>
                                <td>$s->[0][1]</td>
                                <td>$s->[0][2]</td>
                            </tr>
                            <tr>
                                <th colspan="3"></th>
                            </tr>
                            <tr>
                                <td>$s->[1][0]</td>
                                <td>$s->[1][1]</td>
                                <td>$s->[1][2]</td>
                            </tr>
                            <tr>
                                <td>$s->[2][0]</td>
                                <td>$s->[2][1]</td>
                                <td>$s->[2][2]</td>
                            </tr>
                            <tr>
                                <th colspan="3"></th>
                                </tr>
                            <tr>
                                <td colspan="3">$s->[3][0]</td>
                            </tr>
                            <tr>
                                <td>$s->[4][0]</td>
                                <td>$s->[4][1]</td>
                                <td>$s->[4][2]</td>
                            </tr>
                            <tr>
                                <td>$s->[5][0]</td>
                                <td>$s->[5][1]</td>
                                <td>$s->[5][2]</td>
                            </tr>
                            <tr>
                                <td>$s->[6][0]</td>
                                <td>$s->[6][1]</td>
                                <td>$s->[6][2]</td>
                            </tr>
                            <tr>
                                <th colspan="3"></th>
                            </tr>
                            <tr>
                                <td colspan="3">$s->[7][0]</td>
                            </tr>
                            <tr>
                                <td>$s->[8][0]</td>
                                <td>$s->[8][1]</td>
                                <td>$s->[8][2]</td>
                            </tr>
                            <tr>
                                <th colspan="3"></th>
                            </tr>
                            <tr>
                                <td>$s->[9][0]</td>
                                <td>$s->[9][1]</td>
                                <td>$s->[9][2]</td>
                            </tr>
                            <tr>
                                <th colspan="3"></th>
                            </tr>
                            <tr>
                                <td>$s->[10][0]</td>
                                <td>$s->[10][1]</td>
                                <td>$s->[10][2]</td>
                            </tr>
                            <tr>
                                <td>$s->[11][0]</td>
                                <td>$s->[11][1]</td>
                                <td>$s->[11][2]</td>
                            </tr>
                            <tr>
                                <th colspan="3"></th>
                            </tr>
                            <tr>
                                <td>$s->[12][0]</td>
                                <td>$s->[12][1]</td>
                                <td>$s->[12][2]</td>
                            </tr>
                            <tr>
                                <td>$s->[13][0]</td>
                                <td>$s->[13][1]</td>
                                <td>$s->[13][2]</td>
                            </tr>
                        </table>
                    </div>
HTML
	print EXP "                </div>\n" if ($p eq 'right');
}

sub PART_2{
    print EXP <<HTML;
                <div class="bcon">
                    <p>根据与参考基因组比对的Reads的分布，对每个样品分别作插入片段，随机性，在染色体的分布图。</p>
HTML
    for my $sample(@samples) {
        print EXP <<HTML;
                        <p>$sample->[1]：</p>
                        <ul>
                            <li>
                                <a target="_blank" href="../geneExpression/$sample->[1].insertSize.png">插入片段分布图；</a>
                                <a target="_blank" href="../geneExpression/$sample->[1].randcheck_per.png">随机性检验图；</a>
                                <p>目录: <a target="_blank" href="../geneExpression/$sample->[1].png">染色体Reads及depth分布图</a></p>
                            </li>
                        </ul>
HTML
    }
    print EXP <<HTML;

                </div>

            </div>

            <div class="block">

                <h3>基因表达丰度统计</h3>

                <div class="bcon">

                    <p>基因表达量采用RPKM（Read Per Kb per Million Fragments）来衡量，RPKM对测序总量和基因长度基因长度进行了归一化，其计算公式为：</p>
                    <p><img src="images/RPKM.png"></p>
                    <p>结果示例</p>
HTML
}

sub Example{
my $x = shift;
print EXP <<HTML;
                    <table id="customers">

                        <tr>
                            <th>GeneID</th>
                            <th>Length</th>
                            <th>fpkm</th>
                            <th>Gene_position</th>
                            <th>Strand</th>
                            <th>raw_fragment_number</th>
                            <th>normalized_count</th>
                        </tr>
                        <tr>
                            <td>$x->[1][0]</td>
                            <td>$x->[1][1]</td>
                            <td>$x->[1][2]</td>
                            <td>$x->[1][3]</td>
                            <td>$x->[1][4]</td>
                            <td>$x->[1][5]</td>
                            <td>$x->[1][6]</td>
                        </tr>
                        <tr>
                            <td>$x->[2][0]</td>
                            <td>$x->[2][1]</td>
                            <td>$x->[2][2]</td>
                            <td>$x->[2][3]</td>
                            <td>$x->[2][4]</td>
                            <td>$x->[2][5]</td>
                            <td>$x->[2][6]</td>
                        </tr>
                        <tr>
                            <td>$x->[3][0]</td>
                            <td>$x->[3][1]</td>
                            <td>$x->[3][2]</td>
                            <td>$x->[3][3]</td>
                            <td>$x->[3][4]</td>
                            <td>$x->[3][5]</td>
                            <td>$x->[3][6]</td>
                        </tr>
                        <tr>
                            <th>基因ID</th>
                            <th>基因长度</th>
                            <th>基因坐标</th>
                            <th>表达量值</th>
                            <th>正负链</th>
                            <th>read数</th>
                            <th>标准化read数</th>
                        </tr>
                    </table>
HTML
}

sub PART_3{
    print EXP <<HTML;
                    <ul>
HTML
    for my $sample (@samples) {
    print EXP <<HTML;
                        <li><a target="_blank" href = "../geneExpression/$sample->[1].geneExpression.xls">$sample->[1]基因表达量表</a></li>
HTML
    }
    print EXP <<HTML;
                    </ul>

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
		The progarm is used to creat geneExpression.htm
	Version: $ver
	
	perl $0 -g geneExpression -c config -o output
	
	-m	mapped stat dir
	-g	geneExpression dir
	-c	config
	-o	output dir
	-h	Help

	Usage End.
	exit;
}
