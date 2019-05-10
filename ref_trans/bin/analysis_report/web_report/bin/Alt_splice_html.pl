

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
    <title>�ɱ���ӷ���</title>

    <!-- Linking styles -->
    <link rel="stylesheet" href="css/style.css" type="text/css" media="screen">
    <link rel="stylesheet" href="css/menu.css" type="text/css" media="screen">
</head>
<body>

    <!-- Defining the header section of the page -->
    <header id="header">
        <div class="head_container">
            <div class="logo">
                <a href="../index.html" title="������ҳ"><img src="./images/logo.gif" height=57 /></a>
            </div>
            <nav>
                <!-- Defining the navigation menu -->
                <ul class="topnav">

                    <li>
                        <a href="rawdata_genome.html">��������ͳ��</a>
                    </li>

                    <li>
                        <a href="geneExpression.html">�ȶ�Ч�ʼ���������</a>
                    </li>

                    <li>
                        <a href="Alt_splice.html">�ɱ���Ӽ��»������</a>
                    </li>

                    <li>
                        <a href="deg.html">����������</a>
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

                <h3>�ɱ��������</h3>

                <div class="bcon">
                    <p align="letf">&nbsp;&nbsp;&nbsp;&nbsp;�ɱ������mRNA����ǰ��ӹ�Ϊ����mRNA�������ں��Ӽ���λ��Ĳ�ͬ�γɲ�ͬת¼��������ѧ���󡣲�ͬת¼�����벻ͬ���ף�����ܵ����������ﲻͬ��״���ɱ����һ���Ϊ7�����ͣ��ֱ��ǣ�</p>
                    <p>&nbsp;&nbsp;&nbsp;&nbsp;<font size="3" color="#1430E5">Exon Skipping��</font>���ӹ����������ӱ������ļ������ͣ�<ul><li><a target="_blank" href = "images/exon_skipping.png">���ʾ��</a></li></ul></p>
                    <p>&nbsp;&nbsp;&nbsp;&nbsp;<font size="3" color="#1430E5">Intron retention��</font>���ӹ������ں��ӱ������������ͣ�<ul><li><a target="_blank" href = "images/intron_retention.png">���ʾ��</a></li></ul></P>
                    <p>&nbsp;&nbsp;&nbsp;&nbsp;<font size="3" color="#1430E5">Alt3' splice��</font>ĳ�������ӵ�3'����λ�㷢���˸ı䣻<ul><li><a target="_blank" href = "images/alt_3_site.png">���ʾ��</a></li></ul></p>
                    <p>&nbsp;&nbsp;&nbsp;&nbsp;<font size="3" color="#1430E5">Alt5' splice��</font>ĳ�������ӵ�5'����λ�㷢���˸ı䣻<ul><li><a target="_blank" href = "images/alt_5_site.png">���ʾ��</a></li></ul></p>
                    <p>&nbsp;&nbsp;&nbsp;&nbsp;<font size="3" color="#1430E5">Alt First exon��</font>ת¼����5'�������Ӽ���λ�÷����˱仯��<ul><li><a target="_blank" href = "images/alt_first_exon.png">���ʾ��</a></li></ul></p>
                    <p>&nbsp;&nbsp;&nbsp;&nbsp;<font size="3" color="#1430E5">Alt Last exon��</font>ת¼����3'�������Ӽ���λ�÷����˱仯��<ul><li><a target="_blank" href = "images/alt_last_exon.png">���ʾ��</a></li></ul></p>
                    <p>&nbsp;&nbsp;&nbsp;&nbsp;<font size="3" color="#1430E5">Mutually exclusive exon��</font>���������ӷֱ�����γ�������ͬ��ת¼����<ul><li><a target="_blank" href = "images/mutually_exclusive_exon.png">���ʾ��</a></li></ul></p>
                </div>
            </div>

            <div class="block">

                <h3>�ɱ���ӷ��������</h3>

                <div class="bcon">

                    <p>&nbsp;&nbsp;&nbsp;&nbsp;ÿ��������Ʒ��Reads��ο�������ȶԣ�ͨ���ȶ����ȷ������λ�ã��Ƚϼ���λ�õı仯�Ӷ�ȷ����ͬ�Ŀɱ�������ͣ��������tophat<a target="_blank" href="http://tophat.cbcb.umd.edu/index.html">[1]</a>����������бȶԼ����ӷ�����ͨ��splicegrapher<a target="_blank" href="http://splicegrapher.sourceforge.net/">[2]</a>���ͳ�ƿɱ�������ͼ���ͼ������</p>
                    <p>������ʾ����<ul><li><a target="_blank" href = "Alt_splice_info.html">���ʾ��</a></li></ul></p>
                    <p>��ͼ�����<ul><li><a target="_blank" href = "images/altsplice.png">��ͼ���</a></li></ul></p>

                    <p>��Ŀ�ɱ���ӽ��Ŀ¼��Alt_statΪ�ɱ���ӽ��ͳ�Ʊ��Ŀ¼��pdfΪ�ɱ������ͼĿ¼��</p>
HTML
}

sub href{
	my $ref = shift;
	print ALT <<HTML;
                    <p>$ref</p>
                    <ul>
                        <li>
                            <a target="_blank" href="../Alt_splice/$ref">$ref �ɱ���ӽ��Ŀ¼��</a>
                        </li>
                    </ul>
HTML
}

sub newGene {
    print ALT <<HTML;
                </div>

            </div>
            <div class="block">
                <h3>�»��򷢾�</h3>

                <div class="bcon">

                    <p>&nbsp;&nbsp;&nbsp;&nbsp;ÿ��������Ʒ��Reads��ο�������ȶԣ�ͨ�������������ע�ͻ�����������ıȽϣ����򱾴β��������һЩ����������Ϊ���β����ֵ��»�����Ϣ�������ṩ�»����������Ϣ���Լ���ʽ����gffע���ļ�</p>
                    <p>���Ŀ¼��<ul><li><a target="_blank" href = "../NewGene">���Ŀ¼</a></li></ul></p>
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
            <div class="container">����������Ƽ����޹�˾ &copy; 2012 | <a href="#">BMK</a></div>
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
