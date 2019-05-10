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
my ($deg,$outdir,$config);
GetOptions(
				"help|?" =>\&USAGE,
				"d:s"=>\$deg,
				"c:s"=>\$config,
				"o:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($deg and $outdir and $config);

######################################################################################
# ------------------------------------------------------------------
# Main Body
# ------------------------------------------------------------------

mkdir $outdir if (!-d $outdir);
$deg = ABSOLUTE_DIR($deg);
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

open DEG, ">", "$outdir/deg.html" or die "Can't creat $outdir/deg.html:$!\n";

&PART_1();

for my $ref (glob "$deg/*") {
	next if ($ref =~ /All_DEG/);
	&href(basename($ref));
}

&PART_2();

close DEG;

# ------------------------------------------------------------------
# Encode GB2312 to UTF-8
# ------------------------------------------------------------------
`mv $outdir/dge.html $outdir/deg.html.tmp`;
`perl -MEncode -ne 'print encode("utf8",decode("gb2312",\$_))' $outdir/deg.html.tmp > $outdir/deg.html`;
`rm $outdir/deg.html.tmp`;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
# ------------------------------------------------------------------
# subroutine for HTML
# ------------------------------------------------------------------
sub PART_1{
print DEG <<HTML;
<!DOCTYPE HTML>
<html lang="en-US">
<head>
    <meta http-equiv="Content-Type" content="text/html" charset="utf-8" />
    <title>����������</title>

    <!-- Linking styles -->
    <link rel="stylesheet" href="css/style.css" type="text/css" media="screen">
    <link rel="stylesheet" href="css/menu.css" type="text/css" media="screen">
</head>
<body>

    <!-- Defining the header section of the page -->
    <header id="header">
        <div class="head_container">
            <div class="logo">
                <a href="../index.htm" title="������ҳ"><img src="./images/logo.gif" height=57 /></a>
            </div>
            <nav>
                <!-- Defining the navigation menu -->
                <ul class="topnav">

                    <li>
                        <a href="rawdata_genome.html">�������ݼ��ο�������ͳ��</a>
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
                <h3>������������</h3>
                <div class="bcon">
                    <p>���ݲ���Reads��ο��������еıȶԣ��õ���Ӧ��������Ʒ�еı���ȣ������ڲ�ͬ��Ʒ��ı���ȵĲ�ͬ��Ѱ�Ҳ�����Ļ���</p>
                    <p>��ϲ������Ĺ���ע�ͣ�����pathway��KEGG����GO�ĸ����������Լ��Բ�����������ģʽ����ķ�����</p>
                    <p>������������fisher��ȷ���飬ͨ��BonferroniУ��������У�����õ������������������pathway��GO�����࣬�Ա���н�һ���ķ����о���</p>
                </div>
            </div>

            <div class="block">
                <h3>�����������</h3>
                <div class="bcon">
HTML
}

sub href{
	my $ref = shift;
	print DEG <<HTML;
                    <p>$ref</p>
                    <ul>
                        <li>
                            <a target="_blank" href="../DEG_Analysis/$ref/$ref.annotation.xls">��������б�ע�ͣ�</a>
                            <a target="_blank" href="../DEG_Analysis/$ref/$ref.expression.plot.png">Scatterplot��</a>
                            <a target="_blank" href="../DEG_Analysis/$ref/$ref.fold_change.plot.png">log2ratioplot��</a>
                            <a target="_blank" href="../DEG_Analysis/$ref/DEG_Cluster/$ref.DEG.cluster.png">����ͼ��</a>
                        </li>
                        <li>
                            <a target="_blank" href="../DEG_Analysis/$ref/pathway/pathway.html">Kegg pathway������</a>
                        </li>
                        <li>
                            <a target="_blank" href="../DEG_Analysis/$ref/go_enrichment/$ref.GO.png">GOע��ͼ��</a>
                            <a target="_blank" href="../DEG_Analysis/$ref/go_enrichment/${ref}_F1.html">GO(Molecular Function)��</a>
                            <a target="_blank" href="../DEG_Analysis/$ref/go_enrichment/${ref}_P1.html">GO(Biological Process)��</a>
                            <a target="_blank" href="../DEG_Analysis/$ref/go_enrichment/${ref}_C1.html">GO(Cellular Component)��</a>
                        </li>
                        <li>
                            <a target="_blank" href="../DEG_Analysis/$ref/Cog_Anno/$ref.Cog.classfy.png">Cog����ͼ��</a>
                            <a target="_blank" href="../DEG_Analysis/$ref/Cog_Anno/$ref.Cog_class.txt.stat">Cog����ͳ�ƣ�</a>
                        </li>
                    </ul>
HTML
}

sub PART_2{
	if (-d "$deg/All_DEG") {
		print DEG <<HTML;
                    <p>All_DEG</p>
                    <ul>
                        <li>
                            <a target="_blank" href="../DEG_Analysis/All_DEG/pathway/pathway.html">Kegg pathway������</a>
                        </li>
                        <li>
                            <a target="_blank" href="../DEG_Analysis/All_DEG/All_DEG.GO.png">GOע��ͼ��</a>
                            <a target="_blank" href="../DEG_Analysis/All_DEG/go_enrichment/All_DEG_F1.html">GO(Molecular Function)��</a>
                            <a target="_blank" href="../DEG_Analysis/All_DEG/go_enrichment/All_DEG_P1.html">GO(Biological Process)��</a>
                            <a target="_blank" href="../DEG_Analysis/All_DEG/go_enrichment/All_DEG_C1.html">GO(Cellular Component)��</a>
                        </li>
                        <li>
                            <a target="_blank" href="../DEG_Analysis/All_DEG/Cog_Anno/All_DEG.Cog.classfy.png">Cog����ͼ��</a>
                            <a target="_blank" href="../DEG_Analysis/All_DEG/Cog_Anno/All_DEG.Cog_class.txt.stat">Cog����ͳ�ƣ�</a>
                        </li>
                    </ul>
HTML
    }
    print DEG <<HTML;
                </div>

            </div>

            <div class="block">
                <a href="../index.htm" class="button"><span><span>HOME</span></span></a>
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
	
	-d	DEG_Analysis dir
	-c	config
	-o	output dir
	-h	Help

	Usage End.
	exit;
}
