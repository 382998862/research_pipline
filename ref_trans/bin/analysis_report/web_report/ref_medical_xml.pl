#!/usr/bin/perl
use warnings;

use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);
use File::Basename qw(basename dirname);
use XML::Writer;
use IO::File;
use Encode;

my ($old,$output);
GetOptions(
        "h|?"           =>\&USAGE,
	"oldxml:s"	=>\$old,
	"o:s"		=>\$output,
)or &USAGE;
&USAGE unless ($old);


my $inpath=dirname $old;
$output ||= "$inpath/configtest_medical.xml";
my @images=();
my @tables=();
my ($h1,$h2,$h3,$degid);
open(XML,$old)||die $!;
open(OUT,">$output")||die $!;
my @tmp=();
#open(TEST,">$inpath/test.out")||die $!;
while(my $write=<XML>){
	chomp($write);
	$h1=(split(/\"/,$write))[1]	if($write=~/\<h1/);
	$h2=(split(/\"/,$write))[1]	if($write=~/\<h2/);
	$h3=(split(/\"/,$write))[1]	if($write=~/\<h3/);
	

	if($h3 =~/差异表达基因蛋白互作网络/){
#		print TEST "$h3\n";
		print OUT "$write\n";
		$degid= (split(/\s+/,$h3))[0];
#		print TEST "id is .$degid.\n";
		my ($d1,$d2,$d3)=split(/\./,$degid);
		for (my $i=1;$i<5;$i++){
			$write=<XML>;
			print OUT "$write\n";
		}
		my $newid=join(".",($d1,$d2,($d3+1)));
		$h3="$newid 转录因子注释";
		print OUT &pattern("h3","type1","$newid 转录因子注释","$newid 转录因子注释"),"\n";
		print OUT &pattern("p","type1","转录水平调控是基因表达调控的重要环节，其中转录因子（Transcriptionfactor，TF）通过结合基因上游特异性核苷酸序列实现对该基因的转录调控，许多生物学功能都通过调节特异性转录因子实现对下游功能基因的激活或抑制，因此有必要注释不同分组间差异表达基因的转录组因子,使用动物转录因子数据库AnimalTFDB对差异表达基因中的转录因子进行鉴定。"),"\n";
		print OUT &pattern("p","type1","Transcription_factor分析结果"),"\n";
		print OUT &pattern("file_list","xls","注:ID：Ensembl基因ID；Gene_Symbol：基因通用名称；FPKM：相关样本中该转录因子的表达水平；log2FC：样品间差异倍数变化；Description：基因描述;Family:转录因子所属的家族。"),"\n";
		@tables=glob("$inpath/BMK_5_DEG_Analysis/BMK_*_vs_*/BMK_5_DEG_Anno/*_vs_*.DEG.TF.xls");		
		foreach my $i(@tables){
			my $dir=dirname $i;
			my $base=basename $i;
			my $rel_dir=(split(/BMK_5_DEG_Analysis/,$dir))[1];
			print OUT &pattern("file","xls","注:gene: 基因；gene symbol：基因通用名称；FPKM：相关样本中该转录因子的表达水平；log2FC：样品间差异倍数变化；Description：基因描述;Family:转录因子所属的家族。",$base,"BMK_5_DEG_Analysis$rel_dir/$base","xls"),"\n";
		}	
		print OUT '</file_list>',"\n";

		@tables=glob("$inpath/BMK_5_DEG_Analysis/BMK_*_vs_*/BMK_5_DEG_Anno/*_vs_*.DEG.cosmic.xls");	
		my $cancer_num=@tables;
		if($cancer_num>0){
			my $newid=join(".",($d1,$d2,($d3+2)));
			$h3="$newid 癌症基因功能注释";
			print OUT &pattern("h3","type1","$newid 癌症基因功能注释","$newid 癌症基因功能注释"),"\n";
			print OUT &pattern("p","type1","原癌基因（Proto-oncogene）是参与细胞生长、分裂和分化的正常基因，当其发生突变后（如基因序列被改变）就会变成致癌基因（Oncogene）。通常在肿瘤或恶性细胞系中某些特异性癌基因会上调表达，通过了解癌基因在不同实验组的表达情况有助于深入认识疾病的发病机理。Cosmic(https://cancer.sanger.ac.uk/cosmic)是英国Sanger实验室开发并维护的癌基因及相关注释数据库，有较高的权威性和可信度，通过Cosmic数据库，可对差异表达基因中的癌基因部分进行注释。"),"\n";
			print OUT &pattern("p","type1","Oncogenes分析结果"),"\n";

			print OUT &pattern("file_list","xls","注：gene:基因；gene symbol：基因通用名称；log2FC：样品间差异倍数变化；Description：基因描述；Tumor Type(Somatic)：体细胞癌症类型；Tumor Type(Germline)：生殖细胞系癌症类型。"),"\n";
			foreach my $i(@tables){
				my $dir=dirname $i;
				my $rel_dir=(split(/BMK_5_DEG_Analysis/,$dir))[1];
				my $base=basename $i;
				print OUT &pattern("file","xls",$base,$base,"BMK_5_DEG_Analysis$rel_dir/$base","xls"),"\n";
			}	
			print OUT '</file_list>',"\n";
		}

		next;


	}

	if($write=~/参考文献/){
		my $currenth3=(split(/\s+/,$h3))[0];
		@tmp=split(/\./,$currenth3);
		
		my $currenth2=(split(/\s+/,$h2))[0];
		print "当前的h2标题是$h2\nid是$currenth2\n";
		@tmp=split(/\./,$currenth2);
	if(-d "$inpath/BMK_10_Gene_Fusion"){
		my $newid=join(".",($tmp[0],($tmp[1]+1)));
		print OUT &pattern("h2","type1","$newid 融合基因分析","$newid 融合基因分析"),"\n";
		print OUT &pattern("p","type1","融合基因是指将两个或多个基因的编码区首尾相连，置于同一套调控序列(包括启动子、增强子、核糖体结合序列、终止子等)控制之下，构成的嵌合基因。融合基因的表达产物为融合蛋白。我们使用Fusionmap在转录组中研究基因融合事件。Fusionmap首先通过比对到基因组和转录本中双末端(pairend)关系的序列寻找候选的基因融合，然后采用通过与nt等数据库比较，过滤掉假阳性结果。"),"\n";

		$pattern=&pattern("pic_list","type1","注：红色的线代表同一染色体上发生的融合事件，绿色的线代表不同染色体上发生的融合事件。","检测到的基因融合事件");
		print OUT "$pattern\n";
		@images=glob("$inpath/BMK_10_Gene_Fusion/png/*_circos.png");
		foreach my $i(@images){
			my $base=basename $i;
			print OUT &pattern("pic","type1",$base,$base,"BMK_10_Gene_Fusion/png/$base"),"\n";
		}	
		print OUT '</pic_list>',"\n";
		print OUT &pattern("h3","type1","$newid.1 融合基因统计","$newid.1 融合基因统计"),"\n";
		print OUT &pattern("file_list","xls","注：Fusion ID：基因融合事件编号；Unique Mapping Position：映射到唯一位点的基因融合数目；Count：基因融合位点总数；Gene：发生融合的基因；Strand：基因融合发生在正义链“+”还是反义链“-”；Chromosome：融合基因位于的染色体；Start：基因融合起始位点；End：基因融合终止位点；Filter：In Black List，代表数据库中注释的假阳性基因融合。","融合基因事件统计"),"\n";
		@tables=glob("$inpath/BMK_10_Gene_Fusion/report/*_FusionReport.txt.xls");		
		foreach my $i(@tables){
			my $base=basename $i;
			print OUT &pattern("file","xls",$base,$base,"BMK_10_Gene_Fusion/report/$base","xls"),"\n";
		}		
		print OUT '</file_list>',"\n";
		print OUT &pattern("p","type1","注：Fusion ID：基因融合事件编号；Unique Mapping Position：映射到唯一位点的基因融合数目；Count：基因融合位点总数；Gene：发生融合的基因；Strand：基因融合发生在正义链“+”还是反义链“-”；Chromosome：融合基因位于的染色体；Start：基因融合起始位点；End：基因融合终止位点；Filter：基因融合事件类型；InFamilyList代表基因对属于同一家族；InParalogueList代表基因对来自同一旁系同源组；SameEnsembleGene代表基因对有一个super ensemble基因交集；InBlackList代表基因对包含线粒体和核糖体基因或者假基因。"),"\n";

                my $pfam="BMK_10_Gene_Fusion/Fusion_gene_Pfam.anno.xls";
                if(-e "$inpath/$pfam"){
			print OUT &pattern("h3","type1","$newid.2 融合基因结构域分析","$newid.2 融合基因结构域分析"),"\n";
			print OUT &pattern("p","type1","融合基因可能会导致具有新的或不同功能的基因产物生成，这些异常活性的蛋白质会发挥癌基因还可能与一个新的启动子发生融合，从而激活下游癌基因的表达。后者常见于淋巴肿瘤当中。我们对融合基因进行结构域及相关功能预测，有利于进一步探索其与肿瘤发生、发展或者转移过程中的作用机制。Pfam蛋白结构域搜索结果展示如下："),"\n";
			print OUT &pattern("p","type1","融合基因的Pfam注释结果"),"\n";
			print OUT &pattern("file","xls","注：Fusion Gene_ID：fusionmap产生的融合基因ID;Pfam acc：比对到pfam结构域的ID；Pfam name：pfam结构域名称；Pfam description：pfam结构域详细描述;Bit score：比对打分分值；E-value：比对evalue。","Fusion_gene_Pfam.anno.xls",$pfam,"xls"),"\n";
		}
	}

	}
	print OUT "$write\n";
}


close(XML);
close(OUT);
#close(TEST);


########################################################################################################
#
#						Sub function
########################################################################################################



sub pattern{
	my($title,$type,$desc,$name,$path,$action)=@_;
	my $pattern="\<$title type=\"$type\" desc=\"$desc\" name=\"$name\"";
	$pattern .= " path =\"$path\" "		if($path);
	$pattern .= " action=\"$action\" "	if($action);

	if($title =~/list/){
		$pattern .="\>";
	}else{
		$pattern .="\/\>";
	}
	return $pattern;
}



sub readConfig{
	my $configFile=shift;
	my $d=Config::General->new(-ConfigFile => "$configFile");
	my %config=$d->getall;	
	return %config;
}


sub table_html{
	my ($input,$outHtml,$title,$linkColNum,$linkHash,$srcPath,$text)=@_;
	my @inputs=();
	$/="\n";
	open (IN,$input)||die $!;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @tmp=split/\t+/;
		push@inputs,\@tmp;
	}
	$/="》";
	my $titles=basename$input;
	open HH,">$outHtml.tmp" or die "$!";
	print HH <<HTML;
	<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
	<html lang="zh_CN" xmlns="http://www.w3.org/1999/xhtml"><head><title>$title</title>
	<meta charset="UTF-8"></meta>
	<!--[if lt IE 9]>
	<script src="$srcPath/js/html5shiv.min.js"></script><script src="$srcPath/js/respond.min.js"></script>
	<![endif]-->
	<meta content="IE=edge" http-equiv="X-UA-Compatible"></meta>
	<meta content="width=device-width, initial-scale=1" name="viewport"></meta>
	<link href="$srcPath/css/bootstrap.min.css" type="text/css" rel="stylesheet" />
	<link href="$srcPath/css/index.css" type="text/css" rel="stylesheet" />
	<link href="$srcPath/js/fancyBox/jquery.fancybox.css" type="text/css" rel="stylesheet" />
	<link href="$srcPath/css/nav.css" type="text/css" rel="stylesheet" />
	<link href="$srcPath/css/raxus.css" type="text/css" rel="stylesheet" />
	<script src="$srcPath/js/jquery-1.11.3.min.js" type="text/javascript"></script>
	<script src="$srcPath/js/nav.js" type="text/javascript"></script>
	<script src="$srcPath/js/raxus-slider.min.js" type="text/javascript"></script>
	<script src="$srcPath/js/fancyBox/jquery.fancybox.pack.js" type="text/javascript"></script>
	<script src="$srcPath/js/fancyBox/jquery.mousewheel-3.0.6.pack.js" type="text/javascript"></script>
	<script src="$srcPath/js/bootstrap.min.js" type="text/javascript"></script>
	<script src="$srcPath/js/ready.js" type="text/javascript"></script>
	<script src="$srcPath/js/scrolltop.js" type="text/javascript"></script>
	</head>
	<body>
	<div class="container shadow"><header><img src="$srcPath/images/logo.jpg" class="pull-right" />
	<div role="main" ><header><h2 id="title" class="text-center">$title</h2>
	</header>
	</div>
HTML
		if($text){
			print HH "<div class=\"table-responsive\"><p>$text</p></div>\n";
		}
		print HH "<div class=\"table-responsive\"><table class=\"table table-bordered table-hover table-striped\"><thead><tr class=\"bg-info\">\n";
		for (my $i=0;$i<=$#{$inputs[0]};$i++){
			print HH "<th>$inputs[0][$i]</th>\n";	
		}
		print HH "</tr></thead>\n<tbody>\n";
		for (my $k=1;$k<=$#inputs ;$k++) {
			print HH "<tr>";
			for (my $i=0;$i<=$#{$inputs[$k]};$i++){
				if($linkColNum){
					my $j=$i+1;
					if($linkColNum=~/,$j,/){
						#print "out:$outHtml\n$k  $i $inputs[$k][$i]\n"  if(!exists $$linkHash{$inputs[$k][$i]});exit;
						print HH "<td><a href=\"$$linkHash{$titles}{$inputs[$k][$i]}\" target=\"_blank\">$inputs[$k][$i]</a></td>";
					}
					else{
						print HH "<td>$inputs[$k][$i]</td>";
					}
				}
				else{
					print HH "<td>$inputs[$k][$i]</td>";
				}
			}
			print HH "</tr>\n";
		}	
print HH <<XGL;
	</tbody>
	</table>
	</body>
	</html>
XGL
	close HH;
	`iconv -f "GB2312" -t "UTF-8" $outHtml.tmp -o $outHtml `;
	`rm -r  $outHtml.tmp `;
}







sub USAGE{
	my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
	Options:
	-oldxml	<file>	input the xml file produced from the regular analysis, forced
	-o	<str>	output xml file, default the path same with old xml

	-h	Help

Example:

USAGE
	print $usage;
	exit;
}


