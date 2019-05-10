#!/usr/bin/perl
#use warnings;
use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);
use File::Basename qw(basename dirname);
use List::Util qw/max min/;
use XML::Writer;
use IO::File;
use Encode;
my ($data_cfg, $detail_cfg, $output, $id, $cloud);
GetOptions(
        "h|?"           =>\&USAGE,
	"id:s"		=>\$id,
	"cfg1:s"	=>\$data_cfg,
	"cfg2:s"	=>\$detail_cfg,
        "o:s"           =>\$output,
	"cloud"		=>\$cloud,
)or &USAGE;
&USAGE unless ($id and $detail_cfg and $data_cfg);

if(!defined $cloud){
	$output = "$id/configtest_local.xml";
}else{
	$output = "$id/configtest_raw.xml";
}

$id=abs_path($id);
$output=abs_path($output);
my $outpath=dirname($output);
`mkdir -p $outpath/BMK_9_html/BMK_1_template `       unless(-d "$outpath/BMK_9_html/BMK_1_template");
`cp $Bin/template/* $outpath/BMK_9_html/BMK_1_template`;
`rm $outpath/BMK_9_html/BMK_1_template/*.txt`;
my %config=&readConfig($detail_cfg);
my %sample=();
my %info=();	
&extractInfo;
#&extractInfo;
my ($pid,$tid)=(0,0);
&produceMap("$id/BMK_4_geneExpression/BMK_1_Mapped_Statistics","$outpath/BMK_9_html/BMK_1_template/Mappedstat.xls");
	
&run_or_die("perl $Bin/build_qc_xml.pl -id $id")        unless(-e "$id/BMK_9_html/Assess/assess.html");
if(!defined $cloud){
	&run_or_die("perl $Bin/build_deg_xml.pl -id $id -cfg $detail_cfg");
}else{
	&run_or_die("perl $Bin/build_deg_xml.pl -id $id -cfg $detail_cfg -cloud");
}


my $report_time=&gaintime();
my $writer = XML::Writer->new(OUTPUT => 'self');
$writer->xmlDecl('UTF-8');
$writer->startTag('report');
$writer->emptyTag('report_version','value'=>$config{report_version});
my $report_name;
$report_name=$config{Project_name}	if($config{Project_name});
$report_name ||= "转录组分析";
$writer->emptyTag('report_name','value'=>$report_name);
$writer->emptyTag('report_code','value'=>$config{report_code});
$writer->emptyTag('report_user','value'=>$config{report_user});
$writer->emptyTag('report_user_addr','value'=>$config{report_user_addr});
$writer->emptyTag('report_time','value'=>"XXXX/XX/XX;XXXX/XX/XX;XXXX/XX/XX;$report_time");
$writer->emptyTag('report_abstract','value'=>"");
###############################摘要
$writer->emptyTag('h1','name'=>"摘要",'type'=>'type1','desc'=>"");
$writer->emptyTag('h2','desc'=>'合同关键指标','type'=>'type1','name'=>'合同关键指标');
$writer->emptyTag('p','desc'=>"完成$info{Sample} 个样品的转录组测序，每个样品测序产出不少于 6 Gb Clean Data，Q30 碱基百分比达到 85%。完成基因表达量分析和差异表达分析。完成可变剪接预测分析。完成新基因预测及功能注释分析。完成差异表达基因功能注释分析。",'type'=>"type1");
$writer->emptyTag('h2','desc'=>'分析结果概述','type'=>'type1','name'=>'分析结果概述');
$writer->emptyTag('p','desc'=>"完成$info{Sample} 个样品的转录组测序，共获得$info{Total} Gb Clean Data，各样品Clean Data均达到$info{Min} Gb及以上，Q30碱基百分比在$info{Q30} %及以上。分别将各样品的 Clean Reads 与指定的参考基因组进行序列比对，比对效率从 $info{Map_min}% 到 $info{Map_max}% 不等。基于比对结果，进行可变剪接预测分析、基因结构优化分析以及新基因的发掘，发掘新基因 $info{newGene} 个，其中 $info{newAnno} 个得到功能注释。基于比对结果，进行基因表达量分析。根据基因在不同样品中的表达量识别差异表达基因，并对其进行功能注释和富集分析。",'type'=>"type1");

########################实验流程
$writer->emptyTag('h1','name'=>"1 实验流程",'type'=>'type1','desc'=>"1 实验流程");
$writer->emptyTag('p','desc'=>"转录组测序实验流程包括样品检测、文库构建及其质量控制和上机测序。实验流程见下图：",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. 转录组测序实验流程图",'type'=>"type1",'path'=>"BMK_9_html/BMK_1_template/P01_RNA-Seq_experimental_workflow_new.png");
$writer->emptyTag('h2','name'=>"1.1 样品检测",'type'=>'type1','desc'=>"1.1 样品检测");
$writer->emptyTag('p','desc'=>"高质量的RNA是整个项目成功的基础，为保证得到的数据准确性，我们使用以下方法对样品进行检测，检测结果达到要求后方可进行建库：",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(1) Nanodrop检测RNA的纯度（OD260/280）、浓度、核酸吸收峰是否正常；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(2) Agilent 2100精确检测RNA的完整性，检测指标包括：RIN值、28S/18S、图谱基线有无上抬、5S峰。",'type'=>"type1");
$writer->emptyTag('h2','name'=>"1.2 文库构建",'type'=>'type1','desc'=>"1.2 文库构建");
$writer->emptyTag('p','desc'=>"样品检测合格后，进行文库构建，主要流程如下：",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(1) 用带有Oligo（dT）的磁珠富集真核生物mRNA；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(2) 加入Fragmentation Buffer将mRNA进行随机打断；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(3) 以mRNA为模板，用六碱基随机引物（random hexamers）合成第一条cDNA链，然后加入缓冲液、dNTPs、RNase H和DNA polymerase I合成第二条cDNA链，利用AMPure XP beads纯化cDNA；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(4) 纯化的双链cDNA再进行末端修复、加A尾并连接测序接头，然后用AMPure XP beads进行片段大小选择；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(5) 最后通过PCR富集得到cDNA文库。",'type'=>"type1");
$writer->emptyTag('h2','name'=>"1.3 文库质控",'type'=>'type1','desc'=>"1.3 文库质控");
$writer->emptyTag('p','desc'=>"文库构建完成后，对文库质量进行检测，检测结果达到要求后方可进行上机测序，检测方法如下：",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(1) 使用Qubit2.0进行初步定量，使用Agilent 2100对文库的insert size进行检测，insert size符合预期后才可进行下一步实验。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(2) Q-PCR方法对文库的有效浓度进行准确定量（文库有效浓度＞2nM），完成库检。",'type'=>"type1");
$writer->emptyTag('h2','name'=>"1.4 上机测序",'type'=>'type1','desc'=>"1.4 上机测序");
$writer->emptyTag('p','desc'=>"库检合格后，不同文库按照目标下机数据量进行pooling，用IlluminaHiSeq平台进行测序。",'type'=>"type1");
###########
$writer->emptyTag('h1','name'=>"2 生物信息学分析",'type'=>'type1','desc'=>"2 生物信息学分析");
$writer->emptyTag('p','desc'=>"将下机数据进行过滤得到Clean Data，与指定的参考基因组进行序列比对，得到的Mapped Data，进行插入片段长度检验、随机性检验等文库质量评估；进行可变剪接分析、新基因发掘和基因结构优化等结构水平分析；根据基因在不同样品或不同样品组中的表达量进行差异表达分析、差异表达基因功能注释和功能富集等表达水平分析。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"转录组生物信息分析流程见下图：",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. 转录组生物信息分析流程图",'type'=>"type1",'path'=>"BMK_9_html/BMK_1_template/P02_RNA-Seq_analysis_workflow.png");

$writer->emptyTag('h2','name'=>"2.1 测序数据及其质量控制",'type'=>'type1','desc'=>"2.1 测序数据及其质量控制");
$writer->emptyTag('p','desc'=>"该项目各样品数据产出统计见下表：",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>'注：Customer_ID：客户样品编号；SampleID：百迈客样品分析编号；ReadSum：Clean Data中pair-end Reads总数；BaseSum：Clean Data总碱基数；GC(%)：Clean Data GC 含量，即Clean Data中G和C两种碱基占总碱基的百分比；N(%)：Clean Data中N碱基含量；Q20(%)：Clean Data质量值大于或等于20的碱基所占的百分比；Q30(%)：Clean Data质量值大于或等于30的碱基所占的百分比。','type'=>"full",'name'=>"表$tid. 测序数据统计表",'path'=>"BMK_9_html/BMK_1_template/AllSample_GC_Q.stat");
#$writer->startTag('file_list','name'=>"",'desc'=>"",'type'=>"xls");
$writer->emptyTag('p','desc'=>"数据质控网页",'type'=>"type1");
if (!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"assess.html",'action'=>"xls",'path'=>"BMK_9_html/Assess/assess.html",'type'=>"xls");
}else{
	$writer->emptyTag('file','desc'=>"",'name'=>"assess.xml",'action'=>"xls",'path'=>"BMK_9_html/assess.xml",'type'=>"xml");
}
#$writer->endTag('file_list');

$writer->emptyTag('h2','name'=>"2.2 转录组数据与参考基因组序列比对",'type'=>'type1','desc'=>"2.2 转录组数据与参考基因组序列比对");
$writer->emptyTag('p','desc'=>"本项目使用指定的基因组作为参考进行序列比对及后续分析。参考基因组下载地址见：$config{Download}",'type'=>"type1");
$writer->emptyTag('p','desc'=>"HISAT2 <a href=\"#ref1\">[1]</a>是一个来自 RNA 测序实验 reads 的高效比对系统，TopHat2 / Bowtie2 的继任者。HISAT 使用一个基于 Burrows-Wheeler 变换和 Ferragina-Manzini（FM）索引的索引方案，利用了两类索引进行比对：一个全基因组FM索引以定位每个比对，和许多局部 FM 索引用于非常快地扩展这些比对，实现了更快的速度和更少的资源占用。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"利用 StringTie <a href=\"#ref2\">[2]</a> 对比对上的 reads 进行组装，StringTie 是基于最优化理论建立的算法，利用比对信息构建多可变剪切图谱运用构建流量网络从而根据最大流量算法来对 reads 进行组装和评估其表达量，同 Cufflinks 等其他软件相比可以构建更完整的转录本和更好的评估表达量。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"比对分析完成后利用 StringTie 对比对上的 Reads 进行组装和定量，分析流程如下图：",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. HISA2T 分析流程",'type'=>"type1",'path'=>"BMK_9_html/BMK_1_template/hisat2.png");
$writer->emptyTag('h3','name'=>"2.2.1 比对效率统计",'type'=>'type1','desc'=>"2.2.1 比对效率统计");
$writer->emptyTag('p','desc'=>"比对效率指Mapped Reads占Clean Reads的百分比，是转录组数据利用率的最直接体现。比对效率除了受数据测序质量影响外，还与指定的参考基因组组装的优劣、参考基因组与测序样品的生物学分类关系远近（亚种）有关。通过比对效率，可以评估所选参考基因组组装是否能满足信息分析的需求。",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>"注：BMK-ID：百迈客样品分析编号；Total Reads：Clean Reads数目，按单端计；Mapped Reads：比对到参考基因组上的Reads数目及在Clean Reads中占的百分比；Uniq Mapped Reads：比对到参考基因组唯一位置的Reads数目及在Clean Reads中占的百分比；Multiple Map Reads：比对到参考基因组多处位置的Reads数目及在Clean Reads中占的百分比；Reads Map to '+'：比对到参考基因组正链的Reads数目及在Clean Reads中占的百分比；Reads Map to '-'：比对到参考基因组负链的Reads数目及在Clean Reads中占的百分比。",'type'=>"full",'name'=>"表$tid. 样品测序数据与所选参考基因组的序列比对结果统计表",'path'=>"BMK_9_html/BMK_1_template/Mappedstat.xls");
$writer->emptyTag('h3','name'=>"2.2.2 比对结果作图",'type'=>'type1','desc'=>"2.2.2 比对结果作图");
$writer->emptyTag('p','desc'=>"将比对到不同染色体上的Reads进行位置分布统计，绘制Mapped Reads在所选参考基因组上的覆盖深度分布图。",'type'=>"type1");
&piclist("Mapped Reads在参考基因组上的位置及覆盖深度分布图","注：横坐标为染色体位置；纵坐标为覆盖深度以2为底的对数值，以10kb作为区间单位长度，划分染色体成多个小窗口（Window），统计落在各个窗口内的Mapped Reads作为其覆盖深度。蓝色为正链，绿色为负链。","$id/BMK_4_geneExpression/BMK_1_Mapped_Statistics/*.coverage.png","BMK_4_geneExpression");

$writer->emptyTag('p','desc'=>"统计Mapped Reads在指定的参考基因组不同区域（外显子、内含子和基因间区）的数目，绘制基因组不同区域上各样品Mapped Reads的分布图，如下：",'type'=>"type1");
&piclist("基因组不同区域Reads分布图","注：图中将基因组分为外显子区、基因间区、内含子区，区域大小按Map到相应区域的Reads在所有Mapped Reads中所占的百分比。","$id/BMK_4_geneExpression/BMK_1_Mapped_Statistics/*.type.png","BMK_4_geneExpression");
$writer->emptyTag('h3','name'=>"2.2.3 比对结果可视化",'type'=>'type1','desc'=>"2.2.3 比对结果可视化");
$writer->emptyTag('p','desc'=>"转录组测序Reads与参考基因组序列比对结果文件（通常为BAM格式）、物种参考基因组序列和注释文件，推荐使用整合基因组浏览器（IGV <a href=\"#ref3\">[3]</a>，Integrative Genomics Viewer）进行可视化浏览。IGV 具有以下特点：",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(1) 能在不同尺度下显示单个或多个Reads在参考基因组上的位置，包括Reads在各个染色体上的分布情况和在注释的外显子、内含子、剪接接合区、基因间区的分布情况等；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(2) 能在不同尺度下显示不同区域的Reads丰度，以反映不同区域的转录水平；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(3) 能显示基因及其剪接异构体的注释信息；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(4) 能显示其他注释信息；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(5) 既可以从远程服务器端下载各种注释信息，又可以从本地加载注释信息。",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. IGV浏览器界面",'type'=>"type1",'path'=>"BMK_9_html/BMK_1_template/P07_IGV_interface.png");
############
$writer->emptyTag('h2','name'=>"2.3 转录组文库质量评估",'type'=>'type1','desc'=>"2.3 转录组文库质量评估");
$writer->emptyTag('p','desc'=>"合格的转录组文库是转录组测序的必要条件，为确保文库的质量，从以下3个不同角度对转录组测序文库进行质量评估：",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(1) 通过检验插入片段在基因上的分布，评估mRNA片段化的随机性、mRNA的降解情况；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(2) 通过插入片段的长度分布，评估插入片段长度的离散程度；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(3) 通过绘制饱和度图，评估文库容量和Mapped Data是否充足。",'type'=>"type1");

$writer->emptyTag('h3','name'=>"2.3.1 mRNA片段化随机性检验",'type'=>'type1','desc'=>"2.3.1 mRNA片段化随机性检验");
$writer->emptyTag('p','desc'=>"mRNA片段化后的插入片段大小选择，是从mRNA序列中独立随机地抽取子序列，mRNA数目越大、打断方式和时间控制得越合适，目的RNA每个部分被抽取到的可能性就越接近，mRNA片段化随机性越高，mRNA上覆盖的Reads越均匀。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"通过Mapped Reads在各mRNA转录本上的位置分布，模拟mRNA片段化结果，检验mRNA片段化的随机程度。如果mRNA存在严重降解，被降解的碱基序列不能被测序，即无Reads比对上。因此，通过查看Mapped Reads在mRNA转录本上的位置分布可了解mRNA的降解情况。样品Mapped Reads在mRNA转录本上的位置分布如下图：",'type'=>"type1");
&piclist("Mapped Reads在mRNA上的位置分布图","注：横坐标为标准化后的mRNA位置，纵坐标为对应位置区间内Reads在总Mapped Reads中所占百分比。由于参考的mRNA长度不同，作图时对把每个mRNA按照长度划分成100个区间，进而统计每一区间内的Mapped Reads数目及所占的比例，图中反映的是所有mRNA各个区间内的Mapped Reads比例的汇总。","$id/BMK_4_geneExpression/BMK_2_Library_Assessment/*.randcheck.png","BMK_4_geneExpression");
$writer->emptyTag('h3','name'=>"2.3.2 插入片段长度检验",'type'=>'type1','desc'=>"2.3.2 插入片段长度检验");
$writer->emptyTag('p','desc'=>"插入片段长度检验插入片段长度的离散程度能直接反映出文库制备过程中磁珠纯化的效果。通过插入片段两端的Reads在参考基因组上的比对起止点之间的距离计算插入片段长度。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"大部分的真核生物基因为断裂基因，外显子被内含子隔断，而转录组测序得到的是无内含子的成熟mRNA。当mRNA中跨内含子的片段两端的Reads比对到基因组上时，比对起止点之间的距离要大于插入片段长度。因此，在插入片段长度模拟分布图中，主峰右侧形成1个或多个杂峰。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"各样品的插入片段长度模拟分布图如下：",'type'=>"type1");
&piclist("插入片段长度模拟分布图","注：横坐标为双端Reads在参考基因组上的比对起止点之间的距离，范围为0到800bp；纵坐标为比对起止点之间不同距离的双端Reads或插入片段数量。","$id/BMK_4_geneExpression/BMK_2_Library_Assessment/*.insertSize.png","BMK_4_geneExpression");
$writer->emptyTag('h3','name'=>"2.3.3 转录组测序数据饱和度检验",'type'=>'type1','desc'=>"2.3.3 转录组测序数据饱和度检验");
$writer->emptyTag('p','desc'=>"为了评估数据是否充足并满足后续分析，对测序得到的基因数进行饱和度检测。由于一个物种的基因数目是有限的，且基因转录具有时间和空间特异性，因此随着测序量的增加，检测到的基因数目会趋于饱和。对于表达量越高的基因，越容易被检测定量。因此，对于表达量越低的基因，需要更大的数据量才能被准确定量。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"使用各样品的Mapped Data对检测到的不同表达情况的基因数目饱和情况进行模拟，绘制曲线图如下，可查看随着测序数据量的增加，检测到的不同表达量的基因数目是否趋于饱和。",'type'=>"type1");
&piclist("转录组数据饱和度模拟图","注：本图为随机抽取10%、20%、30%……90%的总体测序数据单独进行基因定量分析的结果；横坐标代表抽取数据定位到基因组上的Reads数占总定位的reads数的百分比，纵坐标代表所有抽样结果中表达量差距小于15%的Gene在各个FPKM范围的百分比。","$id/BMK_4_geneExpression/BMK_2_Library_Assessment/*.Saturation.png","BMK_4_geneExpression");
$writer->emptyTag('h2','name'=>"2.4 SNP/InDel分析",'type'=>'type1','desc'=>"2.4 SNP/InDel分析");
$writer->emptyTag('p','desc'=>"SNP (Single Nucleotide Polymorphisms) 是指在基因组上由单个核苷酸变异形成的遗传标记，其数量很多，多态性丰富。百迈客基于各样品 reads 与参考基因组序列的 Hisat 比对结果，使用 GATK <a href=\"#ref4\">[4]</a> 软件识别测序样品与参考基因组间的单碱基错配，识别潜在的 SNP 位点。进而可以分析这些 SNP 位点是否影响了基因的表达水平或者蛋白产物的种类。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"InDel (insertion-deletion) 是指相对于参考基因组，样本中发生的小片段的插入缺失，该插入缺失可能含一个或多个碱基。GATK 也能够检测样品的插入缺失（InDel）。InDel 变异一般比SNP变异少，同样反映了样品与参考基因组之间的差异，并且编码区的 InDel 会引起移码突变，导致基因功能上的变化。GATK 识别标准如下：",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(1) 35bp范围内连续出现的单碱基错配不超过3个； (2) 经过序列深度标准化的 SNP 质量值大于 2.0。 各样品分别按照以上条件筛选，最终获得可靠的 SNP 位点。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"SnpEff <a href=\"#ref5\">[5]</a> 是一款用于注释变异（SNP、InDel）和预测变异影响的软件。根据变异位点在参考基因组上的位置以及参考基因组上的基因位置信息，可以得到变异位点在基因组发生的区域（基因间区、基因区或 CDS 区等），以及变异产生的影响（同义非同义突变等）。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"由于转录完成之后，mRNA 除了需要加帽、加 Ploy(A) 和可变剪接之外，较少 mRNA 会经历 RNA 编辑（RNA editing），从而会产生单碱基的替换、插入、缺失。RNA编辑能使同一基因产生序列多样的 mRNA，但是这种多态性不是基因组固有的多态性。从比对结果来看，SNP 和单碱基替换的 RNA 编辑结果是一样的。因此，通过转录组测序数据识别出 SNP 不免会含有RNA编辑的产物。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"SNP/InDel位点信息:",'type'=>"type1");
&filelist("","$id/BMK_7_SNP_Analysis/final.*.anno.gatk.all.list","BMK_7_SNP_Analysis");
$tid++;
$writer->emptyTag('table','desc'=>'注：Chr：SNP/InDel位点所在染色体编号；Pos：SNP/InDel位点在染色体上的位置；Gene_id：SNP/InDel位点所在的基因或原来未注释的基因区（表中用Intergenic表示）；Ref：所选参考基因组中的SNP/InDel等位；Alt：测序样品中识别到的其他的SNP/InDel等位；T*：样品T*该SNP/InDel位点的分型；Depth：样品T*该SNP/InDel位点的测序深度；AlleDp：样品T*该SNP/InDel位点的各等位测序深度；Effect：SNP/InDel所在区域或类型；Codon_change：编码改变方式，未改变用点表示。核酸编码表见附表3，Effect具体说明详见：http://snpeff.sourceforge.net/SnpEff_manual.html。','type'=>"full",'name'=>"表$tid. 示例: SNP/InDel位点信息",'path'=>"BMK_9_html/BMK_1_template/final.InDel.anno.gatk.all.xls");
$writer->emptyTag('h3','name'=>"2.4.1 SNP位点统计",'type'=>'type1','desc'=>"2.4.1 SNP位点统计");
$writer->emptyTag('p','desc'=>"根据 SNP 位点碱基替换的不同方式，可以将 SNP 位点分为转换（Transition）和颠换（Transversion）两种类型。根据 SNP 位点的等位（Allele）数目，可以将 SNP 位点分为纯合型 SNP 位点（只有一个等位）和杂合型 SNP 位点（两个或多个等位）。不同物种杂合型 SNP 所占的比例存在差异。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"对各样品筛选出的 SNP 位点数目、转换类型比例、颠换类型比例以及杂合型 SNP 位点比例进行统计，如下表：",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>'注：BMK-ID：百迈客样品分析编号；SNP Number：SNP位点总数；Genic SNP：基因区SNP位点总数；Intergenic SNP：基因间区SNP位点总数；Transition：转换类型的SNP位点数目在总SNP位点数目中所占的百分比；Transversion：颠换类型的SNP位点数目在总SNP位点数目中所占的百分比；Heterozygosity：杂合型SNP位点数目在总SNP位点数目中所占的百分比。','type'=>"full",'name'=>"表$tid. SNP位点统计表",'path'=>"BMK_7_SNP_Analysis/AllSample.SNP.stat");
$writer->emptyTag('p','desc'=>"SNP突变类型统计分布如下图所示：",'type'=>"type1");
&piclist("SNP突变类型分布图","注：横轴为SNP突变类型，纵轴为相应的SNP数目。","$id/BMK_7_SNP_Analysis/BMK_3_SNP_type/*.SNP.type.png","BMK_7_SNP_Analysis");
$writer->emptyTag('h3','name'=>"2.4.2 基因的SNP密度分布",'type'=>'type1','desc'=>"2.4.2 基因的SNP密度分布");
$writer->emptyTag('p','desc'=>"将每个基因的SNP位点数目除以基因的长度，得到每个基因的SNP位点密度值，统计所有基因的SNP位点密度值并做密度分布图。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"基因的SNP位点密度分布图如下：",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"注：横轴为基因上平均每1000bp序列中分布的SNP数目，纵轴为基因数。",'name'=>"图$pid. SNP密度分布图",'type'=>"type1",'path'=>"BMK_7_SNP_Analysis/AllSample.SNP_density.png");

$writer->emptyTag('h3','name'=>"2.4.3 SNP/InDel注释",'type'=>'type1','desc'=>"2.4.3 SNP/InDel注释");
$writer->emptyTag('p','desc'=>"采用SNPEff分别对SNP，InDel注释，SNP，InDel的注释结果统计如下所示：",'type'=>"type1");
&piclist("SNP注释分类图","注：纵轴为SNP所在区域或类型，横轴为分类数目。","$id/BMK_7_SNP_Analysis/BMK_2_SNP_anno/*.SNP.anno.stat.png","BMK_7_SNP_Analysis");
&piclist("InDel注释分类图","注：纵轴为InDel所在区域或类型，横轴为分类数目。","$id/BMK_7_SNP_Analysis/BMK_1_InDel_anno/*.InDel.anno.stat.png","BMK_7_SNP_Analysis");
$writer->emptyTag('h2','name'=>"2.5 可变剪接事件预测",'type'=>'type1','desc'=>"2.5 可变剪接事件预测");
$writer->emptyTag('p','desc'=>"基因转录生成的前体mRNA（pre-mRNA），有多种剪接方式，选择不同的外显子，产生不同的成熟mRNA，从而翻译为不同的蛋白质，构成生物性状的多样性。这种转录后的mRNA加工过程称为可变剪接或选择性剪接（Alternative splicing）。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"采用 StringTie 对 Hisat2 的比对结果进行拼接，通过ASprofile <a href=\"#ref6\">[6]</a> 软件获取每个样品存在的可变剪接类型及相应表达量。基因可变剪接类型如下图所示：",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"注：(A) 外显子跳跃和多外显子跳跃；(B) 单内含子保留和多内含子保留；(C) 可变外显子；(D) 可变转录起始位点；(E) 可变转录终止位点；其中红色处为可变剪接类型。",'name'=>"图$pid. 基因可变剪接类型",'type'=>"type1",'path'=>"BMK_9_html/BMK_1_template/asprofile.png");
$writer->emptyTag('p','desc'=>"ASprofile 软件将可变剪接类型细分为12类，分别为：",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(1) TSS: Alternative 5' first exon (transcription start site) 第一个外显子可变剪切；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(2) TTS: Alternative 3' last exon (transcription terminal site) 最后一个外显子可变剪切；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(3) SKIP: Skipped exon(SKIP_ON,SKIP_OFF pair) 单外显子跳跃；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(4) XSKIP: Approximate SKIP (XSKIP_ON,XSKIP_OFF pair) 单外显子跳跃（模糊边界）；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(5) MSKIP: Multi-exon SKIP (MSKIP_ON,MSKIP_OFF pair) 多外显子跳跃；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(6) XMSKIP: Approximate MSKIP (XMSKIP_ON,XMSKIP_OFF pair) 多外显子跳跃（模糊边界）；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(7) IR: Intron retention (IR_ON, IR_OFF pair) 单内含子滞留；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(8) XIR: Approximate IR (XIR_ON,XIR_OFF pair) 单内含子滞留（模糊边界）；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(9) MIR: Multi-IR (MIR_ON, MIR_OFF pair) 多内含子滞留 ；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(10) XMIR: Approximate MIR (XMIR_ON, XMIR_OFF pair) 多内含子滞留（模糊边界）；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(11) AE: Alternative exon ends (5', 3', or both) 可变 5'或3'端剪切；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(12) XAE: Approximate AE 可变 5'或3'端剪切（模糊边界）。",'type'=>"type1");
$writer->emptyTag('h3','name'=>"2.5.1 可变剪接事件数量统计",'type'=>'type1','desc'=>"2.5.1 可变剪接事件数量统计");
$writer->emptyTag('p','desc'=>"各样品中预测的可变剪接事件数量统计见下图：",'type'=>"type1");
&piclist("可变剪接事件数量统计图","注：横轴为该种事件下可变剪切的数量，纵轴为可变剪切事件的分类缩写。","$id/BMK_6_Alt_splice/*png","BMK_6_Alt_splice");
$writer->emptyTag('h3','name'=>"2.5.2 可变剪切事件结构统计",'type'=>'type1','desc'=>"2.5.2 可变剪切事件结构统计");
$writer->emptyTag('p','desc'=>"可变剪切事件结构统计表",'type'=>"type1");
&filelist("","$id/BMK_6_Alt_splice/*xls","BMK_6_Alt_splice");
$tid++;
$writer->emptyTag('table','desc'=>'注：event_id: AS事件编号；event_type: AS事件类型；gene_id: 基因ID；chrom: 染色体编号；event_start: AS事件起始位置；event_end: AS事件结束位置；event_pattern: AS事件特征 ；strand: 基因正负链信息。','type'=>"full",'name'=>"表$tid. 示例: 可变剪切事件结构统计表",'path'=>"BMK_9_html/BMK_1_template/AS.xls");

$writer->emptyTag('h2','name'=>"2.6 基因结构优化分析",'type'=>'type1','desc'=>"2.6 基因结构优化分析");
$writer->emptyTag('p','desc'=>"由于使用的软件或数据本身的局限性，导致所选参考基因组的注释往往不够精确，这样就有必要对原有注释的基因结构进行优化。如果在原有基因边界之外的区域有连续的Mapped Reads支持，将基因的非翻译区（Untranslated Region，UTR）向上下游延伸，修正基因的边界。此项目对463个基因结构进行了优化，基因结构优化结果见下面文件：",'type'=>"type1");
$writer->emptyTag('file','desc'=>"注：GeneID：基因ID；Locus：基因座，格式为“染色体编号:起点坐标-终点坐标”；Strand：正负链；Site：优化的位置，3'或5'UTR；OriginalRegion：原来注释的第一个或最后一个外显子的起止坐标；OptimizedRegion：延伸之后的第一个或最后一个外显子的起止坐标。",'name'=>"基因结构优化表",'action'=>"xls",'path'=>"BMK_9_html/BMK_1_template/$config{Project_key}.geneStructure.optimize.xls",'type'=>"xls");
$writer->emptyTag('h2','name'=>"2.7 新基因分析",'type'=>'type1','desc'=>"2.7 新基因分析");
$writer->emptyTag('h3','name'=>"2.7.1 新基因发掘",'type'=>'type1','desc'=>"2.7.1 新基因发掘");
$writer->emptyTag('p','desc'=>"基于所选参考基因组序列，使用StringTie软件对Mapped Reads进行拼接，并与原有的基因组注释信息进行比较，寻找原来未被注释的转录区，发掘该物种的新转录本和新基因，从而补充和完善原有的基因组注释信息。过滤掉编码的肽链过短（少于50个氨基酸残基）或只包含单个外显子的序列，共发掘6,108个新基因。新基因的GFF格式文件见下面文件：",'type'=>"type1");
$writer->emptyTag('file','desc'=>"注：文件总共有9列。第1列：染色体号；第2列：注释信息的来源，Cufflinks软件；第3列：注释特征（Feature）类型；第4、5列：特征序列的起止位置；第6列：得分，数字，注释信息可能性的说明，“.”表示缺失值；第7列：特征序列所在的正负链；第8列：仅对注释类型为CDS有效，表示起始编码的位置，有效值为0、1、2，“.”表示缺失值；第9列：以多个键值对组成的注释信息描述。",'name'=>"新基因GFF文件",'action'=>"xls",'path'=>"BMK_9_html/BMK_1_template/$config{Project_key}.newGene_final.filtered.gff.xls",'type'=>"xls");
$writer->emptyTag('p','desc'=>"提供基因组注释补充信息的同时，也提供以FASTA格式存储的新基因序列。新基因序列的FASTA文件见下面文件：",'type'=>"type1");
$writer->emptyTag('file','desc'=>"",'name'=>"新基因序列文件",'action'=>"xls",'path'=>"BMK_9_html/BMK_1_template/$config{Project_key}.newGene.longest_transcript.fa.xls",'type'=>"xls");
$writer->emptyTag('h3','name'=>"2.7.2 新基因功能注释",'type'=>'type1','desc'=>"2.7.2 新基因功能注释");
$writer->emptyTag('p','desc'=>"使用BLAST <a href=\"#ref7\">[7]</a> 软件将发掘的新基因与NR <a href=\"#ref8\">[8]</a>，Swiss-Prot <a href=\"#ref9\">[9]</a>，GO <a href=\"#ref10\">[10]</a>，COG <a href=\"#ref11\">[11]</a>，KOG <a href=\"#ref12\">[12]</a>，Pfam <a href=\"#ref13\">[13]</a>，KEGG <a href=\"#ref14\">[14]</a> 数据库进行序列比对，使用KOBAS2.0 <a href=\"#ref15\">[15]</a> 得到新基因的KEGG Orthology结果，预测完新基因的氨基酸序列之后使用HMMER <a href=\"#ref16\">[16]</a> 软件与Pfam数据库比对，获得新基因的注释信息。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"最终得到各数据库注释的新基因数量统计见下表：",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>'注：Annotated databases：用于功能注释的数据库；New Gene Number：获得相应数据库注释信息的新基因数目。','type'=>"full",'name'=>"表$tid. 新基因功能注释结果统计",'path'=>"BMK_9_html/BMK_1_template/NewGene_anno.xls");

$writer->emptyTag('h2','name'=>"2.8 表达定量及差异分析",'type'=>'type1','desc'=>"2.8 表达定量及差异分析");
$writer->emptyTag('p','desc'=>"关于基因表达量的整体分布，差异表达筛选，及功能注释和富集分析详见以下网页:",'type'=>"type1");
#$writer->startTag('file_list','name'=>"基因差异表达及功能分析网页",'desc'=>"基因差异表达及功能分析网页",'type'=>"xls");
#$writer->emptyTag('p','desc'=>"基因差异表达及功能分析网页",'type'=>"type1");
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"deg.html",'action'=>"xls",'path'=>"BMK_9_html/Gene/deg.html",'type'=>"xls");
}else{
	$writer->emptyTag('file','desc'=>"",'name'=>"deg.xml",'action'=>"xls",'path'=>"BMK_9_html/deg.xml",'type'=>"xml");
}
#$writer->endTag('file_list');
$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. 差异及富集分析示意图",'type'=>"type1",'path'=>"BMK_9_html/BMK_1_template/1.geneDEG.png");
$writer->emptyTag('file','desc'=>"",'name'=>"基因表达量文件",'action'=>"xls",'path'=>"BMK_9_html/BMK_1_template/All_gene_fpkm.xls",'type'=>"xls");
my $h2=8;
my @deu=glob("$id/BMK_5_DEG_Analysis/*/BMK_4_DEXSeqReport/*html");
if(scalar(@deu)>0){
$h2++;
$writer->emptyTag('h2','name'=>"2.$h2 DEU分析",'type'=>'type1','desc'=>"2.$h2 DEU分析");
$writer->emptyTag('p','desc'=>"RNA-seq 除了基因水平的差异分析外，还可以进行 exon 水平的差异分析。对于有生物学重复的样本，采用DEXSeq <a href=\"#ref2\">[17]</a> 进行 DEU（differential exon usage）分析。将 FDR < 0.01 作为筛选外显子差异表达的标准。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"DEU分析结果",'type'=>"type1");

my @deu_xls=glob("$id/BMK_5_DEG_Analysis/*/BMK_4_DEXSeqReport/*DEU.final.xls");
$writer->startTag('file_list','name'=>"",'desc'=>"",'type'=>"xls");
foreach my $xls(@deu_xls){
	my $base=basename $xls;
	my $dir=(split(/BMK_5_DEG_Analysis/,$xls))[1];
	$writer->emptyTag('file','desc'=>"",'name'=>"$base",'action'=>"xls",'path'=>"BMK_5_DEG_Analysis/$dir",'type'=>"xls");
}
$writer->endTag('file_list');

$writer->emptyTag('p','desc'=>"差异外显子模式图见下图：",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"注：(A) Fitted expression. The plot represents the expression estimates from a call to testForDEU.Shown in red is the exon that showed significant differential exon usage；(B) Normalized counts. As in Figure A, with normalized count values of each exon in each；(C) Fitted splicing. The plot represents the estimated effects, as in Figure A, but after subtraction of overall changes in gene expression；(D) Transcripts. As in Figure A, but including the annotated transcript models。",'name'=>"图$pid. 差异外显子模式图",'type'=>"type1",'path'=>"BMK_9_html/BMK_1_template/DEU_modle.png");
$writer->emptyTag('p','desc'=>"外显子差异分析结果",'type'=>"type1");

my @deu_html=glob("$id/BMK_5_DEG_Analysis/*/BMK_4_DEXSeqReport/*html");
$writer->startTag('file_list','name'=>"",'desc'=>"",'type'=>"xls");
foreach my $xls(@deu_html){
        my $dir=(split(/BMK_5_DEG_Analysis/,$xls))[1];
	my $base=(split(/\//,(split(/\/BMK_3_/,$dir))[1]))[0];
        $writer->emptyTag('file','desc'=>"",'name'=>"$base.html",'action'=>"xls",'path'=>"BMK_5_DEG_Analysis/$dir",'type'=>"xls");
}
$writer->endTag('file_list');

}

if(-d "$id/BMK_10_Gene_Fusion"){
$h2++;
$writer->emptyTag('h2','name'=>"2.$h2 融合基因分析",'type'=>'type1','desc'=>"2.$h2 融合基因分析");
$writer->emptyTag('p','desc'=>"融合基因是指将两个或多个基因的编码区首尾相连，置于同一套调控序列(包括启动子、增强子、核糖体结合序列、终止子等)控制之下，构成的嵌合基因。融合基因的表达产物为融合蛋白。我们使用 Fusionmap <a href=\"#ref18\">[18]</a> 在转录组中研究基因融合事件。Fusionmap 首先通过比对到基因组和转录本中双末端(pairend)关系的序列寻找候选的基因融合，然后采用通过与nt等数据库比较，过滤掉假阳性结果。",'type'=>"type1");
&piclist("检测到的基因融合事件","注：红色的线代表同一染色体上发生的融合事件，绿色的线代表不同染色体上发生的融合事件。","$id/BMK_10_Gene_Fusion/png/*_circos.png","BMK_10_Gene_Fusion");
$writer->emptyTag('p','desc'=>"融合基因事件统计:",'type'=>"type1");
&filelist("","$id/BMK_10_Gene_Fusion/report/*_FusionReport.txt.xls","BMK_10_Gene_Fusion");
$writer->emptyTag('p','desc'=>"注：Fusion ID：基因融合事件编号；Unique Mapping Position：映射到唯一位置的read count数；Count：基因融合的read count总数；Gene：发生融合的基因；Strand：基因融合发生在正义链还是反义链；Chromosome：融合基因位于的染色体；Start：基因起始位点；End：基因终止位点；Filter：基因融合事件类型；InFamilyList代表基因对属于同一家族；InParalogueList代表基因对来自同旁系同源组；SameEnsembleGene代表基因对有一个super ensemble基因交集；InBlackList代表基因对包含线粒体和核糖体基因或者假基因。上述类型均被过滤，因此本列均为空。",'type'=>"type1");
}

$writer->startTag('ref_list',desc=>"参考文献",type=>"type1",name=>"参考文献");
&reference ($writer,[
        ["Kim D, Langmead B, Salzberg S L. HISAT: a fast spliced aligner with low memory requirements[J]. Nature methods, 2015, 12(4): 357-360.","https://www.nature.com/articles/nmeth.3317","1"],
        ["Pertea M, Pertea G M, Antonescu C M, et al. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads[J]. Nature biotechnology, 2015, 33(3): 290-295.","https://www.nature.com/articles/nbt.3122","2"],
        ["Helga Thorvaldsdóttir, James T. Robinson, Jill P. Mesirov.  Integrative Genomics Viewer (IGV): high-performance genomics data visualization and exploration. Briefings in Bioinformatics 14, 178-192 (2013).","https://academic.oup.com/bib/article/14/2/178/208453","3"],
        ["McKenna A, Hanna M, Banks E, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research. 2010, 20(9): 1297-1303.","https://genome.cshlp.org/content/20/9/1297.short","4"],
        ["Cingolani P, Platts A, Wang L L, et al. A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3. Fly. 2012, 6(2): 80-92.","https://www.tandfonline.com/doi/full/10.4161/fly.19695","5"],
        ["Florea L, Song L, Salzberg S L. Thousands of exon skipping events differentiate among splicing patterns in sixteen human tissues. F1000Research, 2013, 2:188.","https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3892928/","6"],
	["Altschul S F, Madden TL, Zhang J, et al. Gapped BLAST and PSI BLAST: A New Generation of Protein Database Search Programs. Nucleic Acids Research. 1997, 25(17): 3389-3402.","https://academic.oup.com/nar/article/25/17/3389/1061651","7"],
        ["Deng YY, Li JQ, Wu SF, Zhu YP, et al. Integrated nr Database in Protein Annotation System and Its Localization. Computer Engineering. 2006, 32(5):71-74.","http://en.cnki.com.cn/Article_en/CJFDTOTAL-JSJC200605025.htm","8"],
        ["Apweiler R, Bairoch A, Wu CH, et al. UniProt: the universal protein knowledgebase. Nucleic acids research. 2004, 32: D115-D119.","https://www.ncbi.nlm.nih.gov/pubmed/14681372","9"],
        ["Ashburner M, Ball C A, Blake J A, et al. Gene ontology: tool for the unification of biology. Nature genetics. 2000, 25(1): 25-29.","https://www.ncbi.nlm.nih.gov/pubmed/10802651","10"],
	["Tatusov RL, Galperin MY, Natale D A. The COG database: a tool for genome scale analysis of protein functions and evolution. Nucleic Acids Research. 2000, 28(1):33-36.","https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102395/","11"],
	["Koonin EV, Fedorova ND, Jackson JD, et al. A comprehensive evolutionary classification of proteins encoded in complete eukaryotic genomes. Genome biology. 2004, 5(2): R7.","https://www.ncbi.nlm.nih.gov/pmc/articles/PMC395751/","12"],
	["Finn RD, Bateman A, Clements J, et al. Pfam: the protein families database. Nucleic acids research. 2013: gkt1223.","https://www.ncbi.nlm.nih.gov/pubmed/24288371","13"],
	["Kanehisa M, Goto S, Kawashima S, Okuno Y, et al. The KEGG resource for deciphering the genome. Nucleic Acids Research. 2004, 32:D277-D280.","https://www.ncbi.nlm.nih.gov/pubmed/14681412","14"],
	["Xie C, Mao X, Huang J, et al. KOBAS 2.0: a web server for annotation and identification of enriched pathways and diseases. Nucleic acids research, 2011, 39(suppl 2): W316-W322.","https://academic.oup.com/nar/article/39/suppl_2/W316/2507217","15"],
	["Eddy S R. Profile hidden Markov models. Bioinformatics, 1998, 14(9): 755-763.","https://academic.oup.com/bioinformatics/article/14/9/755/259550","16"],
	["Simon Anders, Alejandro Reyes and Wolfgang Huber (2012): Detecting differential usage of exons from RNA-seq data. Genome Research (2012) 22:2008. doi:10.1101/gr.133744.111","https://www.ncbi.nlm.nih.gov/pubmed/22722343/","17"],
	["Ge H, Liu K, Juan T, et al. FusionMap: detecting fusion genes from next-generation sequencing data at base-pair resolution[J]. Bioinformatics, 2011, 27(14):1922-1928.","https://www.ncbi.nlm.nih.gov/pubmed/21593131","18"],
]);
$writer->endTag('ref_list');
#<a href=\"#ref1\">[1]</a><a href=\"#ref2\">[2]</a>

$writer->emptyTag('h1','desc'=>'附录','type'=>'type1','name'=>'附录');
$writer->emptyTag('h2','desc'=>'1 英文材料、方法与数据分析','type'=>'type1','name'=>'1 英文材料、方法与数据分析');
$writer->emptyTag('p','desc'=>"英文描述见以下附件:",'type'=>"type1");
$writer->emptyTag('file','desc'=>"",'name'=>"英文材料、方法与数据分析",'action'=>"xls",'path'=>"BMK_9_html/BMK_1_template/materials_and_methods.pdf",'type'=>"xls");

$writer->endTag('report');
open OUT,">:utf8", "$output";
my $xmlstr=&decorate($writer->to_string);
$xmlstr = Encode::decode("utf-8", $xmlstr);
print OUT $xmlstr;
close(OUT);
$writer->end();

if(!defined $cloud){
	print "/share/nas2/genome/biosoft/Python/2.7.8/bin/python $Bin/htmlConvert/xml2HtmlConverter.py -i $output -o $id\n";
	`/share/nas2/genome/biosoft/Python/2.7.8/bin/python $Bin/htmlConvert/xml2HtmlConverter.py -i $output -o $id`;
}

################################################################################################
#		                        sub function
################################################################################################
sub extractInfo{
	open(CFG,$data_cfg)||die $!;
	while(<CFG>){
		chomp; my @tmp=split(/\s+/,$_); 
		if($_=~/^Sample/) {$sample{$tmp[1]}++;}
	}
	close(CFG);	
	$info{Sample}=scalar(keys %sample);
	open(CFG,$detail_cfg)||die $!;
	while(<CFG>){
		chomp;my @tmp=split(/\s+/,$_);
		if(exists $sample{$tmp[0]}){	$sample{$tmp[0]}=$tmp[1];	}
	}
	close(CFG);
	open(DATA,"$id/BMK_1_rawdata/AllSample_GC_Q.stat")||die $!;
	open(OUT,">$outpath/BMK_9_html/BMK_1_template/AllSample_GC_Q.stat")||die $!;
	my $head=<DATA>;$head=~s/^#//;
	print OUT "#Customer_ID\t$head";
	my $Total=0;
	my $Min;
	my @group_min=();
	my $Q30=100;
	while(<DATA>){
		chomp;my @tmp=split;
		print OUT "$sample{$tmp[0]}\t$_\n";
		$Total+=$tmp[2];
		push @group_min, $tmp[2];
		$Q30=$tmp[-1]	if($tmp[-1]<$Q30);
	}
	close(OUT);
	close(DATA);
	$Min=min @group_min;
	$info{Min}=sprintf("%.2f",$Min/1000000000);
	$info{Total}=sprintf("%.2f",$Total/1000000000);
	$info{Q30}=$Q30;	
	my @stats=glob("$id/BMK_4_geneExpression/BMK_1_Mapped_Statistics/*.mappedStat.xls");
	open(MAP,">$outpath/BMK_9_html/BMK_1_template/Map_stat.xls")||die $!;
	print MAP "BMK-ID\tTotal Reads\tMapped Reads\tUniq Mapped Reads\tMultiple Map Reads\tReads Map to '+'\tReads Map to '-'\n";
	my %maps=();
	my @maprate=();
	foreach my $stat(@stats){
		open(STAT,$stat)||die $!;
		while(<STAT>){
			chomp;my @tmp=split(/\t/,$_);
			$maps{$tmp[0]}="$tmp[1]($tmp[2])";
			if($tmp[0] eq "mapped Reads"){
				$tmp[2]=~s/\%//;
				push @maprate,$tmp[2];
			}
		}
		close(STAT);
	}
	close(MAP);
	$info{Map_max}=max @maprate;
	$info{Map_min}=min @maprate;
	my $line=`grep '>' $id/BMK_2_NewGene/$config{Project_key}.newGene.longest_transcript.fa |wc -l`;
	chomp($line);
	$info{newGene}=(split(/\s+/,$line))[0];
	open(ANNO,"$id/BMK_2_NewGene/BMK_1_NewGene_Anno/Function_Annotation.stat.xls")||die $!;
	open(ANNOS,">$outpath/BMK_9_html/BMK_1_template/NewGene_anno.xls")||die $!;
	my @db=();my %num=();
	while(<ANNO>){
		chomp;	next if($_=~/^#/);	my @tmp=split(/\t/,$_);
		$tmp[0]=~s/_Annotation//;
		push @db,$tmp[0];	$num{$tmp[0]}=$tmp[1];
		$info{newAnno}=$tmp[1]	if($tmp[0]=~/^All/);
	}
	print ANNOS "Annotated databases\tNew Gene Number\n";
	foreach my $d(@db){	print ANNOS "$d\t$num{$d}\n";	}
	close(ANNOS);
	close(ANNO);
	##########
	#extract information for generating report
	#
	`head -n 6 $id/BMK_7_SNP_Analysis/final.InDel.anno.gatk.all.list > $outpath/BMK_9_html/BMK_1_template/final.InDel.anno.gatk.all.xls`;
	print "head -n 6 $id/BMK_7_SNP_Analysis/final.InDel.anno.gatk.all.list > $outpath/BMK_9_html/BMK_1_template/final.InDel.anno.gatk.all.xls\n";
	###rMATs diff AS
	my $AS=(glob("$id/BMK_6_Alt_splice/*.AS.list.xls"))[0];
	print "head -n 6 $AS > $outpath/BMK_9_html/BMK_1_template/AS.xls\n";
	`head -n 6 $AS > $outpath/BMK_9_html/BMK_1_template/AS.xls`;
	my @diff=glob("$id/BMK_5_DEG_Analysis/BMK_3_*_vs_*/");
	`head -n 6 $diff[0]/BMK_5_diff_AS_analysis/SE.MATS.JC.xls > $outpath/BMK_9_html/BMK_1_template/Diff_AS.xls`;
	print "head -n 6 $diff[0]/BMK_5_diff_AS_analysis/SE.MATS.JC.xls > $outpath/BMK_9_html/BMK_1_template/Diff_AS.xls\n";
	###TFBS analysis
	my @tfbs_xls = glob "$id/BMK_5_DEG_Analysis/BMK_4_TF_Analysis/BMK_1_TFBS_Analysis/*Genes_TFBS_predictRes.xls";
	if(@tfbs_xls>0){
        	`head -n 6 $tfbs_xls[0] > $outpath/BMK_9_html/BMK_1_template/tfbs.xls`;
	}
	###TF_activity analysis
	my @tf_activity = glob "$id/BMK_5_DEG_Analysis/BMK_4_TF_Analysis/BMK_2_TF_activity/TFs_*.xls";
	if(@tf_activity>0){
		for (my $i=0;$i<@tf_activity;$i++){
			my $tmp=$tf_activity[$i];
	                if ($tmp=~m/influence/) {
        	                `head -n 6 $tmp > $outpath/BMK_9_html/BMK_1_template/tfs_influence.xls`;
                	}elsif($tmp=~m/cornet/){
                        	`head -n 6 $tmp > $outpath/BMK_9_html/BMK_1_template/tfs_cornet.xls`;
                	}else{
                        	`head -n 6 $tmp > $outpath/BMK_9_html/BMK_1_template/tfs_activity_grn.xls`;
                	}
		}
	`cp $id/BMK_5_DEG_Analysis/BMK_4_TF_Analysis/BMK_2_TF_activity/TFs_influences_heatmap.png $outpath/BMK_9_html/BMK_1_template/`;
	`cp $id/BMK_5_DEG_Analysis/BMK_4_TF_Analysis/BMK_2_TF_activity/TFs_network.png $outpath/BMK_9_html/BMK_1_template/`;
	}
	#
	`head -n 6 $diff[0]/BMK_1_Statistics_Visualization/*.DEG_final.xls > $outpath/BMK_9_html/BMK_1_template/DEG.xls`;
	print "head -n 6 $diff[0]/BMK_1_Statistics_Visualization/*.DEG_final.xls > $outpath/BMK_9_html/BMK_1_template/DEG.xls\n";
	`cp $id/BMK_2_NewGene/$config{Project_key}.newGene_final.filtered.gff $outpath/BMK_9_html/BMK_1_template/$config{Project_key}.newGene_final.filtered.gff.xls`;
	`cp $id/BMK_2_NewGene/$config{Project_key}.newGene.longest_transcript.fa $outpath/BMK_9_html/BMK_1_template/$config{Project_key}.newGene.longest_transcript.fa.xls`;
	`cp $id/BMK_4_geneExpression/BMK_3_Expression_Statistics/All_gene_fpkm.list $outpath/BMK_9_html/BMK_1_template/All_gene_fpkm.xls`;
	`cp $id/BMK_8_Gene_Structure_Optimize/$config{Project_key}.geneStructure.optimize.xls $outpath/BMK_9_html/BMK_1_template/$config{Project_key}.geneStructure.optimize.xls`;
}

sub decorate{
	my $xmlstr=shift;
	$xmlstr=~s/(<[^>]*)\>/$1\>\n/mgo;
	return $xmlstr;
}

sub reference  {
    my ($writer,$list)=@_;
    my ($data,$value,$unit,$item);
    
    foreach $data (@$list)    {		
	$writer->emptyTag('ref','name'=>$$data[0],'link'=>$$data[1],'id'=>$$data[2]);
	#$writer->characters($$data[2]);
    }    
}

sub piclist{
	my ($name,$desc,$pics,$path)=@_;
	$pid++;
	$writer->startTag('pic_list','name'=>"图$pid. $name",'desc'=>"$desc",'type'=>"type1");
	my @images=glob("$pics");
	foreach my $s(@images){
		my $base=basename $s;
		my $dir=dirname $s;
		my $tmp=(split(/$path/,$dir))[1];
		my $new =$path.$tmp;
        	$writer->emptyTag('pic','desc'=>"",'name'=>"$base",'type'=>"type1",'path'=>"$new/$base");
	}
	$writer->endTag('pic_list');
}

sub filelist{
	my ($desc,$pics,$path)=@_;
	$writer->startTag('file_list','name'=>"",'desc'=>"$desc",'type'=>"xls");
	my @images=glob("$pics");
	foreach my $s(@images){
		my $base=basename $s;
		my $dir=dirname $s;
		my $tmp=(split(/$path/,$dir))[1];
		$path .=$tmp;
        	$writer->emptyTag('file','desc'=>"",'name'=>"$base",'action'=>"xls",'path'=>"$path/$base",'type'=>"xls");
	}
	$writer->endTag('file_list');
}

sub gaintime{
	my $timestamp=time(); 
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($timestamp); 
	my $y = $year + 1900; 
	my $m = $mon + 1; 
	$timestamp=sprintf("%4d-%02d-%02d",$y,$m,$mday);
	return $timestamp
}


sub readConfig{
        my $configFile=shift;
        my $d=Config::General->new(-ConfigFile => "$configFile");
        my %config=$d->getall;
        return %config;
}

sub run_or_die{
        my ($cmd) = @_ ;
        &show_log($cmd);
        my $flag = system($cmd) ;
        if ($flag != 0){
                &show_log("Error: command fail: $cmd");
                exit(1);
        }
        &show_log("done.");
        return ;
}
sub show_log{
        my ($txt) = @_ ;
        my $time = time();
        my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime($time);
        $wday = $yday = $isdst = 0;
        my $Time=sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
        print "$Time:\t$txt\n" ;
}
sub produceMap{
        my ($dir,$file)=@_;
        my @stat=glob("$dir/*.mappedStat.xls");
        my %info=();
        foreach my $s(@stat){
                my $base=basename $s;
                my $sam=(split(/\./,$base))[0];
                open(STAT,$s)||die $!;
                while(<STAT>){
                        chomp;my @tmp=split(/\t/,$_);
                        $info{$sam}{total}=&format_figure($tmp[1])                      if($_=~/^Total/);
                        $info{$sam}{map}=&format_figure($tmp[1])."($tmp[2])"            if($_=~/^mapped/);
                        $info{$sam}{Uniq}=&format_figure($tmp[1])."($tmp[2])"           if($_=~/^Uniq/);
                        $info{$sam}{Multiple}=&format_figure($tmp[1])."($tmp[2])"       if($_=~/^Multiple/);
                        $info{$sam}{Plus}=&format_figure($tmp[1])."($tmp[2])"           if($_=~/Plus/);
                        $info{$sam}{Minus}=&format_figure($tmp[1])."($tmp[2])"          if($_=~/Minus/);
                }
                close(STAT);
        }
        open(FILE,">$file")||die $!;
        print FILE "BMK-ID\tTotal Reads\tMapped Reads\tUniq Mapped Reads\tMultiple Map Reads\tReads Map to '+'\tReads Map to '-'\n";
        foreach my $sam(sort{$a cmp $b} keys %info){
                print FILE "$sam\t$info{$sam}{total}\t$info{$sam}{map}\t$info{$sam}{Uniq}\t$info{$sam}{Multiple}\t$info{$sam}{Plus}\t$info{$sam}{Minus}\n";
        }
        close(FILE);
}
sub Integer_Three_Digit{#
        my $interger = shift;
        $interger=~s/(?<=\d)(?=(\d\d\d)+$)/,/g;
        return $interger;
}
sub format_figure{#
        my $figure = shift;
        if (!defined $figure) {
                die;
        }
        if ($figure=~/\./) {
	        if ($figure == 100) {
        	    $figure = 100;
	        } else {
        	    $figure = sprintf("%.2f",$figure);
        	}
        }else{
                $figure = Integer_Three_Digit($figure);
        }
        return $figure;
}



sub USAGE{
        my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
        Options:
	-id	<path>	input images path, forced
        -o      <path>	output xml name, default inpath/configtest_raw.xml
	-cfg1	<file>	data_cfg
	-cfg2	<file>	detail_cfg

        -h      Help

Example:

USAGE
        print $usage;
        exit;
}

