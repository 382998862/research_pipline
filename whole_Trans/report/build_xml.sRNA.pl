#!/usr/bin/perl
#use warnings;
use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);
use File::Basename qw(basename dirname);
use XML::Writer;
use IO::File;
use Encode;
my ($detail_cfg,$data_cfg,$output,$inpath,$cloud);
GetOptions(
        "h|?"           =>\&USAGE,
	"i:s"		=>\$inpath,
	"cfg1:s"	=>\$data_cfg,
	"cfg2:s"	=>\$detail_cfg,
        "o:s"           =>\$output,
	"cloud"		=>\$cloud,
)or &USAGE;
&USAGE unless ($inpath and $data_cfg and $detail_cfg);

$output ||= "$inpath/configtest.xml";

$inpath=abs_path($inpath);
$output=abs_path($output);
my $outpath=dirname($output);

my %sample=();
&relate($data_cfg);
my $sam_num=scalar(keys %sample);
my %config=&readConfig($detail_cfg);
my %info=&readConfig("$inpath/HTML/template/stat.info");
my ($sRNA_dir,$sRNA_relative_dir);
$sRNA_dir = (glob "$inpath/BMK_*_miRNA")[0];
$sRNA_relative_dir = (split /\//,$sRNA_dir)[-1];
my $deg_flag=1;
if(($sam_num ==1 )|| (!exists $config{Sep} && !exists $config{Com})){
	$deg_flag=0;
}

my ($pid,$tid)=(0,0);
my $report_time=&gaintime();
my $writer = XML::Writer->new(OUTPUT => 'self');
$writer->xmlDecl('UTF-8');
$writer->startTag('report');
$writer->emptyTag('report_version','value'=>$config{report_version});
my $report_name;
$report_name=$config{Project_name}	if($config{Project_name});
$report_name ||= "全转录组联合分析";
$writer->emptyTag('report_name','value'=>$report_name);
$writer->emptyTag('report_code','value'=>$config{Project_id});
$writer->emptyTag('report_user','value'=>$config{report_user});
$writer->emptyTag('report_user_addr','value'=>$config{report_user_addr});
$writer->emptyTag('report_time','value'=>"XXXX/XX/XX;XXXX/XX/XX;XXXX/XX/XX;$report_time");
$writer->emptyTag('report_abstract','value'=>"");
###############################摘要
my @abstract;

push @abstract,"合同关键指标：";
push @abstract,"(1) 完成 $sam_num 个样品的small RNA测序，每个样品测序得到不少于 10 M的Clean Reads，各样品Q30≥85%。";
push @abstract,"(2) 鉴定已知miRNA和预测新miRNA。";
push @abstract,"(3) miRNA表达定量分析。";
push @abstract,"(4) miRNA靶基因预测。";
if($deg_flag==1){
	push @abstract,"(5) 筛选差异表达miRNA，差异表达miRNA靶基因的功能注释和富集分析。";
}
push @abstract,"分析结果概述：";
push @abstract,"(1) 完成了 $sam_num 个样品的small RNA测序，共得到 $info{total_clean_reads} M Clean Reads，各样品不少于 $info{min_clean_reads} M Clean Reads 。";
push @abstract,"(2) 检测到 $info{miRNA_num} 个miRNA，其中已知miRNA $info{miRNA_known_num} 个，新预测miRNA $info{miRNA_new_num} 个。";
push @abstract,"(3) 定量各样品miRNA表达丰度，获得差异表达miRNA筛选结果。";
if($deg_flag==1){
	push @abstract,"(4) 预测到miRNA靶基因 $info{target_gene_num} 个，并完成对差异表达miRNA靶基因的功能注释和富集分析。";
}else{
	push @abstract,"(4) 预测到miRNA靶基因 $info{target_gene_num} 个。";
}
my $abstract = join("</p><p class=\"p-abstract\">",@abstract);
$writer->emptyTag('report_abstract','value'=>"<p class= \"p-abstract\">$abstract</p>");

########################实验流程
$writer->emptyTag('h1','name'=>"1 实验流程",'type'=>'type1','desc'=>"1 实验流程");
$writer->emptyTag('p','desc'=>"RNA样品经过一系列严格的质量控制，合格的样品用于文库构建，文库构建严格按照NEB Next Ultra small RNA Sample Library Prep Kit for Illumina 试剂盒操作进行，合格的文库进行高通量测序。百迈客从样品的检测到文库制备的每一过程均经过严格的把关，从根本上获得可靠、高质量的分析数据。实验流程图如下所示：",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. 实验流程图",'type'=>"type1",'path'=>"HTML/template/library_miRNA_sep.png");
$writer->emptyTag('h2','name'=>"1.1 样品检测",'type'=>'type1','desc'=>"1.1 样品检测");
$writer->emptyTag('p','desc'=>"对RNA样品的检测主要包括以下几种方法：",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(1) Nanodrop检测：检测RNA样品的纯度（OD260/280≥1.8；OD260/230≥1.0）",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(2) Qubit 2.0检测：精确定量RNA样本的浓度（总 RNA 浓度≥250 ng/ul）",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(3) Agilent 2100 bioanalyzer检测：检测RNA样本的完整性等，以保证使用合格的样品进行测序（总RNA 的 RIN 值≥8.0， 28S/18S≥1.5；图谱基线无上抬； 5S 峰正常）",'type'=>"type1");
$writer->emptyTag('h2','name'=>"1.2 文库构建",'type'=>'type1','desc'=>"1.2 文库构建");
$writer->emptyTag('p','desc'=>"样品检测合格后，以2.5ng的量作为RNA样本起始量，用水补充体积至6ul，使用small RNA Sample Pre Kit试剂盒进行文库构建。由于Small RNA的5'端有磷酸基团，3'端有羟基，利用T4 RNA Ligase 1和T4 RNA Ligase 2（truncated）分别在small RNA 3'端和5'端连接上接头，反转录合成cDNA，PCR扩增，采用胶分离技术筛选目的片段，切胶回收得到的片段即为small RNA文库。文库构建流程如下：",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. 文库构建流程图",'type'=>"type1",'path'=>"HTML/template/library.png");
$writer->emptyTag('h2','name'=>"1.3 文库质控",'type'=>'type1','desc'=>"1.3 文库质控");
$writer->emptyTag('p','desc'=>"文库构建完成后，对文库质量进行检测，检测结果达到要求后方可进行上机测序，检测方法如下：",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(1) 使用Qubit2.0进行初步定量，使用Agilent 2100对文库的insert size进行检测，insert size符合预期后才可进行下一步实验。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(2) Q-PCR方法对文库的有效浓度进行准确定量（文库有效浓度＞2nM），完成库检。",'type'=>"type1");
$writer->emptyTag('h2','name'=>"1.4 上机测序",'type'=>'type1','desc'=>"1.4 上机测序");
$writer->emptyTag('p','desc'=>"库检合格后，不同文库按照目标下机数据量进行pooling，用Illumina HiSeq平台进行测序",'type'=>"type1");

####################2 生物信息学分析
$writer->emptyTag('h1','name'=>"2 生物信息学分析",'type'=>'type1','desc'=>"2 生物信息学分析");
$writer->emptyTag('p','desc'=>"small RNA项目分析流程示意图如下(如果为单样本，则没有差异分析相关内容)：",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. sRNA测序信息分析流程图",'type'=>"type1",'path'=>"HTML/template/d.workflow_miRNA_sep.png");
	##测序数据及其质量控制
$writer->emptyTag('h2','name'=>"2.1 测序数据及其质量控制",'type'=>'type1','desc'=>"2.1 测序数据及其质量控制");

$writer->emptyTag('p','desc'=>"sRNA文库所测的原始序列质量控制的标准为：",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(1) 对于每个样本，将质量值低的序列去掉。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(2) 去除未知碱基N（N为无法识别的碱基）含量大于等于10%的Reads；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(3) 去除没有3’接头序列的reads；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(4) 剪切掉3’接头序列；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(5) 去除短于$config{min_len}或长于$config{max_len}个核苷酸的序列；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"sRNA测序数据产出统计见下表：",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>"注：Sample_ID：客户样本名称；BMK-ID：百迈客编号；Raw_reads：测序原始数据；Low_quality_reads：质量值低于30的碱基所占比例超过20% 的reads；Containing'N'reads：含至少10% 的未知碱基N的reads；Length<$config{min_len}：去掉接头后，小于$config{min_len}个核苷酸的reads数；Length>$config{max_len}：去掉接头后，大于$config{max_len}个核苷酸的reads数；Clean_reads：质量值大于或等于30的碱基的Reads 数；Q30(%):质量值大于30的比例。",'type'=>"part",'name'=>"表$tid. sRNA测序数据统计表",'path'=>"BMK_1_rawData/BMK_1_Data_Assess/All_sample_filter.stat.xls");
$writer->emptyTag('p','desc'=>"通过质量控制，每个样品的Clean Data均大于 $info{min_clean_reads} M。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"测序数据质量评估详细结果见以下网页：",'type'=>"type1");
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"数据质量网页",'action'=>"xls",'path'=>"HTML/dataassess/dataassess.html",'type'=>"xls");
	$writer->emptyTag('file','desc'=>"",'name'=>"sRNA文库测序数据产出统计结果路径",'action'=>"xls",'path'=>"BMK_1_rawData/BMK_1_Data_Assess/",'type'=>"xls");
}else{
	$writer->emptyTag('file','desc'=>"",'name'=>"数据质量网页",'action'=>"xls",'path'=>"HTML/dataassess.xml",'type'=>"xml");
}
$writer->emptyTag('h2','name'=>"2.2 miRNA分析",'type'=>'type1','desc'=>"2.2 miRNA分析");
$writer->emptyTag('h3','name'=>"2.2.1 sRNA分类注释",'type'=>'type1','desc'=>"2.2.1 sRNA分类注释");
$writer->emptyTag('p','desc'=>"Bowtie<a href=\"#ref1\">[1]</a>软件是一种短序列比对软件，尤其适用于高通量测序获得的reads比对，百迈客利用Bowtie软件，将Clean Reads分别与Silva数据库、GtRNAdb数据库、Rfam数据库和Repbase数据库进行序列比对，过滤核糖体RNA（rRNA）、转运RNA（tRNA）、核内小RNA（snRNA）、核仁小RNA（snoRNA）等ncRNA以及重复序列，获得包含miRNA的Unannotated reads。sRNA注释分类统计见下表：",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>"",'type'=>"part",'name'=>"表$tid. sRNA分类注释统计",'path'=>"BMK_1_rawData/BMK_2_Mapped_Statistics/All_ncRNA_mapped.stat.xls");

$writer->emptyTag('h3','name'=>"2.2.2 与参考基因组比对",'type'=>'type1','desc'=>"2.2.2 与参考基因组比对");
my $ref_gene_addr = "<a href=\"$config{ref_gene_addr}\" title=\"click\" target=\"_blank\">$config{ref_gene_addr}</a>";
$writer->emptyTag('p','desc'=>"使用指定的 $config{ref_gene_name} 作为参考基因组进行序列比对及后续分析。下载地址: $ref_gene_addr 。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"利用Bowtie软件将Unannotated reads与参考基因组进行序列比对，获取在参考基因组上的位置信息，即为Mapped Reads。统计结果见下表：",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>'注：#BMK-ID：百迈客编号；Total_Reads：未被注释的用于与参考基因组比对的reads数目；Mapped_Reads：比对到参考基因组的Clean Reads；Mapped_reads(+)：比对正链上的Clean Reads数目； Mapped_reads(-)：比对到负链上的Clean Reads数目。','type'=>"part",'name'=>"表$tid. 与参考基因组比对信息统计表",'path'=>"BMK_1_rawData/BMK_2_Mapped_Statistics/All_sample_map.stat.xls");
$writer->emptyTag('p','desc'=>"将比对到不同染色体上的reads进行位置分布统计，绘制Mapped Reads在参考基因组上覆盖深度的分布图。样品的Mapped Reads在参考基因组染色体上的覆盖深度分布图如下：",'type'=>"type1");
&piclist("Reads在参考基因组上的位置及覆盖深度分布图","注：横坐标为染色体位置；纵坐标为染色体上对应位置的覆盖深度取以2为底的对数值。","$inpath/BMK_1_rawData/BMK_2_Mapped_Statistics/*/*.chro_distribution.png","BMK_1_rawData");

if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"sRNA分类注释和比对统计结果路径",'action'=>"xls",'path'=>"BMK_1_rawData/BMK_2_Mapped_Statistics/",'type'=>"xls");
}
$writer->emptyTag('h3','name'=>"2.2.3 miRNA预测",'type'=>'type1','desc'=>"2.2.3 miRNA预测");
$writer->emptyTag('p','desc'=>"在已知miRNA鉴定方面，我们将比对上参考基因组的reads序列与已知miRNA数据库miRBase（v21）中的成熟miRNA序列进行比对。序列与已知miRNA完全相同的reads被认为是本项目中鉴定到的已知miRNA。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"miRNA转录起始位点多位于基因间隔区、内含子以及编码序列的反向互补序列上，其前体具有标志性的发夹结构，成熟体的形成是由Dicer/DCL酶的剪切实现的。针对miRNA的生物特征，对于未鉴定到已知miRNA的序列，百迈客利用miRDeep2<a href=\"#ref2\">[2]</a>软件进行新miRNA的预测。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"我们利用miRDeep2软件包,通过reads比对到基因组上的位置信息得到可能的前体序列，基于reads在前体序列上的分布信息（基于miRNA产生特点，mature,star,loop）及前体结构能量信息（RNAfold randfold）采用贝叶斯模型经打分最终实现新miRNA的预测。miRDeep2主要用于动物miRNA的预测，通过参数的调整及打分系统的改变也可以对植物的miRNA进行预测。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"经分析：所有样品共得到 $info{miRNA_num} 个miRNA，其中已知miRNA $info{miRNA_known_num} 个，新预测miRNA $info{miRNA_new_num} 个。统计表格详见如下：");
$tid++;
$writer->emptyTag('table','desc'=>'注：BMK-ID：百迈客编号；Known-miRNA：已知miRNA数目；Novel-miRNAs：新预测的miRNA数量；Total：总的miRNA数量。','type'=>"part",'name'=>"表$tid. 各样品miRNA统计结果",'path'=>"$sRNA_relative_dir/BMK_1_miRNA_Prediction/BMK_4_Summary_Stat/Total_miRNA.stat.xls");
$writer->emptyTag('p','desc'=>"由于Dicer酶和DCL酶的特异性，最终生成的成熟miRNA长度主要集中在20nt到24nt的范围内，动物的miRNA以22nt为主。鉴定出的已知miRNA和新miRNA的长度分布见下图：",'type'=>"type1");
&piclist("miRNA的长度分布图","注：横坐标为miRNA长度，纵坐标为特定长度的miRNA数目。","$inpath/$sRNA_relative_dir/BMK_1_miRNA_Prediction/BMK_2_Len_Distribution/*_miRNA.length.png","$sRNA_relative_dir");

if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"miRNA长度分布统计结果路径",'action'=>"xls",'path'=>"$sRNA_relative_dir/BMK_1_miRNA_Prediction/BMK_2_Len_Distribution/",'type'=>"xls");
}
$writer->emptyTag('p','desc'=>'miRDeep2软件预测得到的新的miRNA的候选前体序列都有相应的pdf文件，如下图：','type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"注：图A为miRDeep2软件对miRNA前体的打分和比对到成熟序列、环状结构、star序列上的reads数及其前体的二级结构预测图：红色为成熟序列，黄色为环状结构，紫色为star序列；图B显示比对到本条前体上的reads分布；图C包含了成熟序列、环状结构、star序列的位置，紫色为miRDeep2软件预测的star序列，亮蓝色为测序reads支持的star序列。",'name'=>"图$pid. miRNA前体结构及序列深度结果文件",'type'=>"type1",'path'=>"HTML/template/miRDeep2.png");
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"miRNA二级结构预测结果路径",'action'=>"xls",'path'=>"$sRNA_relative_dir/BMK_1_miRNA_Prediction/BMK_5_PDF/",'type'=>"xls");
}
################
$writer->emptyTag('h3','name'=>"2.2.4 miRNA靶向关系预测",'type'=>'type1','desc'=>"2.2.4 miRNA靶向关系预测");
$writer->emptyTag('p','desc'=>"根据已知miRNA和新预测的miRNA与对应物种的基因序列信息，植物用 TargetFinder软件<a href=\"#ref3\">[3]</a>进行靶基因预测；动物用miRanda<a href=\"#ref4\">[4]</a>和targetscan<a href=\"#ref5\">[5]</a>进行靶基因预测。miRNA靶基因预测结果如下表：",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>'注：0表示该靶向关系对在该软件中并未预测到，1表示该靶向关系对在该软件中预测到，最终取两款软件的交集作为最终靶向关系对。','type'=>"part",'name'=>"表$tid. miRNA靶基因预测结果",'path'=>"HTML/template/miRNA_gene.xls");
$writer->emptyTag('file','desc'=>"注：第一列是miRNA,第二列是基因，如果靶向多个基因，以逗号分割。",'name'=>"miRNA靶基因关系表",'action'=>"xls",'path'=>"BMK_3_Target_Predict/mir2target.list.xls",'type'=>"xls");
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"miRNA靶向关系结果路径",'action'=>"xls",'path'=>"BMK_3_Target_Predict/",'type'=>"xls");
}
my $HMDD_dir = (split /$sRNA_relative_dir\//,(glob "$inpath/$sRNA_relative_dir/BMK_*_HMDD")[0])[1];
if(-e "$inpath/$sRNA_relative_dir/$HMDD_dir/HMDD_related_miRNA_with_target_gene.xls"){
        $writer->emptyTag('h3','name'=>"2.2.5 HMDD",'type'=>'type1','desc'=>"2.2.5 HMDD");
        $writer->emptyTag('p','desc'=>"HMDD（the Human microRNA Disease Database）<a href=\"#ref6\">[6]</a>：人疾病相关数据库，数据库提供了各种人类疾病中异常调控的microRNA，是一个综合资源数据库。数据库基于3511篇已发表论文共整理收纳了572个人类 miRNA 基因 和378 种疾病的10368对精选关系。自数据库2007年12月首个版本（HMDD v1.0）成功上线以来，目前已经累计更新30次以上。",'type'=>"type1");
        $writer->emptyTag('p','desc'=>"microRNA-疾病关系代表了异常调控microRNA在人类疾病中的病理作用。基于该数据库对本项目的样本进行注释，有利于更好的获的数据库中已收入疾病相关的信息，同时有利于发现与本项目相关miRNA的调控关系。部分结果展示见下表：",'type'=>"type1");
        $tid++;
        $writer->emptyTag('table','desc'=>'','type'=>"part",'name'=>"表$tid. HMDD注释",'path'=>"HTML/template/hmdd.xls");
        $writer->emptyTag('file','desc'=>"",'name'=>"miRNA人疾病相关数据库HMDD注释",'action'=>"xls",'path'=>"$sRNA_relative_dir/$HMDD_dir/HMDD_related_miRNA_with_target_gene.xls",'type'=>"xls");
	if(!defined $cloud){
		$writer->emptyTag('file','desc'=>"",'name'=>"HMDD分析结果路径",'action'=>"xls",'path'=>"BMK_2_miRNA/$HMDD_dir/",'type'=>"xls");
	}
}

	###miRNA 差异分析
$writer->emptyTag('h2','name'=>"2.3 miRNA表达量和差异表达及功能分析",'type'=>'type1','desc'=>"2.3 miRNA表达量和差异表达及功能分析");
$writer->emptyTag('p','desc'=>"关于miRNA表达量的整体分布，差异表达miRNA筛选及差异表达miRNA的靶基因功能注释和富集分析（单样品无差异相关分析）可详见以下网页：",'type'=>"type1");
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"miRNA差异表达及功能分析网页（单样品无差异相关分析）",'action'=>"xls",'path'=>"HTML/miRNA/miRNA.html",'type'=>"xls");
}else{
	$writer->emptyTag('file','desc'=>"",'name'=>"miRNA差异表达及功能分析网页（单样品无差异相关分析）",'action'=>"xls",'path'=>"HTML/miRNA.xml",'type'=>"xml");
}

$pid++;
$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. 差异miRNA及靶基因功能富集分析示意图",'type'=>"type1",'path'=>"HTML/template/2.noncodeDEG.png");

$writer->emptyTag('h2','name'=>"2.4 结构分析",'type'=>'type1','desc'=>"2.4 结构分析");
	###miRNA碱基偏好性
$writer->emptyTag('h3','name'=>"2.4.1 miRNA碱基偏好性",'type'=>'type1','desc'=>"2.4.1 miRNA碱基偏好性");
$writer->emptyTag('p','desc'=>"Dicer酶和DCL酶在识别和切割前体miRNA时，5’端首位碱基对U具有很强的偏向性。通过对miRNA的碱基偏好性分析，获得典型的miRNA碱基比例。新预测的miRNA的5'端首位碱基偏好和各位点的碱基偏好统计见下图：");
&piclist("不同长度miRNA的首位碱基分布图","注：横坐标表示不同长度的序列；纵坐标表示不同长度miRNA首位碱基所占百分比","$inpath/BMK_2_miRNA/BMK_1_miRNA_Prediction/BMK_3_Base_Distribution/*.first_base.png","BMK_2_miRNA");
&piclist("miRNA各位点碱基分布图","注：横坐标表示序列的各位点；纵坐标表示miRNA各位点碱基所占百分比","$inpath/BMK_2_miRNA/BMK_1_miRNA_Prediction/BMK_3_Base_Distribution/*.base_bias.png","BMK_2_miRNA");
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"miRNA碱基偏好性分析结果路径",'action'=>"xls",'path'=>"$sRNA_relative_dir/BMK_1_miRNA_Prediction/BMK_3_Base_Distribution/",'type'=>"xls");
}
	###miRNA碱基编辑
my $Edit_dir = (split /$sRNA_relative_dir\//,(glob "$inpath/$sRNA_relative_dir/BMK_*_miRNA_Edit")[0])[1];
$h3=2;
$writer->emptyTag('h3','name'=>"2.4.$h3 miRNA碱基编辑",'type'=>'type1','desc'=>"2.4.$h3 miRNA碱基编辑");
$writer->emptyTag('p','desc'=>"miRNA存在转录后碱基的编辑，导致种子序列（seed sequence）改变，进而使得作用的靶基因发生改变。百迈客使用isomiRID<a href=\"#ref7\">[7]</a>软件检测发生碱基编辑的miRNA，isomiRID软件首先调用bowtie和前体序列进行第一轮（r0）比对，完全比对的序列作为下轮比对的参考模板。第二轮（r1）比对时，允许一个错配，3'端有一个错配的记为M3, 5'端有一个错配的记为M5,中间有一个错配的记为MM。部分miRNA编辑分析结果如下：",'type'=>"type1");
$writer->emptyTag('file','desc'=>"注：第一列是小RNA及pre-miRNA的序列及匹配情况，第二列表示round数及碱基编辑类型，第三列表示小RNA的长度，第四列表示相对于前体的碱基变异情况，第五列是小RNA在前体的位置，六到七表示小RNA在各样本中的count数。",'name'=>"miRNA碱基编辑结果",'action'=>"xls",'path'=>"$sRNA_relative_dir/$Edit_dir/miRNAEdit_cutoff_20.xls",'type'=>"xls");

	###miRNA家族
my $Family_dir = (split /$sRNA_relative_dir\//,(glob "$inpath/$sRNA_relative_dir/BMK_*_miRNA_Family")[0])[1];
$h3++;
$writer->emptyTag('h3','name'=>"2.4.$h3 miRNA家族分析",'type'=>'type1','desc'=>"2.4.$h3 miRNA家族分析");
$writer->emptyTag('p','desc'=>"miRNA在物种间具有高度保守性，基于序列的相似性对检测到的已知miRNA和新miRNA进行miRNA家族分析，研究miRNA在进化中的保守性。miRNA家族分析部分结果见下表：",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>'','type'=>"part",'name'=>"表$tid. miRNA家族分析",'path'=>"HTML/template/miRNA_family.xls");
$writer->emptyTag('file','desc'=>"注：第一列为miRNA ID，第二列为miRNA家族信息。",'name'=>"miRNA家族分析结果",'action'=>"xls",'path'=>"$sRNA_relative_dir/$Family_dir/Family_In_miR.xls",'type'=>"xls");

############################
$writer->emptyTag('h1','desc'=>'3. 结果文件查看说明','type'=>"type1",'name'=>'3. 结果文件查看说明');
$writer->emptyTag('p','desc'=>"(1) 上传目录中有Readme.txt说明，详细介绍了每个文件所代表的内容。上传的结果数据文件多以文本格式为主(fa文件、txt文件，detail文件，xls文件等)。在Windows系统下查看文件，推荐使用Editplus或UltraEdit或Notepad作为文本浏览程序，否则会因文件过大造成死机。如果文件大小超过50M建议用pilotedit打开；在Unix或Linux系统下可以浏览较大的文本文件，用less等操作命令可以顺利地查看。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(2) 阅读报告时请注意：报告中出现的路径链接，只有所有结果文件都在的时候才可以查看，只有单独的结题报告而没有结果文件时是无法查看的。",'type'=>"type1");

$writer->startTag('ref_list',desc=>"参考文献",type=>"type1",name=>"参考文献");
&reference ($writer,[
	["Langmead, B., Trapnell, C., Pop, M., and Salzberg, S.L. (2009). Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome biology 10,R25.","http://bowtie-bio.sourceforge.net/index.shtml","1"],
	["Friedlander,M.R., Mackowiak,S.D., Li,N., Chen,W. and Rajewsky,N. (2012) miRDeep2 accurately identifies known and hundreds of novel microRNA genes in seven animal clades. Nucleic Acids Res., 40, 37-52.","https://academic.oup.com/nar/article/40/1/37/1275937","2"],
	["Bo X, Wang S. TargetFinder: a software for antisense oligonucleotide target site selection based on MAST and secondary structures of target mRNA[J]. Bioinformatics, 2005, 21(8):1401.","https://academic.oup.com/bioinformatics/article/21/8/1401/250353","3"],
	["Betel D, Wilson M, Gabow A, et al. The microRNA. org resource: targets and expression[J]. Nucleic acids research, 2008, 36(suppl 1): D149-D153.","https://academic.oup.com/nar/article/36/suppl_1/D149/2508385","4"],
	["Lewis B P, Shih I H, Jonesrhoades M W, et al. Prediction of mammalian microRNA targets.[J]. Cell, 2003, 115(7):787-798.","https://www.cell.com/cell/pdf/S0092-8674(03)01018-3.pdf","5"],
	["Dr. Qinghua Cui, 38 Xueyuan Rd, Department of Biomedical Informatics, Peking University Health Science Center, Beijing 100191, China","http://210.73.221.6/hmdd#fragment-4","6"],
	["de Oliveira L F, Christoff A P, Margis R. isomiRID: a framework to identify microRNA isoforms.[J]. Bioinformatics, 2013, 29(20):2521-3.","https://academic.oup.com/bioinformatics/article/29/20/2521/276800","7"],
]);
$writer->endTag('ref_list');
#<a href=\"#ref1\">[1]</a><a href=\"#ref2\">[2]</a>

$writer->emptyTag('h1','desc'=>'附录','type'=>'type1','name'=>'附录');
$writer->emptyTag('h2','desc'=>'1 英文材料、方法与数据分析','type'=>'type1','name'=>'1 英文材料、方法与数据分析');
$writer->emptyTag('p','desc'=>"英文描述见以下附件:",'type'=>"type1");
$writer->emptyTag('file','desc'=>"",'name'=>"英文材料、方法与数据分析",'action'=>"xls",'path'=>"HTML/template/materials_and_methods.pdf",'type'=>"xls");

$writer->endTag('report');
open OUT,">:utf8", "$output";
my $xmlstr=&decorate($writer->to_string);
$xmlstr = Encode::decode("utf-8", $xmlstr);
print OUT $xmlstr;
close(OUT);
$writer->end();


################################################################################################
#
#		                        sub function
################################################################################################
sub decorate{
	my $xmlstr=shift;
	$xmlstr=~s/(<[^>]*)\>/$1\>\n/mgo;
	return $xmlstr;
}

sub reference{
    my ($writer,$list)=@_;
    my ($data,$value,$unit,$item);
    
    foreach $data (@$list)    {		
	$writer->emptyTag('ref','name'=>$$data[0],'link'=>$$data[1],'id'=>$$data[2]);
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
	$writer->startTag('file_list','name'=>"",'desc'=>"$desc",'type'=>"type1");
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
sub relate{    
        my $file=shift;
        open(REL,$file)||die $!;
        while(<REL>){
                chomp;
                next if($_ !~/^SampleID/);
                my @tmp=split(/\s+/,$_);
                $sample{$tmp[1]}=$tmp[2];
        }
	close(REL);
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



sub USAGE{
        my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
        Options:
	-i	<path>	input images path, forced
        -o      <path>	output xml name, default inpath/configtest.xml
	-cfg1	<file>	data_cfg
	-cfg2	<file>	detai_cfg

        -h      Help

Example: perl $0 -i BMK_Result -o BMK_Result/configtest.xml -cfg1 data.cfg -cfg detail.cfg

USAGE
        print $usage;
        exit;
}

