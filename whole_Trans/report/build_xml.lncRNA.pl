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
my ($lnc_relative_dir,$circ_relative_dir,$mRNA_relative_dir);
my $lnc_dir = (glob "$inpath/BMK_*_LncRNA")[0];
my $mRNA_dir = (glob "$inpath/BMK_*_mRNA")[0];
my $circ_dir = (glob "$inpath/BMK_*_circRNA")[0];
my $struc_dir = (glob "$inpath/BMK_*_Structure")[0];

if(defined $lnc_dir){
	$lnc_relative_dir = (split /\//,$lnc_dir)[-1];
	$mRNA_relative_dir = (split /\//,$mRNA_dir)[-1];
}
if(defined $circ_dir){
	$circ_relative_dir =(split /\//,$circ_dir)[-1];
}
if(defined $struc_dir){
	$struc_relative_dir =(split /\//,$struc_dir)[-1];
}
my %config=&readConfig($detail_cfg);
my %info=&readConfig("$inpath/HTML/template/stat.info");

my $deg_flag=1;
if(($sam_num ==1 )|| (!exists $config{Sep} && !exists $config{Com})){
	$deg_flag=0;
}

###摘要描述
my ($disc1,$disc2,$disc3,$disc4,$disc5,$disc6);
if(defined $lnc_dir && defined $circ_dir){
        $disc1 = "去核糖体链特异文库的构建及数据测序，该数据可以同时用于分析长链非编码RNA，circRNA以及mRNA";
        $disc2 = "进行lncRNA，circRNA及新基因预测。共检测到 $info{lncRNA_num} 个lncRNA，其中有 $info{lncRNA_known_num} 个已知lncRNA，$info{lncRNA_new_num} 个新的lncRNA；共检测到 $info{circRNA_num} 个circRNA，其中有 $info{circRNA_known_num} 个已知circRNA，$info{circRNA_new_num} 个新的circRNA；共检测到 $info{gene_num} 个genes，其中有 $info{gene_known_num} 个已知genes，其中有 $info{gene_new} 个新的genes";
	$disc3 = "进行SNP、可变剪切、基因结构优化、融合基因、保守性分析";
	$disc4 = "对差异表达的基因、circRNA的来源基因、lncRNA的靶基因进行功能注释和富集分析";
	$disc5 = "比对效率从$info{lnc_map_min}到$info{lnc_map_max}不等";
	$disc6 = "对每种RNA进行表达定量";
}elsif(defined $lnc_dir && !defined $circ_dir){
        $disc1 = "去核糖体链特异文库的构建及数据测序分析";
	$disc2 = "进行lncRNA及新基因预测。共检测到 $info{lncRNA_num} 个lncRNA，其中有 $info{lncRNA_known_num} 个已知lncRNA，$info{lncRNA_new_num} 个新的lncRNA； $info{gene_num} 个genes，其中有 $info{gene_known_num} 个已知genes，$info{gene_new} 个新的genes";
	$disc3 = "进行SNP、可变剪切、融合基因、基因结构优化分析";
	$disc4 = "对差异表达的基因、lncRNA的靶基因进行功能注释和富集分析";
	$disc5 = "比对效率从$info{lnc_map_min}到$info{lnc_map_max}不等";
	$disc6 = "对每种RNA进行表达定量";
}elsif(!defined $lnc_dir && defined $circ_dir){
        $disc1 = "circRNA文库的构建及数据测序分析";
	$disc2 = "进行circRNA预测，共识别到 $info{circRNA_num} 个circRNA，其中有 $info{circRNA_known_num} 个已知circRNA，$info{circRNA_new_num} 个新的circRNA";
	$disc3 = "完成circRNA的可变剪切分析，环状RNA的保守性分析等";
	$disc4 = "对差异表达circRNA的来源基因进行功能注释和富集分析";
	$disc5 = "比对效率从$info{circ_map_min}到$info{circ_map_max}不等";
	$disc6 = "对circRNA进行表达定量";
}

my ($pid,$tid)=(0,0);

my $report_time=&gaintime();
my $writer = XML::Writer->new(OUTPUT => 'self');
$writer->xmlDecl('UTF-8');
$writer->startTag('report');
$writer->emptyTag('report_version','value'=>$config{report_version});
my $report_name;
$report_name=$config{Project_name}	if($config{Project_name});
$report_name ||= "RNA分析报告";
$writer->emptyTag('report_name','value'=>$report_name);
$writer->emptyTag('report_code','value'=>$config{Project_id});
$writer->emptyTag('report_user','value'=>$config{report_user});
$writer->emptyTag('report_user_addr','value'=>$config{report_user_addr});
$writer->emptyTag('report_time','value'=>"XXXX/XX/XX;XXXX/XX/XX;XXXX/XX/XX;$report_time");
$writer->emptyTag('report_abstract','value'=>"");
###############################摘要
my @abstract;

push @abstract,"(1) 完成 $sam_num 个样品，共 $info{lncRNA_total} Gb Clean Data的$disc1。各样品Clean Data均达到 $info{lncRNA_min} Gb及以上，Q30碱基百分比在 $info{lncRNA_Q30}% 及以上。";
push @abstract,"(2) 各样品分别与参考基因组序列比对，$disc5。";
push @abstract,"(3) 基于比对结果，$disc2。";
push @abstract,"(4) $disc3。";
if($deg_flag==1){
	push @abstract,"(5) $disc6，并进行差异分析，$disc4 。";
}else{
	push @abstract,"(5) $disc6。";
}

my $abstract = join("</p><p class=\"p-abstract\">",@abstract);
$writer->emptyTag('report_abstract','value'=>"<p class= \"p-abstract\">$abstract</p>");

########################实验流程
$writer->emptyTag('h1','name'=>"1 实验流程",'type'=>'type1','desc'=>"1 实验流程");
$writer->emptyTag('h2','name'=>"1.1 RNA样品检测",'type'=>'type1','desc'=>"1.1 RNA样品检测");
$writer->emptyTag('p','desc'=>"RNA样品经过一系列严格的质量控制，合格的样品用于文库构建。样品的检测到文库制备的每一过程我们均经过严格的把关，从根本上获得可靠、高质量的分析数据。实验流程见下图：",'type'=>"type1");
$pid++;
###实验流程图
my $library_png_lncRNA = "HTML/template/library_lncRNA_sep.png";
my $library_png_RNase = "HTML/template/library_RNase.jpg";
if(defined $lnc_dir){
	$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. 测序实验流程图",'type'=>"type1",'path'=>"$library_png_lncRNA");
}else{
	if(defined $config{Lib_type2} && $config{Lib_type2} eq "no"){
		$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. 测序实验流程图",'type'=>"type1",'path'=>"$library_png_lncRNA");
	}else{
		$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. 测序实验流程图",'type'=>"type1",'path'=>"$library_png_RNase");
	}
}

$writer->emptyTag('h2','name'=>"1.2 RNA文库构建",'type'=>'type1','desc'=>"1.2 RNA文库构建");
$writer->emptyTag('p','desc'=>"样品检测合格后，进行文库构建，主要流程如下：",'type'=>"type1");
###建库方法描述
my @experiment_lncRNA = "(1) 利用epicentre Ribo-ZeroTM试剂盒去除样品rRNA；";
push @experiment_lncRNA,"(2) 加入Fragmentation Buffer将rRNA-depleted RNA随机打断；";
push @experiment_lncRNA,"(3) 以rRNA-depleted RNA为模板，用六碱基随机引物（random hexamers）合成第一条cDNA链，然后加入缓冲液、dATP、dUTP、dCTP、dGTP 、RNase H和DNA polymerase I合成第二条cDNA链，利用AMPure XP beads纯化cDNA；";
push @experiment_lncRNA,"(4) 纯化的双链cDNA再进行末端修复、加A并连接测序接头，然后用AMPure XP beads进行片段大小选择；";
push @experiment_lncRNA,"(5) 最后降解含U链，最后通过PCR富集得到cDNA文库。";
my $experiment_lncRNA = join("</p><p class=\"paragraph\">",@experiment_lncRNA);
$experiment_lncRNA = join("</p><p class=\"final-p\">",@experiment_lncRNA) if (defined $cloud);

my @experiment_RNase = "(1) 用rRNA探针去除rRNA；";
push @experiment_RNase,"(2) Rnase R消化去除线性RNA；";
push @experiment_RNase,"(3) 加入片段化试剂(Fragmentation Buffer)将circRNA进行随机打断；";
push @experiment_RNase,"(4) 以片段化的circRNA为模板，用六碱基随机引物(random hexamers)合成第一条链，随后加入缓冲液、dNTPs、RNase H和DNA polymerase I合成cDNA的第二链，再利用AMPure XP beads纯化circRNA；";
push @experiment_RNase,"(5) 利用T4 DNA聚合酶和Klenow DNA聚合酶活性将DNA的粘性末端修复为平末端，3'末端加A尾并连接测序接头，随后用AMPure XP beads进行片段大小选择，之后用USER酶降解含有U的cDNA第二链；";
push @experiment_RNase,"(6) 最后通过PCR扩增得到circRNA文库 。";
my $experiment_RNase = join("</p><p class=\"paragraph\">",@experiment_RNase);
$experiment_RNase = join("</p><p class=\"final-p\">",@experiment_RNase) if(defined $cloud);

if(defined $lnc_dir){
	$writer->emptyTag('p','desc'=>"$experiment_lncRNA",'type'=>"type1");
}else{
	if(defined $config{Lib_type2} && $config{Lib_type2} eq "no"){
		$writer->emptyTag('p','desc'=>"$experiment_lncRNA",'type'=>"type1");
	}else{
		$writer->emptyTag('p','desc'=>"$experiment_RNase",'type'=>"type1");
	}
}

$writer->emptyTag('h2','name'=>"1.3 文库质控",'type'=>'type1','desc'=>"1.3 文库质控");
$writer->emptyTag('p','desc'=>"文库构建完成后，对文库质量进行检测，检测结果达到要求后方可进行上机测序，检测方法如下：",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(1) 使用Qubit2.0进行初步定量，使用Agilent 2100对文库的insert size进行检测，insert size符合预期后才可进行下一步实验。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(2) Q-PCR方法对文库的有效浓度进行准确定量（文库有效浓度＞2nM），完成库检。",'type'=>"type1");
$writer->emptyTag('h2','name'=>"1.4 上机测序",'type'=>'type1','desc'=>"1.4 上机测序");
$writer->emptyTag('p','desc'=>"库检合格后，不同文库按照目标下机数据量进行pooling，用Illumina HiSeq平台进行测序。",'type'=>"type1");
####################2 生物信息学分析
$writer->emptyTag('h1','name'=>"2 生物信息学分析",'type'=>'type1','desc'=>"2 生物信息学分析");
$writer->emptyTag('p','desc'=>"分析流程示意图如下（如果为单样本，则没有差异分析相关内容）：",'type'=>"type1");
$pid++;
$writer->startTag('pic_list','name'=>"图$pid. 生物信息分析流程图",'desc'=>"",'type'=>"type1");
if(defined $lnc_dir && defined $circ_dir){
	$writer->emptyTag('pic','desc'=>"",'name'=>"生物信息分析流程图（mRNA部分）",'type'=>"type1",'path'=>"HTML/template/a.workflow_gene.png");
	$writer->emptyTag('pic','desc'=>"",'name'=>"生物信息分析流程图（lncRNA部分）",'type'=>"type1",'path'=>"HTML/template/b.workflow_lncRNA.png");
	$writer->emptyTag('pic','desc'=>"",'name'=>"生物信息分析流程图（circRNA部分）",'type'=>"type1",'path'=>"HTML/template/c.workflow_circRNA.png");
}elsif(defined $lnc_dir && !defined $circ_dir){
        $writer->emptyTag('pic','desc'=>"",'name'=>"生物信息分析流程图（mRNA部分）",'type'=>"type1",'path'=>"HTML/template/a.workflow_gene.png");
        $writer->emptyTag('pic','desc'=>"",'name'=>"生物信息分析流程图（lncRNA部分）",'type'=>"type1",'path'=>"HTML/template/b.workflow_lncRNA.png");
}elsif(!defined $lnc_dir && defined $circ_dir){
        $writer->emptyTag('pic','desc'=>"",'name'=>"circRNA生物信息分析流程图",'type'=>"type1",'path'=>"HTML/template/c.workflow_circRNA_sep.png");
}
$writer->endTag('pic_list');
	##测序数据及其质量控制
$writer->emptyTag('h2','name'=>"2.1 测序数据及其质量控制",'type'=>'type1','desc'=>"2.1 测序数据及其质量控制");
$writer->emptyTag('p','desc'=>"测序得到的原始序列含有接头序列或低质量序列，为了保证信息分析的准确性，需要对原始数据进行质量控制，得到高质量序列（即Clean Reads），去核糖体文库所测的原始序列质量控制的标准为：",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(1) 去除含接头的reads；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(2) 去除低质量的Reads（包括去除N的比例大于10%的Reads；去除质量值Q≤10的碱基数占整条Read的50%以上的Reads）。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"经过上述一系列的质量控制之后得到的高质量reads，称为Clean Data，Clean Data同样以FASTQ格式提供。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"该项目各样品建库数据产出统计见下表：",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>'注：Sample_ID：客户样品编号；BMK-ID：百迈客对样品的统一编号；ReadSum：Clean Data中pair-end Reads总数；BaseSum：Clean Data总碱基数；GC(%)：Clean Data GC含量，即Clean Data中G和C两种碱基占总碱基的百分比；N(%)：Clean Data中未分辨碱基占总碱基的百分比；Q30(%)：Clean Data质量值大于或等于Q30的碱基所占的百分比。','type'=>"part",'name'=>"表$tid. 测序数据统计表",'path'=>"BMK_1_rawData/BMK_1_Data_Assess/nonRib.AllSample_GC_Q.stat.xls");
$writer->emptyTag('p','desc'=>"经过测序质量控制，共得到$info{lncRNA_total}Gb Clean Data，各样品Q30碱基百分比均不小于$info{lncRNA_Q30}%。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"测序数据质量评估详细结果见以下网页：",'type'=>"type1");
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"数据质量网页",'action'=>"xls",'path'=>"HTML/dataassess/dataassess.html",'type'=>"xls");
	$writer->emptyTag('file','desc'=>"",'name'=>"文库测序数据质控结果路径",'action'=>"xls",'path'=>"BMK_1_rawData/BMK_1_Data_Assess/",'type'=>"xls");
}else{
	$writer->emptyTag('file','desc'=>"",'name'=>"数据质量网页",'action'=>"xls",'path'=>"HTML/dataassess.xml",'type'=>"xml");
}
#######################比对，预测
my ($Path1,$Path2);
if(defined $lnc_dir && defined $circ_dir){
	$Path1 = "BMK_1_rawData/BMK_2_Mapped_Statistics/lncRNA";
	$Path2 = "BMK_1_rawData/BMK_3_Library_Assessment/lncRNA";
}elsif(defined $lnc_dir && !defined $circ_dir){
	$Path1 = "BMK_1_rawData/BMK_2_Mapped_Statistics";
	$Path2 = "BMK_1_rawData/BMK_3_Library_Assessment";
}elsif(!defined $lnc_dir && defined $circ_dir){
	$Path1 = "BMK_1_rawData/BMK_2_Mapped_Statistics";
	$Path2 = "BMK_1_rawData/BMK_3_Library_Assessment";
}

my $h=1;
if(defined $lnc_dir){
	$h++;
	###gene
$writer->emptyTag('h2','name'=>"2.$h 基因分析",'type'=>'type1','desc'=>"2.$h 基因分析");
my $ref_gene_addr = "<a href=\"$config{ref_gene_addr}\" title=\"click\" target=\"_blank\">$config{ref_gene_addr}</a>";
$writer->emptyTag('p','desc'=>"使用指定的基因组$config{ref_gene_name}作为参考进行序列比对及后续分析。$config{ref_gene_name}下载地址：$ref_gene_addr 。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"利用HISAT2<a href=\"#ref2\">[2]</a>软件进行比对，HISAT是一个来自RNA测序实验reads的高效比对系统（TopHat2/Bowtie2的继任者）。HISAT使用一个基于Burrows-Wheeler变换和Ferragina-Manzini（FM）索引的索引方案，利用了两类索引进行比对：一个全基因组FM索引以定位每个比对，和许多局部FM索引用于非常快地扩展这些比对，实现了更快的速度和更少的资源占用的分析。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"比对完成后采用StringTie<a href=\"#ref3\">[3]</a>软件对比对上的read进行组装，StringTie是基于最优化理论建立的算法，利用比对信息构建多可变剪切图谱运用构建流量网络从而根据最大流量算法来对reads进行组装和评估其表达量，同Cufflinks等其他软件相比可以构建更完整的转录本和更好的评估表达量。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"分析流程图<a href=\"#ref3\">[3]</a>如下:");
$pid++;
$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. HISAT和StringTie分析流程",'type'=>"type1",'path'=>"HTML/template/hisat2.png");
	###Map
$writer->emptyTag('h3','name'=>"2.$h.1 比对结果统计及作图",'type'=>'type1','desc'=>"2.$h.1 比对结果统计及作图");
$writer->emptyTag('p','desc'=>"比对效率指Mapped Reads占Clean Reads的百分比，它除了受数据测序质量影响外，还与指定的参考基因组组装的优劣、参考基因组与测序样品的生物学分类关系远近（亚种）有关。通过比对效率，可以评估所选参考基因组是否能满足信息分析的需求。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"各样品测序数据与所选参考基因组的序列比对结果统计见下表：",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>"注：BMK-ID：百迈客对样品的统一编号；Total Reads：Clean Reads数目，按单端计；Mapped Reads：比对到参考基因组上的Reads数目及在Clean Reads中占的百分比；Uniq Mapped Reads：比对到参考基因组唯一位置的Reads数目及在Clean Reads中占的百分比；Multiple Mapped Reads：比对到参考基因组多个位置的数目及在Clean Reads中占的百分比；Reads Map to '+'：比对到参考基因组正链的Reads数目及在Clean Reads中占的百分比。Reads Map to '-'：比对到参考基因组负链的Reads数目及在Clean Reads中占的百分比。",'type'=>"part",'name'=>"表$tid. 样品测序数据与所选参考基因组的序列比对结果统计表",'path'=>"$Path1/All.mappedStat.xls");
$writer->emptyTag('p','desc'=>"长链非编码RNA测序为双端链特异性建库测序，链特异性测序可以确定两条链的转录方向性。利用RSeQC对测序数据和基因在基因组上的位置的相对应的关系进行分析，得到以下统计结果。见下表：",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>'注：第一列sample_id：样品名称；第二列read1_vs_gene：该列数值代表read1比对到参考基因组的位置链信息(+/-)同相应位置基因的链信息一致的read1数量占总的read1的比例；第三列read２_vs_gene：该列数值代表read2比对到参考基因组的位置链信息(+/-)同相应位置基因的链信息一致的read２数量占总的read２的比例；第四列undefined：在基因组中可能含有双向的转录本存在，该区域的read无法判断是否同源基因链信息是否一致，所以命名为undefined．','type'=>"part",'name'=>"表$tid. 数据链特异性建库评估信息表",'path'=>"$Path2/Lib_assess.xls");

#if(!defined $cloud){
#	$writer->emptyTag('file','desc'=>"",'name'=>"LncRNA文库测序数据链特异性评估结果路径",'action'=>"xls",'path'=>"$Path2",'type'=>"xls");
#}
$writer->emptyTag('p','desc'=>"将比对到不同染色体上Reads进行位置分布统计，绘制Mapped Reads在所选参考基因组上的覆盖深度分布图。各样品的Mapped Reads在参考基因组上的覆盖深度分布图如下：",'type'=>"type1");
&piclist("Mapped Reads在参考基因组上的位置及覆盖深度分布图","注：横坐标为染色体位置；纵坐标为覆盖深度以2为底的对数值，以10kb作为区间单位长度，划分染色体成多个小窗口（Window），统计落在各个窗口内的Mapped Reads作为其覆盖深度。蓝色为正链，绿色为负链。","$inpath/$Path1/*.map.png","BMK_1_rawData");
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"lncRNA比对效率统计和覆盖深度分布结果路径",'action'=>"xls",'path'=>"$Path1",'type'=>"xls");
}
$writer->emptyTag('h3','name'=>"2.$h.2 长链非编码文库质量评估",'type'=>'type1','desc'=>"2.$h.2 长链非编码文库质量评估");
$writer->emptyTag('p','desc'=>"mRNA片段化后的插入片段大小选择，是从mRNA序列中独立随机地抽取子序列，mRNA数目越大、打断方式和时间控制得越合适，目的RNA每个部分被抽取到的可能性就越接近，mRNA片段化随机性越高，mRNA上覆盖的Reads越均匀。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"通过Mapped Reads在各mRNA转录本上的位置分布，模拟mRNA片段化结果，检验mRNA片段化的随机程度。如果mRNA存在严重降解，被降解的碱基序列不能被测序，即无Reads比对上。因此，通过查看Mapped Reads在mRNA转录本上的位置分布可了解mRNA的降解情况。样品Mapped Reads在mRNA转录本上的位置分布如下图：",'type'=>"type1");
&piclist("Mapped Reads在mRNA上的位置分布图","注：横坐标为标准化后的mRNA位置，纵坐标为对应位置区间内Reads在总Mapped Reads中所占百分比。由于参考的mRNA长度不同，作图时把每个mRNA按照长度划分成100个区间，进而统计每一区间内的Mapped Reads数目及所占的比例，图中反映的是所有mRNA各个区间内的Mapped Reads比例的汇总。","$inpath/$Path2/*.randcheck.png","BMK_1_rawData");
$writer->emptyTag('p','desc'=>"插入片段长度的离散程度能直接反映出文库制备过程中磁珠纯化的效果。通过插入片段两端的Reads在参考基因组上的比对起止点之间的距离计算插入片段长度。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"大部分的真核生物基因为断裂基因，外显子被内含子隔断，而长链非编码测序得到的是无内含子的成熟mRNA。当mRNA中跨内含子的片段两端的Reads比对到基因组上时，比对起止点之间的距离要大于插入片段长度。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"各样品的插入片段长度模拟分布图如下：",'type'=>"type1");
&piclist("插入片段长度模拟分布图","注：横坐标为双端Reads在参考基因组上的比对起止点之间的距离，范围为0到800bp；纵坐标为比对起止点之间不同距离的双端Reads或插入片段数量。","$inpath/$Path2/*.insertSize.png","BMK_1_rawData");
$writer->emptyTag('p','desc'=>"上图中曲线的主峰落在200bp附近，没有偏离目标区域，说明插入片段长度的离散程度较小，插入片段大小选择正常。",'type'=>"type1");

$writer->emptyTag('p','desc'=>"为了评估数据是否充足并满足后续分析，对测序得到的基因数进行饱和度检测。由于一个物种的基因数目是有限的，且基因转录具有时间和空间特异性，因此随着测序量的增加，检测到的基因数目会趋于饱和。对于表达量越高的基因，越容易被检测定量。因此，对于表达量越低的基因，需要更大的数据量才能被准确定量。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"使用各样品的Mapped Data对检测到的不同表达情况的基因数目饱和情况进行模拟，绘制曲线图如下，可查看随着测序数据量的增加，检测到的不同表达量的基因数目是否趋于饱和。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"使用各样品的Mapped Data对检测到的基因数目的饱和情况进行模拟，绘制曲线图如下：",'type'=>"type1");
&piclist("长链非编码测序数据饱和度模拟图","注：通过将Mapped Reads等量地分成10份，逐渐增加数据查看检测到的基因表达量来绘制饱和度曲线。横坐标为随机抽取的Mapped Reads数目占总的Mapped Reads数的百分比，纵坐标为由随机抽取一定比例数据中检测到在某一表达量区间的基因同该基因的表达量终值(所有数据评估的表达量值)误差在15%之内的基因数比例。纵坐标越趋近去1.0表示表达量越趋于饱和,图中不同颜色的曲线代表不同表达量水平基因额饱和度情况;表达量FPKM不小于0.1的基因定义为表达的基因。","$inpath/$Path2/*.Saturation.png","BMK_1_rawData");
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"LncRNA文库质量及链特异性评估结果路径",'action'=>"xls",'path'=>"$Path2",'type'=>"xls");
}
	###新基因预测
$writer->emptyTag('h3','name'=>"2.$h.3 新基因预测",'type'=>'type1','desc'=>"2.$h.3 新基因预测");
$writer->emptyTag('p','desc'=>"基于所选参考基因组序列，使用StringTie软件对Mapped Reads进行拼接，并与原有的基因组注释信息进行比较，寻找原来未被注释的转录区，发掘该物种的新转录本和新基因，从而补充和完善原有的基因组注释信息。过滤掉编码的肽链过短（少于50个氨基酸残基）或只包含单个外显子的序列，共发现$info{new}个新基因。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"新基因的GFF格式文件部分见下表：",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>'注#Seq_ID：染色体号；Source：注释信息的来源，StringTie软件；Type：注释特征（Feature）类型；Start/End：特征序列的起止位置；Score：得分，数字，注释信息可能性的说明，“.”表示缺失值；Strand：特征序列所在的正负链；Phase：仅对注释类型为CDS有效，表示起始编码的位置，有效值为0、1、2，“.”表示缺失值；Attributes：以多个键值对组成的注释信息描述。','type'=>"part",'name'=>"表$tid. 新基因的GFF文件",'path'=>"HTML/template/new_gene.xls");
$writer->emptyTag('file','desc'=>"",'name'=>"新基因的GFF文件",'action'=>"xls",'path'=>"$mRNA_relative_dir/BMK_1_NewGene/$config{Project_key}.newGene_final.filtered.gff",'type'=>"xls");
$writer->emptyTag('file','desc'=>"",'name'=>"新基因序列FASTA文件",'action'=>"xls",'path'=>"$mRNA_relative_dir/BMK_1_NewGene/$config{Project_key}.newGene.longest_transcript.fa",'type'=>"xls");

$writer->emptyTag('h3','name'=>"2.$h.4 新基因功能注释",'type'=>'type1','desc'=>"2.$h.4 新基因功能注释");
$writer->emptyTag('p','desc'=>"使用BLAST软件将发掘的新基因与NR，Swiss-Prot，GO，COG，KEGG数据库进行序列比对，获得新基因的注释信息。",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>'注：Annotated databases：用于功能注释的数据库；New Gene Number：获得相应数据库注释信息的新基因数目。','type'=>"part",'name'=>"表$tid. 新基因功能注释结果统计",'path'=>"HTML/template/New_gene_anno.stat");

if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"新基因功能注释结果路径",'action'=>"xls",'path'=>"$mRNA_relative_dir/BMK_1_NewGene/BMK_1_NewGene_Anno/",'type'=>"xls");
}
	###DEG 分析
$writer->emptyTag('h3','name'=>"2.$h.5 基因表达定量、差异分析及功能注释",'type'=>'type1','desc'=>"2.$h.5 基因表达定量、差异分析及功能注释");
$writer->emptyTag('p','desc'=>"关于基因表达量的整体分布，差异表达筛选，及功能注释和富集分析详见以下网页:",'type'=>"type1");
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"基因差异表达及功能分析网页",'action'=>"xls",'path'=>"HTML/gene/gene.html",'type'=>"xls");
}else{
	$writer->emptyTag('file','desc'=>"",'name'=>"基因差异表达及功能分析网页",'action'=>"xls",'path'=>"HTML/gene.xml",'type'=>"xml");
}
$pid++;
$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. 差异及富集分析示意图",'type'=>"type1",'path'=>"HTML/template/1.geneDEG.png");
$writer->emptyTag('file','desc'=>"",'name'=>"基因表达量文件",'action'=>"xls",'path'=>"$mRNA_relative_dir/BMK_2_geneExpression/gene_expression.xls",'type'=>"xls");
	
	###lncRNA
$h++;
$writer->emptyTag('h2','name'=>"2.$h lncRNA分析",'type'=>'type1','desc'=>"2.$h lncRNA分析");
$writer->emptyTag('h3','name'=>"2.$h.1 lncRNA鉴定",'type'=>'type1','desc'=>"2.$h.1 lncRNA鉴定");
$writer->emptyTag('p','desc'=>"LncRNA (长链非编码RNA) 是一类转录本长度超过200nt，不能够编码蛋白的 RNA 分子。我们基于转录组拼接结果，首先与已知的lncRNA进行比较，去除组装出来的已知转录本（class_code为=,c），对余下的转录本进行lncRNA预测。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"新lncRNA的预测包含基本筛选和潜在编码能力筛选两部分。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(1) 基本筛选包含：长度≥200bp，Exon个数≥2，FPKM≥0.1以及与已知基因比较的class_code为\"i\",\"x\",\"u\",\"o\",\"e\"的转录本；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(2) 编码能力筛选进一步去掉基本筛选后中具有潜在编码能力的转录本，余下的即为新预测的lncRNA。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"因lncRNA不编码蛋白，因此，通过对转录本进行编码潜能筛选，判断其是否具有编码潜能，从而可以判定该转录本是否为lncRNA。百迈客综合目前应用最广泛的编码潜能分析方法对以上候选lncRNA进行进一步的筛选，主要包括：CPC<a href=\"#ref4\">[4]</a>分析、CNCI<a href=\"#ref5\">[5]</a>分析、CPAT<a href=\"#ref6\">[6]</a>分析、pfam<a href=\"#ref7\">[7]</a>蛋白结构域分析四种方法。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"将各软件所识别出的lncRNA转录本进行统计，采用维恩图直观的展示各个方法预测出的lncRNA转录本共有和特有的数目。取几种工具预测结果的交集，作为本次分析预测得到新的lncRNA 数据集进行后续的分析。",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"注：每个圆圈代表一种预测lncRNA的方法，圈中数字代表该方法所预测出lncRNA个数，最后项目取四款软件的交集作为预测结果。",'name'=>"图$pid. 预测方法维恩图",'type'=>"type1",'path'=>"$lnc_relative_dir/BMK_2_LncRNA_Prediction/BMK_1_Software_Result/venn.png");
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"LncRNA编码能力预测结果路径",'action'=>"xls",'path'=>"$lnc_relative_dir/BMK_2_LncRNA_Prediction/BMK_1_Software_Result/",'type'=>"xls");
}
$writer->emptyTag('p','desc'=>"经过四款软件预测得到的lncRNA其分类结果见柱形图：",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"注：横坐标为4种不同类型的LncRNA（lincRNA：基因间区长链非编码RNA；Antisense-lncRNA：反义长链非编码RNA；Intronic-lncRNA：内含子长链非编码RNA；sense_lncRNA：正义长链非编码RNA），纵坐标为对应的LncRNA数量。",'name'=>"图$pid. 预测长链非编码RNA统计图",'type'=>"type1",'path'=>"$lnc_relative_dir/BMK_2_LncRNA_Prediction/LncRNA_classification.png");

$writer->emptyTag('file','desc'=>"",'name'=>"LncRNA分类结果文件路径",'action'=>"xls",'path'=>"$lnc_relative_dir/BMK_2_LncRNA_Prediction/",'type'=>"xls");

$writer->emptyTag('h3','name'=>"2.$h.2 lncRNA表达定量",'type'=>'type1','desc'=>"2.$h.2 lncRNA表达定量");
$writer->emptyTag('p','desc'=>"对新预测的lncRNA和已知的lncRNA进行定量分析，得到每个样品的lncRNA表达量情况如下表所示：",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>'','type'=>"part",'name'=>"表$tid. lncRNA表达量结果文件",'path'=>"HTML/template/lncRNA_exp.xls");
$writer->emptyTag('file','desc'=>"",'name'=>"lncRNA表达量结果文件",'action'=>"xls",'path'=>"$lnc_relative_dir/BMK_3_LncRNA_Expression/lncRNA_expression.xls",'type'=>"xls");

$writer->emptyTag('h3','name'=>"2.$h.3 lncRNA靶基因预测",'type'=>'type1','desc'=>"2.$h.3 lncRNA靶基因预测");
$writer->emptyTag('p','desc'=>"第一种，lncRNA调控其邻近基因的表达，主要根据lncRNA与gene的位置关系预测，lncRNA 100kb范围内的邻近基因为其靶基因；",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>'','type'=>"part",'name'=>"表$tid. 基于位置关系的LncRNA靶基因预测结果",'path'=>"HTML/template/cis.xls");
$writer->emptyTag('file','desc'=>"",'name'=>"LncRNA顺式靶基因",'action'=>"xls",'path'=>"$lnc_relative_dir/BMK_4_LncRNA_Target/Cis_target_gene.xls",'type'=>"xls");
$writer->emptyTag('p','desc'=>"第二种，通过样本间 lncRNA 与 mRNA 的表达量相关性分析方法来预测 lncRNA的反式靶基因。此方法需要一定样品数量保证其预测的准确性，因此当样本量小于5时不进行此项分析，结果为空。采用 Pearson 相关系数法分析样本间 lncRNA 与 mRNA 的相关性，取相关性绝对值大于$config{Trans_cor}且显著性p值小于$config{Trans_p}的基因作为lncRNA的反式靶基因。",'type'=>"type1");
if(-e "$inpath/HTML/template/trans.xls"){
	$tid++;
	$writer->emptyTag('table','desc'=>'','type'=>"part",'name'=>"表$tid. 基于表达相关性的LncRNA靶基因预测结果",'path'=>"HTML/template/trans.xls");
	$writer->emptyTag('file','desc'=>"",'name'=>"LncRNA反式靶基因",'action'=>"xls",'path'=>"$lnc_relative_dir/BMK_4_LncRNA_Target/Trans_target_gene.xls",'type'=>"xls");
}

	###lncRNA差异分析
$writer->emptyTag('h3','name'=>"2.$h.4 lncRNA差异表达及功能分析",'type'=>'type1','desc'=>"2.$h.4 lncRNA差异表达及功能分析");
$writer->emptyTag('p','desc'=>"关于lncRNA表达量的整体分布，差异表达lncRNA筛选，及差异表达lncRNA的顺/反式靶基因功能注释和富集分析可详见以下网页:",'type'=>"type1");
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"lncRNA差异表达及功能分析网页",'action'=>"xls",'path'=>"HTML/lncRNA/lncRNA.html",'type'=>"xls");
}else{
	$writer->emptyTag('file','desc'=>"",'name'=>"lncRNA差异表达及功能分析网页",'action'=>"xls",'path'=>"HTML/lncRNA.xml",'type'=>"xml");
}
$pid++;
$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. 差异lncRNA及靶基因功能富集分析示意图",'type'=>"type1",'path'=>"HTML/template/2.noncodeDEG.png");

	###lncRNA disease分析
if(-e "$inpath/$lnc_relative_dir/BMK_6_Lnc_disease/LncRNA_disease.xls"){
	$writer->emptyTag('h3','name'=>"2.$h.5 lncRNA疾病注释",'type'=>'type1','desc'=>"2.$h.5 lncRNA疾病注释");
	$writer->emptyTag('p','desc'=>"虽然目前对于lncRNA的了解有限，但新兴研究显示，lncRNA在很多的生物过程中发挥关键作用，并与许多疾病相关，如癌症，心血管疾病和神经变性疾病。lncRNA对于解读生命科学是至关重要的，尤其是癌症。LncRNADisease<a href=\"#ref8\">[8]</a>数据库中存储了实验验证的与疾病相关的一些lncRNA。对项目中的lncRNA，采用序列比对的方式判定是否与疾病相关。结果如下：",'type'=>"type1");
	$tid++;
	$writer->emptyTag('table','desc'=>'注：#LncRNA_id，项目中lncRNA ID；DB_lncRNA，数据库中比对上的lncRNA；Source，lncRNA来源；disease：相关的疾病；regulate：调控方式；description：疾病相关描述；chr：染色体；start：起始；end：终止；strand：正负链；pmid：文章。','type'=>"part",'name'=>"表$tid. lncRNA疾病注释",'path'=>"HTML/template/lncRNA_disease.xls");
	$writer->emptyTag('file','desc'=>"",'name'=>"lncRNA疾病注释信息",'action'=>"xls",'path'=>"$lnc_relative_dir/BMK_6_Lnc_disease/LncRNA_disease.xls",'type'=>"xls");
}
	###lncRNA 前体分析
if(-e "$inpath/$lnc_relative_dir/BMK_7_Lnc_precursor/lncRNA_precursor.xls"){
        $writer->emptyTag('h3','name'=>"2.$h.6 lncRNA中miRNA前体分析",'type'=>'type1','desc'=>"2.$h.6 lncRNA中miRNA前体分析");
        $writer->emptyTag('p','desc'=>"lncRNA可以作为miRNA的前体分子，利用miRbase<a href=\"#ref9\">[9]</a>数据库中的miRNA序列与lncRNA比对预测，统计作为miRNA前体的lncRNA。结果见下表：",'type'=>"type1");
	$tid++;
	$writer->emptyTag('table','desc'=>"注：#LncRNA_ID：LncRNA ID；miRNA_hp：miRNA前体ID；LncRNA_start：LncRNA比对起始位置；LncRNA_end：LncRNA比对终止位置。",'type'=>"part",'name'=>"表$tid. lncRNA前体注释",'path'=>"HTML/template/lncRNA_precursor.xls");
        $writer->emptyTag('file','desc'=>"",'name'=>"miRNA前体的lncRNA统计",'action'=>"xls",'path'=>"$lnc_relative_dir/BMK_7_Lnc_precursor/lncRNA_precursor.xls",'type'=>"xls");
	}
}
if(defined $circ_dir){
	$h++;
	my ($path1,$path2,$path3,$path4,$Disc1);
	if(defined $lnc_dir){
		$path1 = "BMK_1_rawData/BMK_2_Mapped_Statistics/circRNA";
		$Disc1 = "使用指定的基因组作为参考进行序列比对及后续分析";
	}else{
		my $ref_link = "<a href=\"$config{ref_gene_addr}\" title=\"click\" target=\"_blank\">$config{ref_gene_addr}</a>";
		$path1 = "BMK_1_rawData/BMK_2_Mapped_Statistics";
		$Disc1 = "使用指定的基因组$config{ref_gene_name}作为参考进行序列比对及后续分析。$config{ref_gene_name} 下载地址：$ref_link";
	}

	####circRNA
	$writer->emptyTag('h2','name'=>"2.$h circRNA分析",'type'=>'type1','desc'=>"2.$h circRNA分析");
	$writer->emptyTag('p','desc'=>"$Disc1。",'type'=>"type1");
	$writer->emptyTag('p','desc'=>"获得Clean Reads后，将其与参考基因组进行序列比对，获取在参考基因组或基因上的位置信息，以及测序样品特有的序列特征信息。使用BWA<a href=\"#ref10\">[10]</a>进行比对，我们将比对到指定的参考基因组上的reads称为Mapped Reads，对应的数据称为Mapped Data。",'type'=>"type1");
	$writer->emptyTag('h3','name'=>"2.$h.1 比对结果统计及作图",'type'=>'type1','desc'=>"2.$h.1 比对结果统计及作图");
	if(!defined $lnc_dir){
		$writer->emptyTag('p','desc'=>"比对效率指Mapped Reads占Clean Reads的百分比，是数据利用率的最直接体现。比对效率除了受数据测序质量影响外，还与指定的参考基因组组装的优劣、参考基因组与测序样品的生物学分类关系远近（亚种）有关。因此，通过比对效率，可以评估所选参考基因组组装是否能满足信息分析的需求，及后期数据分析的可靠性。",'type'=>"type1");
	}
	$writer->emptyTag('p','desc'=>"各样品测序数据与所选参考基因组的序列比对结果统计见下表：",'type'=>"type1");
	$tid++;
	$writer->emptyTag('table','desc'=>"注：BMK-ID：百迈客对样品的统一编号；Total Reads：Clean Reads数目，按单端计；Mapped Reads：比对到参考基因组上的Reads数目及在Clean Reads中占的百分比；Uniq Mapped Reads：比对到参考基因组唯一位置的Reads数目及在Clean Reads中占的百分比；Multiple Map Reads：比对到参考基因组多处位置的Reads数目及在Clean Reads中占的百分比；Reads Map to '+'：比对到参考基因组正链的Reads数目及在Clean Reads中占的百分比；Reads Map to '-'：比对到参考基因组负链的Reads数目及在Clean Reads中占的百分比。",'type'=>"part",'name'=>"表$tid. Clean Data与参考基因组比对结果统计表",'path'=>"$path1/All.mappedStat.xls");

	my $hh3=1;
	if(!defined $lnc_dir){
		$writer->emptyTag('p','desc'=>"将比对到不同染色体上Reads进行位置分布统计，绘制Mapped Reads在所选参考基因组上的覆盖深度分布图。各样品的Mapped Reads在参考基因组上的覆盖深度分布图如下：",'type'=>"type1");
		&piclist("Mapped Reads在参考基因组上的位置及覆盖深度分布图","注：横坐标为染色体位置；纵坐标为覆盖深度以2为底的对数值，以10kb作为区间单位长度，划分染色体成多个小窗口（Window），统计落在各个窗口内的Mapped Reads作为其覆盖深度。蓝色为正链，绿色为负链。","$inpath/BMK_1_rawData/BMK_2_Mapped_Statistics/*.readDensity.png","BMK_1_rawData");
		$writer->emptyTag('p','desc'=>"统计Mapped Reads在指定的参考基因组不同区域（外显子、内含子和基因间区）的数目，绘制基因组不同区域上各样品Mapped Reads的分布直方图，如下：",'type'=>"type1");
		$pid++;
		$writer->emptyTag('pic','desc'=>"注：图中每个直方柱表示一个样品，粉色区域为外显子区、绿色区域为基因间区、蓝色区域为内含子区，区域的高度表示比对到该区域的Mapped Reads在所有Mapped Reads中所占的百分比。",'name'=>"图$pid. 基因组不同区域Reads分布直方图",'type'=>"type1",'path'=>"BMK_1_rawData/BMK_2_Mapped_Statistics/Total_Type.png");
		if(!defined $cloud){
			$writer->emptyTag('file','desc'=>"",'name'=>"circRNA比对统计结果路径",'action'=>"xls",'path'=>"BMK_1_rawData/BMK_2_Mapped_Statistics/",'type'=>"xls");
		}
		$hh3++;
		$writer->emptyTag('h3','name'=>"2.$h.$hh3 文库质量评估",'type'=>'type1','desc'=>"2.$h.$hh3 文库质量评估");
		$writer->emptyTag('p','desc'=>"合格的测序文库是测序的必要条件，为确保文库的质量，需要对测序文库进行质量评估。",'type'=>"type1");
		$writer->emptyTag('p','desc'=>"我们通过插入片段的长度分布，评估插入片段长度的离散程度。插入片段长度的离散程度能直接反映出文库制备过程中磁珠纯化的效果。通过插入片段两端的Reads在参考基因组上的比对起止点之间的距离计算插入片段长度。",'type'=>"type1");
		$writer->emptyTag('p','desc'=>"大部分的真核生物基因为断裂基因，外显子被内含子隔断，而转录组测序得到的是无内含子的成熟CircRNA。当CircRNA中跨内含子的片段两端的Reads比对到基因组上时，比对起止点之间的距离要大于插入片段长度。因此，在插入片段长度模拟分布图中，主峰右侧形成1个或多个次峰。",'type'=>"type1");
		$writer->emptyTag('p','desc'=>"各样品的插入片段长度模拟分布图如下：",'type'=>"type1");
		&piclist("插入片段长度模拟分布图","注：横坐标为双端Reads在参考基因组上的比对起止点之间的距离，范围为0到800bp；纵坐标为比对起止点之间不同距离的双端Reads或插入片段数量。","$inpath/BMK_1_rawData/BMK_3_Library_Assessment/*.insertSize.png","BMK_1_rawData");
		if(!defined $cloud){
			$writer->emptyTag('file','desc'=>"",'name'=>"circRNA文库质量评估结果路径",'action'=>"xls",'path'=>"BMK_1_rawData/BMK_3_Library_Assessment/",'type'=>"xls");
		}
	}

	$hh3++;
	$writer->emptyTag('h3','name'=>"2.$h.$hh3 circRNA预测",'type'=>'type1','desc'=>"2.$h.$hh3 circRNA预测");
	$writer->emptyTag('p','desc'=>"鉴于circRNA鉴定存在较高的假阳性率，采用两款软件鉴定circRNA，取二者的交集。 两款软件鉴定原理如下：",'type'=>"type1");
	$writer->emptyTag('p','desc'=>"CIRI预测环状：CIRI<a href=\"#ref11\">[11]</a>软件使用BWA软件与参考基因序列比对，生成SAM文件，并对SAM文件中的CIGAR值进行分析的，从SAM文件中扫描PCC信号（paired chiastic clipping signals）。 CIGAR值在junction read的特征是xS/HyM或者xMyS/H，其中x,y代表碱基数目，M是mapping上的，S是soft clipping，H是hard clipping。对于双端Reads，CIRI算法会考虑一对reads，其中一条可以mapping到CircRNA上，另一条也需要mapping到CircRNA的区间内。对于单外显子成环，或者“长外显子1-短外显子-长外显子2”形成的环形结构，CIGAR值应该是xS/HyMzS/H以及(x+y)S/HzM或者xM(y+z)S/H，CIRI软件可以很好的将这两种情况分开。对于splicing 信号(GT,AG) CIRI也会考虑其他弱splicing 信息例如（AT-AC），算法会从GTF/GFF文件中抽取外显子边界位置，并用已知的边界来过滤假阳性。",'type'=>"type1");
	if(!defined $cloud){
		$writer->emptyTag('file','desc'=>"",'name'=>"CIRI软件预测结果路径",'action'=>"xls",'path'=>"$circ_relative_dir/BMK_1_circRNA_Prediction/BMK_1_CIRC_info/BMK_1_CIRI/",'type'=>"xls");
	}
	$writer->emptyTag('p','desc'=>"find_circ预测环状：由于环状RNA成环的剪接位点不能直接比对到基因组上，find_circ<a href=\"#ref12\">[12]</a>软件首先将不能和基因组比对上的reads两端各取20bp作为锚点，再将锚点作为独立的reads往基因组上比对并寻找唯一匹配位点，如果两个锚点的比对位置在线性上方向呈反向，就延长锚点的reads，直至找到环状RNA的接合位置，若此时两侧的序列分别为GT/AG剪接信号，则判断为环状RNA。候选环状RNA必须满足以下条件：",'type'=>"type1");
	$writer->emptyTag('p','desc'=>"1.GU/AG 在剪接位点的两侧出现；",'type'=>"type1");
	$writer->emptyTag('p','desc'=>"2.可以检测到清晰的breakpoint；",'type'=>"type1");
	$writer->emptyTag('p','desc'=>"3.只支持2个mismatch；",'type'=>"type1");
	$writer->emptyTag('p','desc'=>"4.breakpoint不能在anchor 2个核苷酸之外的地方出现，也就是最多2nt；",'type'=>"type1");
	$writer->emptyTag('p','desc'=>"5.至少有两条reads支持这个junction；",'type'=>"type1");
	$writer->emptyTag('p','desc'=>"6.比对到正确的一个短序列的位置要比它比对到其他位置的分值高35分以上。",'type'=>"type1");
	if(!defined $cloud){
		$writer->emptyTag('file','desc'=>"",'name'=>"find_circ软件预测结果路径",'action'=>"xls",'path'=>"$circ_relative_dir/BMK_1_circRNA_Prediction/BMK_1_CIRC_info/BMK_2_find_circ/",'type'=>"xls");
	}
	$writer->emptyTag('p','desc'=>"对两款软件预测的结果取交集，共得到 $info{circRNA_num} 个circRNA。根据circRNA的来源基因对circRNA进行重命名，并整合circRNA来源基因的位置信息、注释信息，最终得到的结果如下表所示。",'type'=>"type1");
	$tid++;
	$writer->emptyTag('table','desc'=>"注：New_ID：根据circRNA来源基因对环状RNA的重新命名；circRNA_ID：本项目中所使用的circRNA ID；circRNA_chr：circRNA所在染色体；circRNA_start：circRNA开始位置；circRNA_end：cicrRNA结束位置；Source_gene_id：circRNA来源基因ID；Source_gene_start：来源基因开始位置；Source_gene_end：来源基因结束位置；之后几列分别为circRNA的在每个样品中的junction read数，geneLength为环状长度，之后是circRNA来源基因注释结果。",'type'=>"part",'name'=>"表$tid. 环状RNA预测结果",'path'=>"HTML/template/circRNA_new.xls");

	$writer->emptyTag('file','desc'=>"",'name'=>"CircRNA预测结果",'action'=>"xls",'path'=>"$circ_relative_dir/BMK_1_circRNA_Prediction/circRNA_newname.xls",'type'=>"xls");
	
	$writer->emptyTag('p','desc'=>"各样品预测到的circRNA的数目如下表所示：",'type'=>"type1");
	$tid++;
	$writer->emptyTag('table','desc'=>'注：Total_reads：最后得到的junction reads总数；Total_circRNA：最后得到的CircRNA数目。','type'=>"part",'name'=>"表$tid. circRNA数目统计表",'path'=>"$circ_relative_dir/BMK_1_circRNA_Prediction/BMK_2_statistics/statistic.xls");
	$hh3++;
	$writer->emptyTag('h3','name'=>"2.$h.$hh3 circRNA来源统计及长度和表达量分布",'type'=>'type1','desc'=>"2.$h.$hh3 circRNA来源统计及长度和表达量分布");
	$writer->emptyTag('p','desc'=>"各样品的CircRNA 在基因组上的分布情况，结果见下图：",'type'=>"type1");
	&piclist("每个样品的CircRNA来源分布图","","$inpath/$circ_relative_dir/BMK_1_circRNA_Prediction/BMK_2_statistics/*.circRNA.type.png","$circ_relative_dir");
	$writer->emptyTag('p','desc'=>"统计各样品的CircRNA在不同染色体上的reads分布，结果见下图：",'type'=>"type1");
	&piclist("每个样品的CircRNA染色体reads分布图","注：横坐标表示染色体，纵坐标表示对应染色体上的环状RNA junction reads数量。","$inpath/$circ_relative_dir/BMK_1_circRNA_Prediction/BMK_2_statistics/*.chromosome_distrbution.png","$circ_relative_dir");
	
	$writer->emptyTag('p','desc'=>"对预测得到的所有 CircRNA 长度进行分布统计，如下图所示：",'type'=>"type1");
	$pid++;
	$writer->emptyTag('pic','desc'=>"注：横坐标代表circRNA的长度区间；纵坐标代表该长度区间内的环状RNA数目。",'name'=>"图$pid. CircRNA长度分布图",'type'=>"type1",'path'=>"$circ_relative_dir/BMK_1_circRNA_Prediction/BMK_2_statistics/All_circRNA.circRNA.length.distribution.png");
	
	if(-e "$inpath/$circ_relative_dir/BMK_1_circRNA_Prediction/BMK_2_statistics/circRNA_circos.png"){
		$writer->emptyTag('p','desc'=>"用$config{normalization}算法对circRNA表达量进行归一化处理，具体公式可见“circRNA差异表达及功能分析”网页。circRNA在基因组的位置和表达量分布如下：",'type'=>"type1");
		$pid++;
		$writer->emptyTag('pic','desc'=>"注：CircRNA分布在基因组的统计以及各个样品表达量的circos图。其中，柱状图表示各个位置的CircRNA的个数，热图表示每个样品的表达量分布，0~1：蓝色，1~10：绿色，10~100：橘色，100～：红色，热图由外向内依次表示T01、T02、T03等。",'name'=>"图$pid. circRNA表达量circos图",'type'=>"type1",'path'=>"$circ_relative_dir/BMK_1_circRNA_Prediction/BMK_2_statistics/circRNA_circos.png");
	}
	if(!defined $cloud){
		$writer->emptyTag('file','desc'=>"",'name'=>"预测得到的circRNA相关统计结果路径",'action'=>"xls",'path'=>"$circ_relative_dir/BMK_1_circRNA_Prediction/BMK_2_statistics/",'type'=>"xls");
	}
	######circRNA 差异分析
	$hh3++;
	$writer->emptyTag('h3','name'=>"2.$h.$hh3 circRNA差异表达及功能分析",'type'=>'type1','desc'=>"2.$h.$hh3 circRNA差异表达及功能分析");
	$writer->emptyTag('p','desc'=>"关于circRNA表达量的整体分布，差异表达circRNA筛选，及差异表达circRNA其来源基因的功能注释和富集分析详见以下网页：",'type'=>"type1");
	if(!defined $cloud){
		$writer->emptyTag('file','desc'=>"",'name'=>"circRNA差异表达及功能分析网页",'action'=>"xls",'path'=>"HTML/circRNA/circRNA.html",'type'=>"xls");
	}else{
		$writer->emptyTag('file','desc'=>"",'name'=>"circRNA差异表达及功能分析网页",'action'=>"xls",'path'=>"HTML/circRNA.xml",'type'=>"xml");
	}
	$pid++;
	$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. 差异circRNA及来源基因功能富集分析示意图",'type'=>"type1",'path'=>"HTML/template/2.noncodeDEG.png");
	$writer->emptyTag('file','desc'=>"",'name'=>"circRNA表达量文件",'action'=>"xls",'path'=>"$circ_relative_dir/BMK_2_circRNA_Expression/circRNA_expression.xls",'type'=>"xls");

	if(-e "$inpath/$circ_relative_dir/BMK_1_circRNA_Prediction/BMK_3_Known/circRNA_known.xls"){
		$hh3++;
		$writer->emptyTag('h3','name'=>"2.$h.$hh3 已知circRNA预测",'type'=>'type1','desc'=>"2.$h.$hh3 已知circRNA预测");
		$writer->emptyTag('p','desc'=>"采用序列比对方法与circBase<a href=\"#ref13\">[13]</a>数据库中的环状序列进行比较，如果项目中预测出的circRNA与数据库中的circRNA的相似度超过90%，则认为其是已知的环状RNA。",'type'=>"type1");
		$pid++;
		$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. 已知环状RNA和预测到的环状RNA的统计图",'type'=>"type1",'path'=>"$circ_relative_dir/BMK_1_circRNA_Prediction/BMK_3_Known/Known_Unknown_pie.png");
		$writer->emptyTag('file','desc'=>"注：第一列是项目中鉴定出的circRNA，第二列是已知的circRNA，第三列是比对上的长度占已知环状长度的比例。",'name'=>"已知circRNA关系表",'action'=>"xls",'path'=>"$circ_relative_dir/BMK_1_circRNA_Prediction/BMK_3_Known/circRNA_known.xls",'type'=>"xls");

		if(-e "$inpath/$circ_relative_dir/BMK_1_circRNA_Prediction/BMK_3_Known/circRNA_Circ2disease_annotation.xls"){
			$hh3++;
	                $writer->emptyTag('h3','name'=>"2.$h.$hh3 已知circRNA注释",'type'=>'type1','desc'=>"2.$h.$hh3 已知circRNA注释");
        	        $writer->emptyTag('p','desc'=>"<a href=\"http://bioinformatics.zju.edu.cn/Circ2Disease/index.html\" title=\"click\" target=\"_blank\">Circ2disease</a>中收录了与疾病相关的人类circRNA的信息，且这些与疾病的相关关系是经过实验验证的，此外数据库还预测了这些经过实验验证的circRNA对应的miRNA sponges并整合了CircInteratome数据库中的RBP（RNA结合蛋白）信息。利用数据库中整合的信息，将预测得到的已知circRNA在Circ2disease中进行注释，结果如下表：",'type'=>"type1");
                	$tid++;
	                $writer->emptyTag('table','desc'=>"注：#circRNA_id：预测得到的circRNA id；circBase_id：circBase数据库中circRNA id；circRNA_length：预测得到的circRNA序列长度；circBase_length：circBase数据库中circRNA序列长度；circRNA_name：Circ2disease数据库中circRNA name；Synonyms：Circ2disease数据库中circRNA别名；Genome_Version：Circ2disease数据库中circRNA基因组版本；Gene_Symbol：Circ2disease数据库中circRNA Hostgene Symbol；Type：Circ2disease数据库中circRNA type；Disease_Name：Circ2disease数据库中circRNA 经实验验证的与circRNA相关疾病名称；Methods：验证方法；Sample：使用样品；Location：细胞中位置；Expression_Pattern：表达模式(上/下调)；miRNA_Validated：实验验证相关miRNA；miRNA_Target_Validated：实验验证相关miRNA的靶基因；RBP_Validated：实验验证相关RNA结合蛋白；RBP_Target_Validated：实验验证相关RNA结合蛋白靶；Functional_Description：功能描述；PubMed_ID：来源文章PubMed ID；Year：来源文章发表时间；Journal：来源文章发表期刊；Title：来源文章Title；Note：其他信息。",'type'=>"part",'name'=>"表$tid. 已知circRNA的Circ2disease数据库注释",'path'=>"HTML/template/circ2disease.xls");
        	        $writer->emptyTag('file','desc'=>"",'name'=>"已知circRNA的Circ2disease注释文件",'action'=>"xls",'path'=>"$circ_relative_dir/BMK_1_circRNA_Prediction/BMK_3_Known/circRNA_Circ2disease_annotation.xls",'type'=>"xls");
                	$writer->emptyTag('p','desc'=>"<a href=\"http://www.circbank.cn/index.html\" title=\"click\" target=\"_blank\">CircBank</a>是一款整合的人类circRNA数据库，数据库共收录了140,790条人类circRNA的记录，每一条circRNA记录都单独做了一个详细信息的页面。针对每个circRNA的信息主要包括：该circRNA的详细序列；在小鼠中同源性较高的circRNA及其对应的序列；miRNA结合的预测分析；ORF预测分析；COSMIC记录的突变和多态性位点汇总；m6A修饰信息。此外circBank还专门开发了一套专用的ID号，致力于规范化circRNA命名。利用数据库中整合的信息，将预测得到的已知circRNA在circBank中进行注释，结果如下表：",'type'=>"type1");
	                $tid++;
        	        $writer->emptyTag('table','desc'=>"#circRNA_id：预测得到的circRNA id；circBase_id：circBase数据库中circRNA id；circRNA_length：预测得到的circRNA序列长度；circBase_length：circBase数据库中circRNA序列长度；circBank_id：已知circRNA对应circBank数据库id（来自circBank）；annotation：circRNA注释（来自circBank）；best_transcript：最优转录本（来自circBank）；gene_symbol：circRNA hostgene Symbol（来自circBank）；mm9_circRNA_id：conserved mm9 circRNA（circBase id）；mm9_circbase_info：circBase中mm9 circRNA详细信息；mm9_circbase_sequence：circBase mm9 circRNA序列；orf_size：ORF大小；fickett_score：fickett打分；hexamer_score：hexamer打分；coding_prob：circRNA编码能力；m6A_name：circRNA相关m6A名称。",'type'=>"part",'name'=>"表$tid. 已知circRNA的circBank数据库注释",'path'=>"HTML/template/circbank.xls");
                	$writer->emptyTag('file','desc'=>"",'name'=>"已知circRNA的circBank注释文件",'action'=>"xls",'path'=>"$circ_relative_dir/BMK_1_circRNA_Prediction/BMK_3_Known/circRNA_circBank_annotation.xls",'type'=>"xls");
        	}
	}
}
######################结构分析
if(defined $lnc_dir && defined $circ_dir){
	$h++;
	$writer->emptyTag('h2','name'=>"2.$h miRNA靶向关系预测",'type'=>'type1','desc'=>"2.$h miRNA靶向关系预测");
	$writer->emptyTag('p','desc'=>"根据已知miRNA与对应物种的基因序列信息，植物用 TargetFinder软件<a href=\"#ref16\">[16]</a>进行靶基因预测；动物用miRanda<a href=\"#ref17\">[17]</a>和targetscan<a href=\"#ref18\">[18]</a>进行靶基因预测。miRNA靶基因预测结果如下表：",'type'=>"type1");
	$tid++;
	$writer->emptyTag('table','desc'=>"注：0表示该靶向关系对在该软件中并未预测到，1表示该靶向关系对在该软件中预测到，最终取两款软件的交集作为最终靶向关系对。",'type'=>"part",'name'=>"表$tid. miRNA靶基因预测结果",'path'=>"HTML/template/miRNA_gene.xls");
	$writer->emptyTag('file','desc'=>"注：第一列是miRNA,第二列是基因，如果靶向多个基因，以逗号分割。",'name'=>"miRNA靶基因关系表",'action'=>"xls",'path'=>"$struc_relative_dir/BMK_2_Target_Predict/miRNA-mRNA/mir2target.list.xls",'type'=>"xls");
	$writer->emptyTag('p','desc'=>"根据已知miRNA和新预测的miRNA与lncRNA序列信息，软件同上。miRNA靶向lncRNA预测结果如下表：",'type'=>"type1");
	$tid++;
	$writer->emptyTag('table','desc'=>'注：0表示该靶向关系对在该软件中并未预测到，1表示该靶向关系对在该软件中预测到，最终取两款软件的交集作为最终靶向关系对。','type'=>"part",'name'=>"表$tid. miRNA靶向lncRNA预测结果",'path'=>"HTML/template/miRNA_lncRNA.xls");
	$writer->emptyTag('file','desc'=>"注：第一列是miRNA，第二列是lncRNA，如果靶向多个lncRNA，以逗号分割。",'name'=>"miRNA靶向lncRNA关系表",'action'=>"xls",'path'=>"$struc_relative_dir/BMK_2_Target_Predict/miRNA-lncRNA/mir2target.list.xls",'type'=>"xls");
	$writer->emptyTag('p','desc'=>"根据已知miRNA和新预测的miRNA与对应物种的circRNA序列信息，软件同上。miRNA靶circRNA预测结果如下表：",'type'=>"type1");
	$tid++;
	$writer->emptyTag('table','desc'=>'注：0表示该靶向关系对在该软件中并未预测到，1表示该靶向关系对在该软件中预测到，最终取两款软件的交集作为最终靶向关系对。','type'=>"part",'name'=>"表$tid. miRNA靶向circRNA预测结果",'path'=>"HTML/template/miRNA_circRNA.xls");
	$writer->emptyTag('file','desc'=>"注：第一列是miRNA，第二列是circRNA，如果靶向多个circRNA，以逗号分割。",'name'=>"miRNA靶向circRNA关系表",'action'=>"xls",'path'=>"$struc_relative_dir/BMK_2_Target_Predict/miRNA-circRNA/mir2target.list.xls",'type'=>"xls");
	if(!defined $cloud){
		$writer->emptyTag('file','desc'=>"",'name'=>"miRNA靶向关系结果路径",'action'=>"xls",'path'=>"$struc_relative_dir/BMK_2_Target_Predict/",'type'=>"xls");
	}
}elsif(defined $lnc_dir && !defined $circ_dir){
	$h++;
        $writer->emptyTag('h2','name'=>"2.$h miRNA靶向关系预测",'type'=>'type1','desc'=>"2.$h miRNA靶向关系预测");
        $writer->emptyTag('p','desc'=>"根据已知miRNA与对应物种的基因序列信息，植物用 TargetFinder软件<a href=\"#ref16\">[16]</a>进行靶基因预测；动物用miRanda<a href=\"#ref17\">[17]</a>和targetscan<a href=\"#ref18\">[18]</a>进行>靶基因预测。miRNA靶基因预测结果如下表：",'type'=>"type1");
        $tid++;
        $writer->emptyTag('table','desc'=>'注：0表示该靶向关系对在该软件中并未预测到，1表示该靶向关系对在该软件中预测到，最终取两款软件的交集作为最终靶向关系对。','type'=>"part",'name'=>"表$tid. miRNA靶基因预测结果",'path'=>"HTML/template/miRNA_gene.xls");
        $writer->emptyTag('file','desc'=>"注：第一列是miRNA,第二列是基因，如果靶向多个基因，以逗号分割。",'name'=>"miRNA靶基因关系表",'action'=>"xls",'path'=>"$struc_relative_dir/BMK_2_Target_Predict/miRNA-mRNA/mir2target.list.xls",'type'=>"xls");
        $writer->emptyTag('p','desc'=>"根据已知miRNA与lncRNA序列信息，软件同上。miRNA靶向lncRNA预测结果如下表：",'type'=>"type1");
        $tid++;
        $writer->emptyTag('table','desc'=>'注：0表示该靶向关系对在该软件中并未预测到，1表示该靶向关系对在该软件中预测到，最终取两款软件的交集作为最终靶向关系对。','type'=>"part",'name'=>"表$tid. miRNA靶向lncRNA预测结果",'path'=>"HTML/template/miRNA_lncRNA.xls");
        $writer->emptyTag('file','desc'=>"注：第一列是miRNA，第二列是lncRNA，如果靶向多个lncRNA，以逗号分割。",'name'=>"miRNA靶向lncRNA关系表",'action'=>"xls",'path'=>"$struc_relative_dir/BMK_2_Target_Predict/miRNA-lncRNA/mir2target.list.xls",'type'=>"xls");
	if(!defined $cloud){
		$writer->emptyTag('file','desc'=>"",'name'=>"miRNA靶向关系结果路径",'action'=>"xls",'path'=>"$struc_relative_dir/BMK_2_Target_Predict/",'type'=>"xls");
	}
}elsif(!defined $lnc_dir && defined $circ_dir){
	$h++;
	$writer->emptyTag('h2','name'=>"2.$h circRNA-miRNA靶向关系预测",'type'=>'type1','desc'=>"2.$h circRNA-miRNA靶向关系预测");
	$writer->emptyTag('p','desc'=>"miRNA具有通过直接的碱基互补配结合靶向mRNA的非编码区，进而影响其翻译的功能。",'type'=>"type1");
	$writer->emptyTag('p','desc'=>"circRNA分子富含microRNA（miRNA）结合位点，在细胞中起到miRNA海绵（miRNA sponge）的作用，进而解除miRNA对其靶基因的抑制作用，升高靶基因的表达水平。这一作用机制被称为竞争性内源RNA（ceRNA）机制。通过与疾病关联的miRNA相互作用，circRNA在疾病中发挥着重要的调控作用。",'type'=>"type1");
	$writer->emptyTag('p','desc'=>"根据已知miRNA与对应物种的基因序列信息，植物用TargetFinder软件<a href=\"#ref16\">[16]</a>进行靶基因预测；动物用miRanda<a href=\"#ref17\">[17]</a>和targetscan<a href=\"#ref18\">[18]</a>进行靶基因预测。miRNA靶基因预测结果如下表：",'type'=>"type1");
	$tid++;
        $writer->emptyTag('table','desc'=>'注：0表示该靶向关系对在该软件中并未预测到，1表示该靶向关系对在该软件中预测到，最终取两款软件的交集作为最终靶向关系对。','type'=>"part",'name'=>"表$tid. miRNA靶基因预测结果",'path'=>"HTML/template/miRNA_circRNA.xls");
        $writer->emptyTag('file','desc'=>"注：第一列是miRNA,第二列是基因，如果靶向多个基因，以逗号分割。",'name'=>"miRNA靶基因关系表",'action'=>"xls",'path'=>"BMK_3_Target_Predict/miRNA-circRNA/mir2target.list.xls",'type'=>"xls");
}



$h++;
$writer->emptyTag('h2','name'=>"2.$h 结构分析",'type'=>'type1','desc'=>"2.$h 结构分析");
$H3=0;
if(defined $lnc_dir){
$H3++;
$writer->emptyTag('h3','name'=>"2.$h.$H3 SNP/InDel分析",'type'=>'type1','desc'=>"2.$h.$H3 SNP/InDel分析");
$writer->emptyTag('p','desc'=>"SNP（Single Nucleotide Polymorphisms）是指在基因组上由单个核苷酸变异形成的遗传标记，其数量很多，多态性丰富。百迈客基于各样品reads与参考基因组序列的比对结果，使用GATK<a href=\"#ref20\">[20]</a>软件识别测序样品与参考基因组间的单碱基错配，识别潜在的SNP位点。进而可以分析这些SNP位点是否影响了基因的表达水平或者蛋白产物的种类。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"InDel(insertion-deletion)是指相对于参考基因组，样本中发生的小片段的插入缺失，该插入缺失可能含一个或多个碱基。GATK也能够检测样品的插入缺失（InDel）。InDel变异一般比SNP变异少，同样反映了样品与参考基因组之间的差异，并且编码区的InDel会引起移码突变，导致基因功能上的变化。GATK识别标准如下：",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(1) 35bp范围内连续出现的单碱基错配不超过3个；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(2) 经过序列深度标准化的SNP质量值大于2.0。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"各样品分别按照以上条件筛选，最终获得可靠的SNP位点。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"SnpEff<a href=\"#ref21\">[21]</a>是一款用于注释变异（SNP、InDel）和预测变异影响的软件。根据变异位点在参考基因组上的位置以及参考基因组上的基因位置信息，可以得到变异位点在基因组发生的区域（基因间区、基因区或CDS区等），以及变异产生的影响（同义非同义突变等）。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"由于转录完成之后，mRNA除了需要加帽、加Ploy(A)和可变剪接之外，较少mRNA会经历RNA编辑（RNA editing），从而会产生单碱基的替换、插入、缺失。RNA编辑能使同一基因产生序列多样的mRNA，但是这种多态性不是基因组固有的多态性。从比对结果来看，SNP和单碱基替换的RNA编辑结果是一样的。因此，通过转录组测序数据识别出SNP不免会含有RNA编辑的产物。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"SNP/InDel位点信息如下：",'type'=>"type1");
&filelist("注：Chr：SNP/InDel位点所在染色体编号；Pos：SNP/InDel位点在染色体上的位置；Gene_id：SNP/InDel位点所在的基因或原来未注释的基因区（表中用Intergenic表示）；Ref：所选参考基因组中的SNP/InDel等位；Alt：测序样品中识别到的其他的SNP/InDel等位；L**：样品L**中SNP/InDel位点的分型；Depth：样品L**中SNP/InDel位点的测序深度；AlleDp：样品L**中SNP/InDel位点的各等位测序深度；Effect：SNP/InDel所在区域或类型；Codon_change：编码改变方式，未改变用点表示。核酸编码表见附表3，Effect具体说明详见：http://snpeff.sourceforge.net/SnpEff_manual.html。","$inpath/$struc_relative_dir/BMK_1_SNP_Analysis/final.*.anno.gatk.all.list.xls","$struc_relative_dir");

$writer->emptyTag('p','desc'=>"根据SNP位点碱基替换的不同方式，可以将SNP位点分为转换（Transition）和颠换（Transversion）两种类型。根据SNP位点的等位（Allele）数目，可以将SNP位点分为纯合型SNP位点（只有一个等位）和杂合型SNP位点（两个或多个等位）。不同物种杂合型SNP所占的比例存在差异。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"对各样品筛选出的SNP位点数目、转换类型比例、颠换类型比例以及杂合型SNP位点比例进行统计，如下表：",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>'注：BMK-ID：百迈客对样品的统一编号；SNP Number：SNP位点总数；Genic SNP：基因区SNP位点总数；Intergenic SNP：基因间区SNP位点总数；Transition：转换类型的SNP位点数目在总SNP位点数目中所占的百分比；Transversion：颠换类型的SNP位点数目在总SNP位点数目中所占的百分比；Heterozygosity：杂合型SNP位点数目在总SNP位点数目中所占的百分比。','type'=>"part",'name'=>"表$tid. SNP位点统计表",'path'=>"$struc_relative_dir/BMK_1_SNP_Analysis/AllSample.snp.stat.xls");

$writer->emptyTag('p','desc'=>"SNP突变类型统计分布如下图所示：",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"注：横轴为SNP突变类型，纵轴为相应的SNP数目。",'name'=>"图$pid. SNP突变类型分布图",'type'=>"type1",'path'=>"$struc_relative_dir/BMK_1_SNP_Analysis/BMK_3_SNP_type/All.snp.type.png");
$writer->emptyTag('p','desc'=>"将每个基因的SNP位点数目除以基因的长度，得到每个基因的SNP位点密度值，统计所有基因的SNP位点密度值并做密度分布图。",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"注：横轴为基因上平均每1000bp序列中分布的SNP数目，纵轴为基因数。",'name'=>"图$pid. SNP密度分布图",'type'=>"type1",'path'=>"$struc_relative_dir/BMK_1_SNP_Analysis/AllSample.SNP_density.png");
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"SNP识别结果文件、位点统计结果及突变类型统计结果路径",'action'=>"xls",'path'=>"$struc_relative_dir/BMK_1_SNP_Analysis/",'type'=>"xls");
}
$writer->emptyTag('p','desc'=>"采用SnpEff分别对SNP，InDel注释。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"常见SNP类型可分为：",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(1) INTERGENIC: intergenic_region 突变发生在基因间区；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(2) INTRAGENIC: intragenic_variant 突变发生在基因区，但是在转录本所有的属性区域外；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(3) INTRON: intron_variant 突变发生在内含子区；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(4) UPSTREAM: upstream_gene_variant 突变发生在基因上游（默认长度：5K bases）；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(5) DOWNSTREAM: downstream_gene_variant 突变发生在基因下游（默认长度：5K bases）；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(6) UTR_5_PRIME: 5_prime_UTR_variant 突变发生在5'UTR区；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(7) UTR_3_PRIME: 3_prime_UTR_variant 突变发生在3'UTR区；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(8) SPLICE_SITE_ACCEPTOR: splice_acceptor_variant 突变发生在可变剪切受体位点（一般认为外显子起始位置的前两个碱基，第一个外显子除外）；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(9) SPLICE_SITE_DONOR: splice_donor_variant 突变发生在可变剪切施体位点（一般认为编码外显子结束位置的后两个碱基，最后一个外显子除外）；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(10) SPLICE_SITE_REGION: splice_region_variant 突变发生在可变剪切区域，外显子的1-3个碱基或者内含子的3-8个碱基；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(11) START_GAINED: 5_prime_UTR_premature start_codon_gain_variant 突变发生在5‘UTR产生起始密码子的3碱基序列中；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(12) START_LOST: start_lost 突变发生在起始密码子中，并使得起始密码子突变为非起始密码子。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(13) SYNONYMOUS_CODING: synonymous_variant 同义突变；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(14) NON_SYNONYMOUS_CODING: missense_variant 错义突变；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(15) START_LOST: start_lost 突变发生在起始密码子中，并使得起始密码子突变为非起始密码子；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(16) SYNONYMOUS_STOP：stop_retained_variant 突变使得终止密码子突变为另一个终止密码子；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(17) STOP_GAINED：stop_gained 突变产生一个终止密码子；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(18) STOP_LOST：stop_lost 突变使得终止密码子突变为非终止密码子。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"SNP/Indel的注释结果统计如下所示：",'type'=>"type1");
&piclist("SNP/Indel注释分类图","注：纵轴为SNP/Indel所在区域或类型，横轴为分类数目。","$inpath/$struc_relative_dir/BMK_1_SNP_Analysis/BMK_*_anno/all.*.anno.stat.png","$struc_relative_dir");
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"SNP注释结果路径",'action'=>"xls",'path'=>"$struc_relative_dir/BMK_1_SNP_Analysis/BMK_2_SNP_anno/",'type'=>"xls");
	$writer->emptyTag('file','desc'=>"",'name'=>"Indel注释结果路径",'action'=>"xls",'path'=>"$struc_relative_dir/BMK_1_SNP_Analysis/BMK_1_InDel_anno/",'type'=>"xls");
}
	#####新可变剪接事件预测
$H3++;
$writer->emptyTag('h3','name'=>"2.$h.$H3 新可变剪接事件预测",'type'=>'type1','desc'=>"2.$h.$H3 新可变剪接事件预测");
$writer->emptyTag('p','desc'=>"基因转录生成的前体mRNA（pre-mRNA），有多种剪接方式，选择不同的外显子，产生不同的成熟mRNA，从而翻译为不同的蛋白质，构成生物性状的多样性。这种转录后的mRNA加工过程称为可变剪接或选择性剪接（Alternative Splicing）。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"用ASprofile<a href=\"#ref22\">[22]</a>软件对StringTie预测出的基因模型对每个样品的12类可变剪切事件分别进行分类和数量统计。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"12类可变剪切事件可分为：",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(1) TSS: Alternative 5' first exon (transcription start site) 第一个外显子可变剪切；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(2) TTS: Alternative 3' last exon (transcription terminal site) 最后一个外显子可变剪切；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(3) SKIP: Skipped exon(SKIP_ON,SKIP_OFF pair) 单外显子跳跃；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(4) XSKIP: Approximate SKIP (XSKIP_ON,XSKIP_OFF pair) 单外显子跳跃（模糊边界）；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(5) MSKIP: Multi-exon SKIP (MSKIP_ON,MSKIP_OFF pair) 多外显子跳跃；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(6) XMSKIP: Approximate MSKIP (XMSKIP_ON,XMSKIP_OFF pair) 多外显子跳跃（模糊边界）；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(7) IR: Intron retention (IR_ON, IR_OFF pair) 单内含子滞留；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(8) XIR: Approximate IR (XIR_ON,XIR_OFF pair) 单内含子滞留（模糊边界）；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(9) MIR: Multi-IR (MIR_ON, MIR_OFF pair) 多内含子滞留；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(10) XMIR: Approximate MIR (XMIR_ON, XMIR_OFF pair) 多内含子滞留（模糊边界）；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(11) AE: Alternative exon ends (5', 3', or both) 可变 5'或3'端剪切；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(12) XAE: Approximate AE 可变 5'或3'端剪切（模糊边界）。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"各个样品中预测的12类可变剪接事件数量统计见下图",'type'=>"type1");
$pid++;
&piclist("可变剪切类型统计图","注：图中横坐标表示属于一种可变剪接的转录本数，纵坐标表示12种可变剪切类型。","$inpath/$struc_relative_dir/BMK_4_Alt_splice/As_event_stat*.png","$struc_relative_dir");
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"可变剪切分析结果路径",'action'=>"xls",'path'=>"$struc_relative_dir/BMK_4_Alt_splice/",'type'=>"xls");
}
}
if(defined $circ_dir){
$H3++;
$writer->emptyTag('h3','name'=>"2.$h.$H3 circRNA可变剪接分析",'type'=>'type1','desc'=>"2.$h.$H3 circRNA可变剪接分析");
$writer->emptyTag('p','desc'=>"统计位于同一基因上，序列间相互部分重叠的circRNA，见下表：",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>'注：CircRNA_ID：CircRNA名称；geneLength：CircRNA长度；gene_id：CircRNA来源基因，其余为各个样品的junction reads值。','type'=>"part",'name'=>"表$tid. CircRNA可变剪>接分析",'path'=>"HTML/template/alitisplice.xls");
$writer->emptyTag('file','desc'=>"",'name'=>"CircRNA可变剪接分析文件",'action'=>"xls",'path'=>"$circ_relative_dir/BMK_1_circRNA_Prediction/overlap_alitisplice.xls",'type'=>"xls");
}

if(defined $lnc_dir){
	###基因结构优化
$H3++;
$writer->emptyTag('h3','name'=>"2.$h.$H3 基因结构优化分析",'type'=>'type1','desc'=>"2.$h.$H3 基因结构优化分析");
$writer->emptyTag('p','desc'=>"由于使用的软件或数据本身的局限性，导致所选参考基因组的注释往往不够精确，这样就有必要对原有注释的基因结构进行优化。如果在原有基因边界之外的区域有连续的Mapped Reads支持，将基因的非翻译区（Untranslated Region，UTR）向上下游延伸，修正基因的边界。此项目对$info{optimization}个基因结构进行了优化，基因结构优化结果见下表：",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>"注：GeneID：基因ID；Locus：基因座，格式为“染色体编号:起点坐标-终点坐标”；Strand：正负链；Site：优化的位置，3'或5'UTR；OriginalRegion：原来注释的第一个或最后一个外显子的起止坐标；OptimizedRegion：延伸之后的第一个或最后一个外显子的起止坐标。",'type'=>"part",'name'=>"表$tid. 基因结构优化表",'path'=>"HTML/template/geneStructure.optimize.xls");
$writer->emptyTag('file','desc'=>"",'name'=>"基因结构优化结果",'action'=>"xls",'path'=>"$struc_relative_dir/BMK_3_Gene_Structure_Optimize/$config{Project_key}.geneStructure.optimize.xls",'type'=>"xls");

	###融合基因
if(defined $config{medical}){
        $H3++;
        $writer->emptyTag('h3','name'=>"2.$h.$H3 融合基因",'type'=>'type1','desc'=>"2.$h.$H3 融合基因");
        $writer->emptyTag('p','desc'=>"融合基因是指将两个或多个基因的编码区首尾相连，置于同一套调控序列(包括启动子、增强子、核糖体结合序列、终止子等)控制之下，构成的嵌合基因。融合基因的表达产物为融合蛋白。我们使用Fusionmap<a href=\"#ref23\">[23]</a>在转录组中研究基因融合事件。Fusionmap首先通过比对到基因组和转录本中双末端(pairend)关系的序列寻找候选的基因融合，然后采用通过与nt等数据库比较，过滤掉假阳性结果。",'type'=>"type1");
        &piclist("检测到的基因融合事件","注：红色的线代表同一染色体上发生的融合事件，绿色的线代表不同染色体上发生的融合事件。","$inpath/$struc_relative_dir/BMK_5_Gene_Fusion/*_circos.png","$struc_relative_dir");
        $writer->emptyTag('p','desc'=>"融合基因事件统计",'type'=>"type1");
        &filelist("注：Fusion ID：基因融合事件编号；Unique Mapping Position：映射到唯一位点的基因融合数目；Count：基因融合位点总数；Gene：发生融合的基因；Strand：基因融合发生在正义链“+”还是反义链“-”；Chromosome：融合基因位于的染色体；Start：融合的基因的起始位点；End：融合的基因的终止位点；Filter：空。","$inpath/$struc_relative_dir/BMK_5_Gene_Fusion/*FusionReport.xls","$struc_relative_dir");
	if(!defined $coud){
		$writer->emptyTag('file','desc'=>"",'name'=>"融合基因分析结果路径",'action'=>"xls",'path'=>"$struc_relative_dir/BMK_5_Gene_Fusion/",'type'=>"xls");
	}
}
}

if(defined $config{medical} ||defined $config{phastCons}){
	$H3++;
	if(defined $lnc_dir && defined $circ_dir){
		$writer->emptyTag('h3','name'=>"2.$h.$H3 保守性分析",'type'=>'type1','desc'=>"2.$h.$H3 保守性分析");
		$writer->emptyTag('p','desc'=>"LncRNA 的序列保守性相对 mRNA 要低，使用 phastCons<a href=\"#ref24\">[24]</a> 软件 (http://compgen.bscb.cornell.edu/phast/) 分别对 mRNA 和 lncRNA 以及circRNA进行保守性打分，得到 lncRNA 和 mRNA 以及circRNA 的保守性分值累积分布图如下所示：",'type'=>"type1");
		$writer->emptyTag('p','desc'=>"不同类型RNA的保守性比较如下：",'type'=>"type1");
		&piclist("不同类型RNA保守性比较","注：上述三张图分别是保守性分值柱状图，保守性分值累积分布图,保守性分值概率密度曲线。","$inpath/$struc_relative_dir/BMK_6_conservation/All*.png","$struc_relative_dir");
		if(!defined $cloud){
			$writer->emptyTag('file','desc'=>"",'name'=>"保守性分析结果路径",'action'=>"xls",'path'=>"$struc_relative_dir/BMK_6_conservation/",'type'=>"xls");
		}
	}elsif(defined $lnc_dir && !defined $circ_dir){
		$writer->emptyTag('h3','name'=>"2.$h.$H3 保守性分析",'type'=>'type1','desc'=>"2.$h.$H3 保守性分析");
		$writer->emptyTag('p','desc'=>"LncRNA 的序列保守性相对 mRNA 要低，使用 phastCons<a href=\"#ref24\">[24]</a> 软件 (http://compgen.bscb.cornell.edu/phast/) 分别对 mRNA 和 lncRNA 进行保守性打分，得到 lncRNA 和 mRNA 的保守性分值累积分布图如下所示：",'type'=>"type1");
		$writer->emptyTag('p','desc'=>"不同类型RNA的保守性比较如下：",'type'=>"type1");
		&piclist("不同类型RNA保守性比较","注：上述三张图分别是保守性分值柱状图，保守性分值累积分布图,保守性分值概率密度曲线。","$inpath/$struc_relative_dir/BMK_6_conservation/All*.png","$struc_relative_dir");
		if(!defined $cloud){
			$writer->emptyTag('file','desc'=>"",'name'=>"保守性分析结果路径",'action'=>"xls",'path'=>"$struc_relative_dir/BMK_6_conservation/",'type'=>"xls");
		}
	}elsif(!defined $lnc_dir && defined $circ_dir){
                $writer->emptyTag('h3','name'=>"2.$h.$H3 circRNA保守性分析",'type'=>'type1','desc'=>"2.$h.$H3 circRNA保守性分析");
                $writer->emptyTag('p','desc'=>"使用 phastCons<a href=\"#ref24\">[24]</a> 软件 (http://compgen.bscb.cornell.edu/phast/) 对circRNA 进行保守性打分，得到circRNA 的保守性分值累积分布图如下所示：
",'type'=>"type1");
                &piclist("不同类型RNA保守性比较","注：上述三张图分别是保守性分值柱状图，保守性分值累积分布图,保守性分值概率密度曲线。","$inpath/$circ_relative_dir/BMK_1_circRNA_Prediction/BMK_4_conservation/All*.png","$circ_relative_dir");
		if(!defined $cloud){
	                $writer->emptyTag('file','desc'=>"",'name'=>"保守性分析结果路径",'action'=>"xls",'path'=>"$circ_relative_dir/BMK_1_circRNA_Prediction/BMK_4_conservation/",'type'=>"xls");
		}
        }
}

############################
$writer->emptyTag('h1','desc'=>'3. 结果文件查看说明','type'=>"type1",'name'=>'3. 结果文件查看说明');
$writer->emptyTag('p','desc'=>"(1) 上传目录中有Readme.txt说明，详细介绍了每个文件所代表的内容。上传的结果数据文件多以文本格式为主(fa文件、txt文件，xls文件等)。在Windows系统下查看文件，推荐使用Editplus或UltraEdit或Notepad作为文本浏览程序，否则会因文件过大造成死机。如果文件大小超过50M建议用pilotedit打开；在Unix或Linux系统下可以浏览较大的文本文件，用less等操作命令可以顺利地查看。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"(2) 阅读报告时请注意：报告中出现的路径链接，只有所有结果文件都在的时候才可以查看，只有单独的结题报告而没有结果文件时是无法查看的。",'type'=>"type1");
$writer->startTag('ref_list',desc=>"参考文献",type=>"type1",name=>"参考文献");
&reference ($writer,[
	["Kim D, Langmead B, Salzberg S L. HISAT: a fast spliced aligner with low memory requirements.[J]. Nature Methods, 2015, 12(4):357-360.","https://www.nature.com/articles/nmeth.3317","1"],
	["Pertea M, Kim D, Pertea G, et al. Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie, and Ballgown[J]. Nature Protocols, 2016, 11(9):1650.","https://www.nature.com/articles/nprot.2016.095","2"],
	["Kong L, Zhang Y, Ye Z Q, et al. CPC: assess the protein-coding potential of transcripts using sequence features and support vector machine[J]. Nucleic acids research, 2007, 35(suppl 2): W345-W349.","https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1933232/?tool=pubmed","3"],
	["Sun L, Luo H, Bu D, et al. Utilizing sequence intrinsic composition to classify protein-coding and long non-coding transcripts[J]. Nucleic acids research, 2013: gkt646.","https://academic.oup.com/nar/article/41/17/e166/2411728","4"],
	["Wang L, Park H J, Dasari S, et al. CPAT: Coding-Potential Assessment Tool using an alignment-free logistic regression model[J]. [Nucleic acids research]. 2013, 41(6): e74-e74.","https://academic.oup.com/nar/article/41/6/e74/2902455","5"],
	["Finn R D, Bateman A, Clements J, et al. Pfam: the protein families database[J]. Nucleic acids research, 2013: gkt1223.","https://academic.oup.com/nar/article/42/D1/D222/1062431","6"],
	["Chen G, Wang Z, Wang D, et al. LncRNADisease: a database for long-non-coding RNA-associated diseases[J]. Nucleic Acids Research, 2013, 41(Database issue):983-6.","https://www.ncbi.nlm.nih.gov/pubmed/23175614","7"],
	["Kozomara A, Griffithsjones S. miRBase: annotating high confidence microRNAs using deep sequencing data[J]. Nucleic Acids Research, 2014, 42(Database issue):D68.","https://academic.oup.com/nar/article/42/D1/D68/1057911","8"],
	["Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM[J]. arXiv preprint arXiv:1303.3997, 2013.","https://arxiv.org/pdf/1303.3997.pdf","9"],
	["Gao Y, Wang J, Zhao F. CIRI: an efficient and unbiased algorithm for de novo circular RNA identification[J]. Genome biology, 2015, 16(1): 1.","https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0571-3","10"],
	["Memczak S, Jens M, Elefsinioti A, et al. Circular RNAs are a large class of animal RNAs with regulatory potency.[J]. Nature, 2014, 495(7441):333-338.","https://www.nature.com/articles/nature11928","11"],
	["Glažar P, Papavasileiou P, Rajewsky N. circBase: a database for circular RNAs. RNA. 2014 Sep 18. ","http://rnajournal.cshlp.org/content/early/2014/09/18/rna.043687.113.abstract","12"],
	["Langmead, B., Trapnell, C., Pop, M., and Salzberg, S.L. (2009). Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome biology 10,R25.","http://bowtie-bio.sourceforge.net/index.shtml","13"],
	["Friedlander,M.R., Mackowiak,S.D., Li,N., Chen,W. and Rajewsky,N. (2012) miRDeep2 accurately identifies known and hundreds of novel microRNA genes in seven animal clades. Nucleic Acids Res., 40, 37-52.","https://academic.oup.com/nar/article/40/1/37/1275937","14"],
	["Bo X, Wang S. TargetFinder: a software for antisense oligonucleotide target site selection based on MAST and secondary structures of target mRNA[J]. Bioinformatics, 2005, 21(8):1401.","https://academic.oup.com/bioinformatics/article/21/8/1401/250353","15"],
	["Betel D, Wilson M, Gabow A, et al. The microRNA. org resource: targets and expression[J]. Nucleic acids research, 2008, 36(suppl 1): D149-D153.","https://academic.oup.com/nar/article/36/suppl_1/D149/2508385","16"],
	["Lewis B P, Shih I H, Jonesrhoades M W, et al. Prediction of mammalian microRNA targets.[J]. Cell, 2003, 115(7):787-798.","https://www.cell.com/cell/pdf/S0092-8674(03)01018-3.pdf","17"],
	["Dr. Qinghua Cui, 38 Xueyuan Rd, Department of Biomedical Informatics, Peking University Health Science Center, Beijing 100191, China","http://210.73.221.6/hmdd#fragment-4","18"],
	["McKenna A, Hanna M, Banks E, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data[J]. Genome research, 2010, 20(9): 1297-1303.","https://genome.cshlp.org/content/20/9/1297.short","19"],
	["Cingolani P, Platts A, Wang L L, et al. A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3[J]. Fly, 2012, 6(2): 80-92.","https://www.tandfonline.com/doi/abs/10.4161/fly.19695#.Vk6TaPRAUwo","20"],
	["Florea L, Song L, Salzberg S L. Thousands of exon skipping events differentiate among splicing patterns in sixteen human tissues [v1; ref status: indexed[J]. 2013.","https://f1000research.com/articles/2-188/v2","21"],
	["Ge H, Liu K, Juan T, Fang F, Newman M, Hoeck W (2011) FusionMap: detecting fusion genes from next-generation sequencing data at base-pair resolution. Bioinformatics","https://academic.oup.com/bioinformatics/article/27/14/1922/194689","22"],
	["de Oliveira L F, Christoff A P, Margis R. isomiRID: a framework to identify microRNA isoforms.[J]. Bioinformatics, 2013, 29(20):2521-3.","https://academic.oup.com/bioinformatics/article/29/20/2521/276800","23"],
	["Hubisz M J, Pollard K S, Siepel A. PHAST and RPHAST: phylogenetic analysis with space/time models[J]. Briefings in Bioinformatics, 2011, 12(1):41-51.","https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3030812/","24"],
]);
$writer->endTag('ref_list');

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

sub get_data_num{
        my @B;
        my $cfg = shift;
        open (IN,"<$cfg");
        while(<IN>){
                if($_=~/^Sample\s+/){
                        my @A = split /\s+/,$_;
                        push @B,$A[1];
                }
        }
        close IN;
        return scalar(@B);
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
		$_=~s/\s+$//;
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

Example: perl $0 -i Web_Report -o Web_Report/configtest.xml -cfg1 data.cfg -cfg2 detail.cfg

USAGE
        print $usage;
        exit;
}

