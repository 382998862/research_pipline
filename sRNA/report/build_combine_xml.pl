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
my ($detail_cfg,$combine,$data_cfg,$type,$output,$cloud);
GetOptions(
        "h|?"           =>\&USAGE,
	"data_cfg:s"	=>\$data_cfg,
	"combine:s"	=>\$combine,
	"type:s"	=>\$type,
	"detail_cfg:s"	=>\$detail_cfg,
        "o:s"           =>\$output,
	"cloud"		=>\$cloud,
)or &USAGE;
&USAGE unless ($combine and $detail_cfg and $data_cfg and $output and $type);

$data_cfg= abs_path($data_cfg);
$combine=abs_path($combine);
$cfg=abs_path($detail_cfg);
$type||="combine";

my $template_dir="template";
my $od=dirname $output;
my %info=&readConfig("$od/$template_dir/stat.info");
my $relate_dir = (split /\//,$combine)[-2];  #BMK_Result/Web_Report
my $combine_dir = (split /\//,$combine)[-1]; #BMK_7_Combine
print "$combine_dir\n";

my %config=&readcfg($detail_cfg);
`mkdir -p $od/$type`	unless(-d "$od/$type");

my $relate_cloud;
if(defined $cloud){
	$relate_cloud = "../";
}else{
	$relate_cloud = "../../";
}

my ($pid,$tid)=(0,0);

my $report_time=&gaintime();
my $writer = XML::Writer->new(OUTPUT => 'self');
$writer->xmlDecl('UTF-8');
$writer->startTag('report');
$writer->emptyTag('report_version','value'=>$config{report_version});
my $report_name = "联合分析";
$writer->emptyTag('report_name','value'=>$report_name);
$writer->emptyTag('report_code','value'=>$config{report_code});
$writer->emptyTag('report_user','value'=>$config{report_user});
$writer->emptyTag('report_user_addr','value'=>$config{report_user_addr});
$writer->emptyTag('report_time','value'=>"XXXX/XX/XX;XXXX/XX/XX;XXXX/XX/XX;$report_time");
$writer->emptyTag('report_abstract','value'=>'');
###
$h1=1;
$writer->emptyTag('h1','name'=>"$h1 样品表达量全局展示",'type'=>'type1','desc'=>"$h1 样品表达量全局展示");
$writer->emptyTag('p','desc'=>"利用高通量测序得到的FASTQ数据，分别对几种RNA进行预测和鉴定，结合已知的RNA，统计几种RNA的数量。结果如下图所示：",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. 不同RNA数量示意图 ",'type'=>"type1",'path'=>"../$combine_dir/BMK_1_Global/BMK_1_All/Total.pie.png");
$writer->emptyTag('p','desc'=>"对每个样本，分别提取mRNA、lncRNA、circRNA、miRNA的表达水平，根据RNA的位置信息，用circos图展示样本中RNA的表达水平。这样，我们可以全局观察每个样本多组学之间的关联。",'type'=>"type1");
&piclist("表达量全局展示","注：最外圈为染色体信息，之后依次是mRNA、lncRNA、circRNA、miRNA的表达水平。这里，表达量均取对数，miRNA位置采用其前体RNA的位置信息。","$combine/BMK_1_Global/BMK_1_All/*/*circos.png","$relate_dir");
$writer->emptyTag('p','desc'=>"通常情况下越长的染色体拥有的基因数目越多，因此展示不同类型的RNA在不同染色体表达量的分布情况，如下所示：",'type'=>"type1");
&piclist("RNA在不同染色体的表达分布","注：不同颜色代表不同的RNA类型，横坐标代表不同染色体，纵坐标代表RNA在该染色体上表达量的比例。同一种RNA在不同染色体上的总分布和为1。这里，表达量均取对数。","$combine/BMK_1_Global/BMK_1_All/*/*exp_barplot.png","$relate_dir");
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"样品表达量全局展示结果路径",'action'=>"xls",'path'=>"../../$combine_dir/BMK_1_Global/BMK_1_All/",'type'=>"xls");
}

$h1++;
$writer->emptyTag('h1','name'=>"$h1 样品差异表达全局展示",'type'=>'type1','desc'=>"$h1 样品差异表达全局展示");
$writer->emptyTag('p','desc'=>"全基因组范围内的差异分析可以发现多组样品的差异mRNA、circRNA、miRNA、lncRNA。对于每种RNA类型的数据，对差异分组之间进行差异RNA筛选，同一差异组合不同类型RNA所筛选出的差异RNA数目
如下图所示：",'type'=>"type1");
&piclist("差异表达的RNA数目分布","","$combine/BMK_1_Global/BMK_2_DEG/*/*.pie.png","$relate_dir");
$writer->emptyTag('p','desc'=>"全基因组范围内观察多组学数据的差异情况有助于发现差异调控机制，差异RNA显著性circos图如下：",'type'=>"type1");
my $flag="FDR";
$flag="pvalue" if(exists $config{PValue} && !exists $config{FDR});
&piclist("差异RNA显著性circos图","注：最外圈为染色体信息，之后依次是mRNA、lncRNA、circRNA、miRNA。每个差异分组的圈图中红色代表上调，蓝色代表下调，高度代表显著性（-log10($flag)）","$combine/BMK_1_Global/BMK_2_DEG/*_vs_*/*.Circos.png","$relate_dir");
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"样品差异表达全局展示结果路径",'action'=>"xls",'path'=>"../../$combine_dir/BMK_1_Global/BMK_2_DEG/",'type'=>"xls");
}

my $ref=0;
my $sample_num=&get_data_num("$data_cfg");
if($sample_num>=5){
	$h1++;
	$ref++;
        $writer->emptyTag('h1','name'=>"$h1 共表达分析",'type'=>'type1','desc'=>"$h1 共表达分析");
	$writer->emptyTag('p','desc'=>"具有相同表达模式的基因更趋向于拥有相同的功能。",'type'=>"type1");
        $writer->emptyTag('p','desc'=>"样品数量不少于5时，采用皮尔森相关性分析，分别对差异RNA之间构建mRNA-lncRNA, miRNA-circRNA共表达网络<a href=\"#ref1\">[$ref]</a>，筛选条件是共表达的RNA对之间相关性的绝对值大于$config{coexp_cor}，显著性p值小于$config{coexp_p}。样品数目小于5时，不进行此项分析，此部分分析内容为空。",'type'=>"type1");
        $writer->emptyTag('p','desc'=>"差异共表达关系见下表：",'type'=>"type1");
        &filelist("","$combine/BMK_5_coexp/Diff.Sig.*xls","$combine_dir");
	if(!defined $cloud){
		$writer->emptyTag('file','desc'=>"",'name'=>"共表达分析结果路径",'action'=>"xls",'path'=>"../../$combine_dir/BMK_5_coexp/",'type'=>"xls");
	}
}

$h1++;
$ref++;
$writer->emptyTag('h1','name'=>"$h1 竞争性内源RNA网络",'type'=>'type1','desc'=>"$h1 竞争性内源RNA网络");
$writer->emptyTag('p','desc'=>"竞争性内源RNA(ceRNA)<a href=\"#ref1\">[$ref]</a>是近几年来的研究热点，是一种全新的转录调控模式。这些ceRNA能够通过miRNA应答元件（microRNA Response Element, MRE）竞争性的结合相同的miRNA来调节彼此的表达水平。ceRNA假说揭示了一种RNA间相互作用的新机制。ceRNA可以通过竞争性的结合miRNA来调控彼此的表达，如：microRNA可以导致基因沉默，而lncRNA竞争性的结合了microRNA，从而影响了microRNA导致基因沉默这一功能。这里，mRNA、lncRNA、circRNA等均可作为ceRNA。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"ceRNA通常具有以下特点：",'type'=>"type1");
$writer->emptyTag('p','desc'=>"a）ceRNA均受miRNA调控；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"b）ceRNA之间需具有相同的miRNA结合位点；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"c）ceRNA之间存在相互调控关系，且表达水平变化趋势一致。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"本项目通过miRNA的靶向关系来获取候选的ceRNA关系对，满足如下条件：",'type'=>"type1");
$writer->emptyTag('p','desc'=>"1) ceRNA之间相同的miRNA数目应大于$config{ceRNA_common}；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"2) 超几何检验p值小于 $config{ceRNA_pvalue}，校正后的fdr值小于 $config{ceRNA_fdr}；",'type'=>"type1");
$writer->emptyTag('p','desc'=>"3) 当样品数量大于等于5时，可以考虑ceRNA的表达水平，将共表达分析的结果考虑进来，即得到ceRNA之间两两均共表达的ceRNA共表达网络。否则不做考虑。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"样品数量不少于5时，此处展示的是考虑共表达后的ceRNA结果。否则，展示的只是ceRNA的结果。",'type'=>"type1");

if($sample_num>=5){
        $writer->emptyTag('p','desc'=>"ceRNA网络（考虑共表达结果）中共包含 $info{ceRNA_relate} 条边，$info{ceRNA_node} 个点，其中lncRNA $info{ceRNA_lncRNA} 个，mRNA $info{ceRNA_gene} 个，circRNA $info{ceRNA_circRNA} 个。",'type'=>"type1");
}else{
        $writer->emptyTag('p','desc'=>"ceRNA网络（不考虑共表达结果）中共包含 $info{ceRNA_relate} 条边，$info{ceRNA_node} 个点，其中lncRNA $info{ceRNA_lncRNA} 个，mRNA $info{ceRNA_gene} 个，circRNA $info{ceRNA_circRNA} 个。",'type'=>"type1");
}

$tid++;
$writer->emptyTag('table','desc'=>'注：#ceRNA1为ceRNA1的ID；ceRNA2为ceRNA2的ID；ceRNA1_num为ceRNA1所结合的miRNA数目；ceRNA2_num为ceRNA2所结合的miRNA数目；overlap_num为ceRNA1和ceRNA2共同结合的miRNA数目；ceRNA1_miRNA_uniq为ceRNA1所结合的特有miRNA；ceRNA2_miRNA_uniq为ceRNA2所结合的特有miRNA；overlap_miRNA为ceRNA1和ceRNA2共同结合的miRNA；pvalue为超几何检验p值；FDR为错误发现率FDR值；coexp_cor为ceRNA关系对之间的皮尔森相关系数；coexp_pvalue为ceRNA关系对之间的皮尔森相关性的显著性p值（若样品数目小于5则不存在最后两列信息）。这里，ceRNA1、ceRNA2可由lncRNA、circRNA、mRNA充当。','type'=>"part",'name'=>"表$tid. ceRNA关系表",'path'=>"$template_dir/ceRNA.xls");
my $ceRNA_file = (split /$relate_dir\//,(glob "$combine/BMK_3_ceRNA_random/*ceRNA_pair_adjust_p_Sig.xls"))[1];
$writer->emptyTag('file','desc'=>"",'name'=>"ceRNA关系表",'action'=>"xls",'path'=>"../../$ceRNA_file",'type'=>"xls");
$writer->emptyTag('p','desc'=>"从ceRNA关系对中对每个差异组合提取差异RNA的一步近邻网络，得到的差异ceRNA关系示例表如下：",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>'注：#ceRNA1为ceRNA1的ID；ceRNA1_regulated为（差异）ceRNA1的上下调，若不是差异则用"--"代替；ceRNA2为ceRNA2的ID；ceRNA2_regulated为（差异）ceRNA2的上下调，若不是差异则用"--"代替；ceRNA1_num为ceRNA1所结合的miRNA数目；ceRNA2_num2为ceRNA2所结合的miRNA数目；overlap_num为ceRNA1和ceRNA2共同结合的miRNA数目；ceRNA1_miRNA_uniq为ceRNA1所结合的特有miRNA；ceRNA2_miRNA_uniq为ceRNA2所结合的特有miRNA；overlap_miRNA为ceRNA1和ceRNA2共同结合的miRNA；pvalue为超几何检验p值；FDR为错误发现率FDR值；coexp_cor为ceRNA关系对之间的皮尔森相关系数；coexp_pvalue为ceRNA关系对之间的皮尔森相关性的显著性p值（若样品数目小于5则不存在最后两列信息）。这里，ceRNA1、ceRNA2可由lncRNA、circRNA、mRNA充当。','type'=>"part",'name'=>"表$tid. 差异ceRNA关系表（示例）",'path'=>"$template_dir/sub_network.xls");
&filelist("","$combine/*/*/sub_*_network.xls","$combine_dir");

if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"ceRNA分析结果路径",'action'=>"xls",'path'=>"../../$combine_dir/BMK_3_ceRNA_random",'type'=>"xls");
}

$h1++;
$writer->emptyTag('h1','name'=>"$h1 关键基因通路整合分析",'type'=>'type1','desc'=>"$h1 关键基因通路整合分析");
$writer->emptyTag('p','desc'=>"对于这种复杂的网络，提取网络中的关键信息是非常必要的。基于随机游走可以对网络中的节点的关键性进行排序，获取网络中的关键RNA。采用随机游走中的经典算法PageRank来获取网络中所有节点的打分（即重要性）。选择网络中打分排布在前$config{keyratio} 的点，即关键RNA，作为重点的研究对象。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"对关键RNA中的基因进行通路富集分析，选择最显著富集的前5个通路。提取这5个通路中基因与基因之间的关系整合成通路网络，并将关键基因map到通路网络如下：",'type'=>"type1");
my @network_png = glob "$combine/BMK_3_ceRNA_random/*/*network.png";
if(scalar(@network_png)>=1){
        &piclist("KEGG整合通路网络","注：每个圆点代表一个基因，每个长方形代表一个通路，线代表通路中基因和基因或其他通路之间的关系。不同线的颜色代表关系来自不同的通路。红色的点为关键基因。本报告中所有网络图
均采用R语言的iGraph包通过脚本进行绘制，同时也会提供绘制网络图所需文件，研究者可采用Cytoscape绘制自己所感兴趣的样式。","$combine/BMK_3_ceRNA_random/*/*network.png","$relate_dir");
}else{
        $pid++;
        $writer->emptyTag('pic','desc'=>"注：每个圆点代表一个基因，每个长方形代表一个通路，线代表通路中基因和基因或其他通路之间的关系。不同线的颜色代表关系来自不同的通路。红色的点为关键基因。本报告中所有网
络图均采用R语言的iGraph包通过脚本进行绘制，同时也会提供绘制网络图所需文件，研究者可采用Cytoscape绘制自己所感兴趣的样式。",'name'=>"图$pid. KEGG整合通路网络(示例)",'type'=>"type1",'path'=>"$template_dir/demo.network.png");
}
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"关键基因通路整合分析结果路径",'action'=>"xls",'path'=>"../../$combine_dir/BMK_3_ceRNA_random/",'type'=>"xls");
}

$h1++;
$writer->emptyTag('h1','name'=>"$h1 差异RNA靶向关系联合分析",'type'=>'type1','desc'=>"$h1 差异RNA靶向关系联合分析");
$writer->emptyTag('p','desc'=>"阅读前必看！！！",'type'=>'type1');
$writer->emptyTag('p','desc'=>"此部分分析是基于差异RNA的结果做的数据整理以及统计，所以当差异结果较少或者靶向关系较少时，部分结果可能为空。",'type'=>'type1');
$writer->emptyTag('p','desc'=>"6.1-6.4，展示表格均为本项目的实际结果，若有多个差异组合，仅以其中一个差异组合为例展示我们的结果，表格中的\"--\"表示没有靶向关系。展示表格对应的结果文件为*_vs_*.DEG.xls。基于结果对应整理了靶向关系的Cytoscape输入文件：*Interaction.list和*Attribution.list文件分别是Cytoscape输入的互作文件和对应的属性文件，若这些文件为空，表示对应的靶向结果为空。具体结果可到结果文件的BMK_7_Combine/BMK_4_cytoscape/*RNA/中查看，如果结果为空，则没有对应的图片。",'type'=>'type1');
$writer->emptyTag('p','desc'=>"6.5-6.6，若结果不为空，展示表格是本项目的实际结果；若结果为空，展示表格为示例表格，对应的Cytoscape输入文件亦为空。Cytoscape展示图为示例图。具体结果可到结果文件的$combine_dir/BMK_4_cytoscape/*RNA-miRNA-mRNA/",'type'=>'type1');
$writer->emptyTag('h2','name'=>"$h1.1 差异LncRNA相关靶向RNA分析",'type'=>'type1','desc'=>"$h1.1 差异LncRNA相关靶向RNA分析");
$writer->emptyTag('p','desc'=>"以差异LncRNA为中心进行靶向关系联合分析，将差异LncRNA及与其具有靶向关系的差异gene和差异miRNA，整理成表格。结果如下表所示：",'type'=>"type1");
$tid++;
if(-e "$combine/../HTML/$template_dir/diff_lncRNA.xls"){
	$writer->emptyTag('table','desc'=>"注：#ID：差异LncRNA的ID；*_Count：count数；*_FPKM：FPKM值；Pvalue：Pvalue值；log2FC：log2FC值；regulated：上下调（up/down）；diff_Cis_target_gene：与LncRNA有顺式靶向关系的差异gene ID；diff_Trans_target_gene（当样品数目不少于5时）：与LncRNA有反式靶向关系的差异gene ID；diff_target_miRNA：与LncRNA有靶向关系的差异miRNA ID。",'type'=>"part",'name'=>"表$tid. 差异LncRNA靶向关系结果",'path'=>"$template_dir/diff_lncRNA.xls");
}else{
	$writer->emptyTag('table','desc'=>"注：示例表格。#ID：差异LncRNA的ID；L*_Count：count数；L*_FPKM：FPKM值；Pvalue：Pvalue值；log2FC：log2FC值；regulated：上下调（up/down）；diff_Cis_target_gene：与LncRNA有顺式靶向关系的差异gene ID；diff_Trans_target_gene（当样品数目不少于5时）：与LncRNA有反式靶向关系的差异gene ID；diff_target_miRNA：与LncRNA有靶向关系的差异miRNA ID。",'type'=>"part",'name'=>"表$tid. 差异LncRNA靶向关系结果（示例）",'path'=>"$template_dir/demo/demo.diff_lncRNA.xls");
}
$writer->emptyTag('p','desc'=>"各组合差异LncRNA靶向关系结果：",'type'=>"type1");
&filelist("","$combine/BMK_4_cytoscape/lncRNA/*_vs_*.DEG.xls","$combine_dir");
my @lnc_png = glob "$combine/BMK_4_cytoscape/lncRNA/*venn.png";
if(@lnc_png==0){
	if($sample_num>=5){
		$writer->emptyTag('p','desc'=>"差异LncRNA与差异基因所靶向的所有LncRNA（当样品数目不少于5时，分顺式和反式分别展示），差异miRNA所靶向的所有LncRNA的交集以Venn图的形式展示。示例图如下：",'type'=>"type1");
		&piclist("差异LncRNA与差异基因所靶向的所有LncRNA、差异miRNA所靶向的所有LncRNA的交集（示例）","注：示例图。DE_LncRNA、DE_Cis.mRNA_TargetLncRNA、DE_Trans.mRNA_TargetLncRNA（样品数目不少于5）、DE_miRNA_TargetLncRNA分别指所有差异LncRNA、差异基因所靶向的顺式/反式LncRNA、差异miRNA所靶向的所有LncRNA","$combine/../HTML/$template_dir/demo.lncRNA.*.venn.png","$relate_dir");
	}else{
		$pid++;
		$writer->emptyTag('p','desc'=>"差异LncRNA与差异基因所靶向的所有LncRNA（当样品数目少于5时，只有顺式的结果），差异miRNA所靶向的所有LncRNA的交集以Venn图的形式展示。示例图如下:",'type'=>"type1");
		$writer->emptyTag('pic','desc'=>"注：示例图。DE_LncRNA、DE_Cis.mRNA_TargetLncRNA、DE_miRNA_TargetLncRNA分别指差异LncRNA、差异基因所靶向的顺式LncRNA、差异miRNA所靶向的所有LncRNA",'name'=>"图$pid. 差异LncRNA与差异基因所靶向的所有LncRNA，差异miRNA所靶向的所有LncRNA的交集（示例）",'type'=>"type1",'path'=>"HTML/$template_dir/demo.lncRNA.Cis.venn.png");
	}
}else{
	$writer->emptyTag('p','desc'=>"差异LncRNA与差异基因所靶向的所有LncRNA（当样品数目不少于5时，分顺式和反式分别展示），差异miRNA所靶向的所有LncRNA的交集以Venn图的形式展示。如下:",'type'=>"type1");
	&piclist("差异LncRNA与差异基因所靶向的所有LncRNA、差异miRNA所靶向的所有LncRNA的交集","注：DE_LncRNA、DE_Cis.mRNA_TargetLncRNA、DE_Trans.mRNA_TargetLncRNA（样品数目不少于5）、DE_miRNA_TargetLncRNA分别指差异LncRNA、差异基因所靶向的顺式/反式LncRNA、差异miRNA所靶向的所有LncRNA","$combine/BMK_4_cytoscape/lncRNA/*venn.png","$relate_dir");
}
if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"差异LncRNA相关靶向RNA分析结果路径",'action'=>"xls",'path'=>"../../$combine_dir/BMK_4_cytoscape/lncRNA",'type'=>"xls");
}

$writer->emptyTag('h2','name'=>"$h1.2 差异circRNA的Host gene及相关靶向RNA分析",'type'=>'type1','desc'=>"$h1.2 差异circRNA的Host gene及相关靶向RNA分析");
$writer->emptyTag('p','desc'=>"以差异circRNA为中心进行靶向关系联合分析，将差异circRNA及与其Host gene和与其具有靶向关系的差异miRNA，整理成表格。结果如下表所示：",'type'=>"type1");
$tid++;
if (-e "$combine/../HTML/$template_dir/diff_circRNA.xls"){
	$writer->emptyTag('table','desc'=>"注：#ID：差异circRNA的ID；*_Count：count数；*_SRPBM/TPM：表达量值；Pvalue：Pvalue值；log2FC：log2FC值；regulated：上下调（up/down）；diff_host_gene：该circRNA对应Host gene ID；diff_target_miRNA：与circRNA有靶向关系的差异miRNA ID",'name'=>"表$tid. 差异circRNA的Host gene及靶向miRNA关系结果",'type'=>"part",'path'=>"$template_dir/diff_circRNA.xls");
}else{
	$writer->emptyTag('table','desc'=>"注：示例表格。#ID：差异circRNA的ID；L*_Count：count数；L*_SRPBM/TPM：表达量值；Pvalue：Pvalue值；log2FC：log2FC值；regulated：上下调（up/down）；diff_host_gene：该circRNA对应Host gene ID；diff_target_miRNA：与circRNA有靶向关系的差异miRNA ID",'name'=>"表$tid. 差异circRNA的Host gene及靶向miRNA关系结果（示例）",'type'=>"part",'path'=>"$template_dir/demo/demo.diff_circRNA.xls");
}
$writer->emptyTag('p','desc'=>"各组合差异差异circRNA的Host gene及靶向miRNA关系结果：",'type'=>"type1");
&filelist("","$combine/BMK_4_cytoscape/circRNA/*_vs_*.DEG.xls","$combine_dir");
my @circ_png = glob "$combine/BMK_4_cytoscape/circRNA/*venn.png";
if(@circ_png==0){
	$pid++;
        $writer->emptyTag('p','desc'=>"差异circRNA，以差异基因作为Host gene的所有circRNA，差异miRNA所靶向的所有circRNA的交集以Venn图的形式展示。示例图如下:",'type'=>"type1");
	$writer->emptyTag('pic','desc'=>"注：示例图。DE_circRNA、DE_Hostgene_circRNA、DE_miRNA_TargetcircRNA分别指差异circRNA、以差异基因作为Host gene的所有circRNA、差异miRNA所靶向的所有circRNA。",'name'=>"图 $pid. 差异circRNA与以差异基因作为Host gene的所有circRNA、差异miRNA所靶向的所有circRNA的交集（示例）",'type'=>"type1",'path'=>"$template_dir/demo.circ.venn.png");
}else{
	$writer->emptyTag('p','desc'=>"差异circRNA，以差异基因作为Host gene的所有circRNA，差异miRNA所靶向的所有circRNA的交集以Venn图的形式展示。如下:",'type'=>"type1");
	&piclist("差异circRNA与以差异基因作为Host gene的所有circRNA、差异miRNA所靶向的所有circRNA的交集","注:DE_circRNA、DE_Hostgene_circRNA、DE_miRNA_TargetcircRNA分别指差异circRNA、以差异基因作为Host gene的所有circRNA、差异miRNA所靶向的所有circRNA","$combine/BMK_4_cytoscape/circRNA/*venn.png","$relate_dir");
}

if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"差异circRNA相关靶向RNA分析结果路径",'action'=>"xls",'path'=>"../../$combine_dir/BMK_4_cytoscape/circRNA",'type'=>"xls");
}

$writer->emptyTag('h2','name'=>"$h1.3 差异miRNA相关靶向RNA分析",'type'=>'type1','desc'=>"$h1.3 差异miRNA相关靶向RNA分析");
$writer->emptyTag('p','desc'=>"以差异miRNA为中心进行靶向关系联合分析，将差异miRNA及与其具有靶向关系的差异gene、差异LncRNA、差异circRNA，整理成表格。结果如下表所示：",'type'=>"type1");
$tid++;
if (-e "$combine/../HTML/$template_dir/diff_miRNA.xls"){
	$writer->emptyTag('table','desc'=>"注：#ID：差异miRNA的ID；*_Count：count数；*_TPM：FPKM值；Pvalue：Pvalue值；log2FC：log2FC值；regulated：上下调（up/down）；diff_target_gene：与miRNA有靶向关系的差异gene ID；diff_target_lncRNA：与miRNA有靶向关系的差异lncRNA ID；diff_target_circRNA：与miRNA有靶向关系的差异circRNA ID。",'name'=>"表$tid. 差异miRNA靶向关系结果",'type'=>"part",'path'=>"$template_dir/diff_miRNA.xls");
}else{
	$writer->emptyTag('table','desc'=>"注：示例表格。#ID：差异miRNA的ID；L*_Count：count数；L*_TPM：FPKM值；Pvalue：Pvalue值；log2FC：log2FC值；regulated：上下调（up/down）；diff_target_gene：与miRNA有靶向关系的差异gene ID；diff_target_lncRNA：与miRNA有靶向关系的差异lncRNA ID；diff_target_circRNA：与miRNA有靶向关系的差异circRNA ID。",'name'=>"表$tid. 差异miRNA靶向关系结果（示例）",'type'=>"part",'path'=>"$template_dir/demo/demo.diff_miRNA.xls");
}
$writer->emptyTag('p','desc'=>"各组合差异miRNA靶向关系结果：",'type'=>"type1");
&filelist("","$combine/BMK_4_cytoscape/sRNA/*_vs_*.DEG.xls","$combine_dir");

my @sRNA_png = glob "$combine/BMK_4_cytoscape/sRNA/*venn.png";
if(@sRNA_png==0){
	$pid++;
	$writer->emptyTag('p','desc'=>"差异miRNA，差异基因所靶向的所有miRNA，差异lncRNA所靶向的所有miRNA，差异circRNA所靶向的所有miRNA的交集以Venn图的形式展示。示例图如下：",'type'=>"type1");
	$writer->emptyTag('pic','desc'=>"注：示例图。DE_miRNA、DE_mRNA_TargetmiRNA、DE_LncRNA_TargetmiRNA、DE_circRNA_TargetmiRNA分别指差异miRNA、差异基因所靶向的所有miRNA、差异lncRNA所靶向的所有miRNA、差异circRNA所靶向的所有miRNA。",'name'=>"图$pid. 差异miRNA、差异基因所靶向的所有miRNA、差异lncRNA所靶向的所有miRNA、差异circRNA所靶向的所有miRNA的交集（示例）",'type'=>"type1",'path'=>"$template_dir/demo.miRNA.venn.png");
}else{
	$writer->emptyTag('p','desc'=>"差异miRNA，差异基因所靶向的所有miRNA，差异lncRNA所靶向的所有miRNA，差异circRNA所靶向的所有miRNA的交集以Venn图的形式展示。如下：",'type'=>"type1");
	&piclist("差异miRNA、差异基因所靶向的所有miRNA、差异lncRNA所靶向的所有miRNA、差异circRNA所靶向的所有miRNA的交集","注:DE_miRNA、DE_mRNA_TargetmiRNA、DE_LncRNA_TargetmiRNA、DE_circRNA_TargetmiRNA分别指差异miRNA、差异基因所靶向的所有miRNA、差异lncRNA所靶向的所有miRNA、差异circRNA所靶向的所有miRNA","$combine/BMK_4_cytoscape/sRNA/*venn.png","$relate_dir");
}

if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"差异miRNA相关靶向RNA分析结果路径",'action'=>"xls",'path'=>"../../$combine_dir/BMK_4_cytoscape/sRNA",'type'=>"xls");
}

$writer->emptyTag('h2','name'=>"$h1.4 差异gene相关靶向RNA分析",'type'=>'type1','desc'=>"$h1.4 差异gene相关靶向RNA分析");
$writer->emptyTag('p','desc'=>"以差异gene为中心进行靶向关系联合分析，将差异gene及与其具有靶向关系的差异miRNA、差异LncRNA、差异circRNA，整理成表格。结果如下表所示：",'type'=>"type1");
$tid++;
if(-e "$combine/../HTML/$template_dir/diff_gene.xls"){
	$writer->emptyTag('table','desc'=>"注：#ID：差异gene的ID；*_Count：count数；*_FPKM：FPKM值；Pvalue：Pvalue值；log2FC：log2FC值；regulated：上下调（up/down）；diff_target_miRNA：有gene有靶向关系的差异miRNA ID；Cis.diff_target_LncRNA：与gene有顺式靶向关系的差异LncRNA ID；Trans.diff_target_LncRNA：与gene有反式靶向关系的差异LncRNA ID（当样品数目不少于5个时）；diff_target_circRNA：与gene有靶向关系的差异circRNA ID",'name'=>"表$tid. 差异gene靶向关系结果",'type'=>"part",'path'=>"$template_dir/diff_gene.xls");
}else{
	$writer->emptyTag('table','desc'=>"注：示例表格。#ID：差异gene的ID；*_Count：count数；*_FPKM：FPKM值；Pvalue：Pvalue值；log2FC：log2FC值；regulated：上下调（up/down）；diff_target_miRNA：有gene有靶向关系的差异miRNA ID；Cis.diff_target_LncRNA：与gene有顺式靶向关系的差异LncRNA ID；Trans.diff_target_LncRNA：与gene有反式靶向关系的差异LncRNA ID（当样品数目不少于5个时）；diff_target_circRNA：与gene有靶向关系的差异circRNA ID",'name'=>"表$tid. 差异gene靶向关系结果（示例）",'type'=>"part",'path'=>"$template_dir/demo/demo.diff_gene.xls");
}
$writer->emptyTag('p','desc'=>"各组合差异gene靶向关系结果：",'type'=>"type1");
&filelist("","$combine/BMK_4_cytoscape/gene/*_vs_*.DEG.xls","$combine_dir");

my @gene_png = glob "$combine/BMK_4_cytoscape/gene/*venn.png";
if (@gene_png==0){
	if($sample_num>=5){
		$writer->emptyTag('p','desc'=>"差异mRNA，差异lncRNA的所有靶基因（当样品数目不少于5时，分顺式和反式分别展示），差异miRNA的所有靶基因，差异circRNA的所有来源基因的交集以Venn图的形式展示。示例图如下:",'type'=>"type1");
	        &piclist("差异mRNA、差异lncRNA的所有靶基因、差异miRNA的所有靶基因、差异circRNA的所有来源基因的交集（示例）","注：示例图。DE_mRNA、DE_LncRNA_Target.CismRNA、DE_LncRNA_Target.TransmRNA(样品数目不少于5)、DE_miRNA_TargetmRNA、DE_circRNA_Hostgene分别指差异mRNA、差异lncRNA的顺式/反式靶基因、差异miRNA的所有靶基因、差异circRNA的所有来源基因","$combine/../HTML/$template_dir/demo.mRNA.*.venn.png","$relate_dir");
	}else{
		$writer->emptyTag('p','desc'=>"差异mRNA，差异lncRNA的所有靶基因（当样品数目少于5时，只有顺式的结果），差异miRNA的所有靶基因，差异circRNA的所有来源基因的交集以Venn图的形式展示。示例图如下：",'type'=>"type1");
		$pid++;
	        $writer->emptyTag('pic','desc'=>"注：示例图。DE_mRNA、DE_LncRNA_Target.CismRNA、DE_miRNA_TargetmRNA、DE_circRNA_Hostgene分别指差异mRNA、差异lncRNA的顺式靶基因、差异miRNA的所有靶基因、差异circRNA的所有来源基因",'name'=>"图$pid. 差异mRNA、差异lncRNA的所有靶基因、差异miRNA的所有靶基因、差异circRNA的所有来源基因的交集（示例）",'type'=>"type1",'path'=>"$template_dir/demo.mRNA.Cis.venn.png");
	}
}else{
	$writer->emptyTag('p','desc'=>"差异mRNA，差异lncRNA的所有靶基因（当样品数目不少于5时，分顺式和反式分别展示），差异miRNA的所有靶基因，差异circRNA的所有来源基因的交集以Venn图的形式展示。如下：",'type'=>"type1");
	&piclist("差异mRNA、差异lncRNA的所有靶基因、差异miRNA的所有靶基因、差异circRNA的所有来源基因的交集","注：DE_mRNA、DE_LncRNA_Target.CismRNA、DE_LncRNA_Target.TransmRNA（样品数目不少于5）、DE_miRNA_TargetmRNA、DE_circRNA_Hostgene分别指差异mRNA、差异lncRNA的所有靶基因、差异miRNA的所有靶基因、差异circRNA的所有来源基因","$combine/BMK_4_cytoscape/gene/*venn.png","$relate_dir");
}

if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"差异gene相关靶向RNA分析结果路径",'action'=>"xls",'path'=>"../../$combine_dir/BMK_4_cytoscape/gene",'type'=>"xls");
}

$writer->emptyTag('h2','name'=>"$h1.5 基于ceRNA网络的差异lncRNA-miRNA-mRNA关系对分析",'type'=>'type1','desc'=>"$h1.5 基于ceRNA网络的差异lncRNA-miRNA-mRNA关系对分析");
$writer->emptyTag('p','desc'=>"基于ceRNA网络，对每个差异组合提取三者都差异的关系对，整理成表。结果如下：",'type'=>"type1");
$tid++;
if(-e "$combine/../HTML/$template_dir/lncRNA-miRNA-mRNA.xls"){
	$writer->emptyTag('table','desc'=>"注：#第1列为ceRNA1，第2列为ceRNA2，第3列为ceRNA1所结合的差异miRNA数目，第4列为ceRNA2所结合的差异miRNA数目，第5列为ceRNA1和ceRNA2共同结合的差异miRNA数目，第6列为ceRNA1所结合的差异miRNA，第7列为ceRNA2所结合的差异miRNA，第8列为ceRNA1和ceRNA2共同结合的差异miRNA数目。这里，ceRNA1、ceRNA2可由差异lncRNA、差异mRNA充当。",'name'=>"表$tid. 基于ceRNA网络的差异lncRNA-miRNA-mRNA关系对",'type'=>"part",'path'=>"$template_dir/lncRNA-miRNA-mRNA.xls");
}else{
	$writer->emptyTag('table','desc'=>"注：示例表格。#第1列为ceRNA1，第2列为ceRNA2，第3列为ceRNA1所结合的差异miRNA数目，第4列为ceRNA2所结合的差异miRNA数目，第5列为ceRNA1和ceRNA2共同结合的差异miRNA数目，第6列为ceRNA1所结合的差异miRNA，第7列为ceRNA2所结合的差异miRNA，第8列为ceRNA1和ceRNA2共同结合的差异miRNA数目。这里，ceRNA1、ceRNA2可由差异lncRNA、差异mRNA充当。",'type'=>"part",'name'=>"表$tid. 基于ceRNA网络的差异lncRNA-miRNA-mRNA关系对（示例）",'path'=>"$template_dir/demo/demo.lncRNA-miRNA-mRNA.xls");
}
&filelist("","$combine/BMK_4_cytoscape/lncRNA-miRNA-mRNA/*_vs_*.lncRNA-miRNA-mRNA.ceRNA_pair_adjust_p_Sig_diff.xls","$combine_dir");
$writer->emptyTag('p','desc'=>"Cytoscape示例图如下:",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"注：示例图。图中，lncRNA、miRNA、mRNA分别用圆、菱形、箭头表示，上下调分别用红色和绿色表示。",'name'=>"图$pid. 差异lncRNA-miRNA-mRNA关系对（示例） ",'type'=>"type1",'path'=>"$template_dir/demo.lncRNA-miRNA-mRNA.png");

if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"基于ceRNA网络的差异lncRNA-miRNA-mRNA关系对分析结果路径",'action'=>"xls",'path'=>"../../$combine_dir/BMK_4_cytoscape/lncRNA-miRNA-mRNA",'type'=>"xls");
}

$writer->emptyTag('h2','name'=>"$h1.6 基于ceRNA网络的差异circRNA-miRNA-mRNA关系对分析",'type'=>'type1','desc'=>"$h1.5 基于ceRNA网络的差异circRNA-miRNA-mRNA关系对分析");
$writer->emptyTag('p','desc'=>"基于ceRNA网络，对每个差异组合提取三者都差异的关系对，整理成表。结果如下：",'type'=>"type1");
$tid++;
if(-e "$combine/../HTML/$template_dir/circRNA-miRNA-mRNA.xls"){
	$writer->emptyTag('table','desc'=>"注：#第1列为ceRNA1，第2列为ceRNA2，第3列为ceRNA1所结合的差异miRNA数目，第4列为ceRNA2所结合的差异miRNA数目，第5列为ceRNA1和ceRNA2共同结合的差异miRNA数目，第6列为ceRNA1所结合的差异miRNA，第7列为ceRNA2所结合的差异miRNA，第8列为ceRNA1和ceRNA2共同结合的差异miRNA数目。这里，ceRNA1、ceRNA2可由差异circRNA、差异mRNA充当。",'name'=>"表$tid. 基于ceRNA网络的差异circRNA-miRNA-mRNA关系对",'type'=>"part",'path'=>"$template_dir/circRNA-miRNA-mRNA.xls");
}else{
	$writer->emptyTag('table','desc'=>"注：示例表格。#第1列为ceRNA1，第2列为ceRNA2，第3列为ceRNA1所结合的差异miRNA数目，第4列为ceRNA2所结合的差异miRNA数目，第5列为ceRNA1和ceRNA2共同结合的差异miRNA数目，第6列为ceRNA1所结合的差异miRNA，第7列为ceRNA2所结合的差异miRNA，第8列为ceRNA1和ceRNA2共同结合的差异miRNA数目。这里，ceRNA1、ceRNA2可由差异circRNA、差异mRNA充当。",'type'=>"part",'name'=>"表$tid. 基于ceRNA网络的差异circRNA-miRNA-mRNA关系对（示例）",'path'=>"$template_dir/demo/demo.circRNA-miRNA-mRNA.xls");
}
&filelist("","$combine/BMK_4_cytoscape/circRNA-miRNA-mRNA/*_vs_*.circRNA-miRNA-mRNA.ceRNA_pair_adjust_p_Sig_diff.xls","$combine_dir");
$writer->emptyTag('p','desc'=>"Cytoscape示例图如下:",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"注：示例图。图中，circRNA、miRNA、mRNA分别用圆、菱形、箭头表示，上下调分别用红色和绿色表示。",'name'=>"图$pid. 差异circRNA-miRNA-mRNA关系对（示例）",'type'=>"type1",'path'=>"$template_dir/demo.circRNA-miRNA-mRNA.png");

if(!defined $cloud){
	$writer->emptyTag('file','desc'=>"",'name'=>"基于ceRNA网络的差异circRNA-miRNA-mRNA关系对分析结果路径",'action'=>"xls",'path'=>"../../$combine_dir/BMK_4_cytoscape/circRNA-miRNA-mRNA",'type'=>"xls");
}

$writer->startTag('ref_list',desc=>"参考文献",type=>"type1",name=>"参考文献");
if($sample_num<4){
&reference ($writer,[
	["Salmena L, Poliseno L, Tay Y, Kats L, Pandolfi PP. A ceRNA hypothesis: the Rosetta Stone of a hidden RNA language? Cell. 2011;146(3):353–358.","https://www.ncbi.nlm.nih.gov/pubmed/?term=A%20ceRNA%20hypothesis%3A%20the%20Rosetta%20Stone%20of%20a%20hidden%20RNA%20language%3F","1"],
]);
}else{
	&reference ($writer,[
	["Dou C, Cao Z, Yang B, et al. Changing expression profiles of lncRNAs, mRNAs, circRNAs and miRNAs during osteoclastogenesis.[J]. Sci Rep, 2016, 6:21499.","https://www.ncbi.nlm.nih.gov/pubmed/?term=Changing%20expression%20profiles%20of%20lncRNAs%2C%20mRNAs%2C%20circRNAs%20and%20miRNAs%20during%20osteoclastogenesis","1"],
        ["Salmena L, Poliseno L, Tay Y, Kats L, Pandolfi PP. A ceRNA hypothesis: the Rosetta Stone of a hidden RNA language? Cell. 2011;146(3):353–358.","https://www.ncbi.nlm.nih.gov/pubmed/?term=A%20ceRNA%20hypothesis%3A%20the%20Rosetta%20Stone%20of%20a%20hidden%20RNA%20language%3F","2"],
]);
}

$writer->endTag('ref_list');

$writer->endTag('report');
open OUT,">:utf8", "$output";
my $xmlstr=&decorate($writer->to_string);
$xmlstr = Encode::decode("utf-8", $xmlstr);
print OUT $xmlstr;
close(OUT);
$writer->end();

################################################################################################
#		                        sub function
################################################################################################
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
		my $tmp=(split(/$path\//,$dir))[1];
        	$writer->emptyTag('pic','desc'=>"",'name'=>"$base",'type'=>"type1",'path'=>"../$tmp/$base");
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
		my $new_path = "$path"."$tmp";
       	$writer->emptyTag('file','desc'=>"",'name'=>"$base",'action'=>"xls",'path'=>"$relate_cloud/$new_path/$base",'type'=>"xls");
	}
	$writer->endTag('file_list');
}

sub readcfg{
	my %config=();
	open(CFG,$cfg)||die $!;
	while(<CFG>){
		chomp;next if($_=~/^#|^\s+|^$/);
		my @tmp=split(/\s+/,$_);
		$config{$tmp[0]}=$tmp[1];				
	}
	close(CFG);
	my @para=("fold","FDR","PValue","Method_RE","Method_NR");
	foreach my $pa(@para){
		$config{"$type.$pa"}=$config{$pa}	if(exists $config{$pa} && !exists $config{"$type.$pa"});
	}
	return %config;
}

sub readConfig{
        my $configFile=shift;
        my $d=Config::General->new(-ConfigFile => "$configFile");
        my %config=$d->getall;
        return %config;
}

sub gaintime{
	my $timestamp=time(); 
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($timestamp); 
	my $y = $year + 1900; 
	my $m = $mon + 1; 
	$timestamp=sprintf("%4d-%02d-%02d",$y,$m,$mday);
	return $timestamp
}
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
sub table_html{
        my ($input,$outHtml,$title,$linkColNum,$linkHash,$srcPath,$text)=@_;
        my @inputs=();
        $/="\n";
        my $i=0;
        open (IN,$input)||die $!;
        while(<IN>){
                chomp;
                next if /^\s*$/;
                my @tmp=split/\t+/;
                push@inputs,\@tmp;
                $i++;
        }
        $/="》";
        open HH,">$outHtml" or die "$!";
	my $titles=basename $outHtml;
	$titles=~s/\.html//;
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
#print "$titles\t$inputs[$k][$i]\t$$linkHash{$titles}{$inputs[$k][$i]}\n";
                                                print HH "<td><a href=\"$$linkHash{$titles}{$inputs[$k][$i]}\" target=\"_blank\">$inputs[$k][$i]</a></td>";
                                        }else{
                                                print HH "<td>$inputs[$k][$i]</td>";
                                        }
                                }else{
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
}



sub USAGE{
        my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
        Options:
	-combine	<path>	input BMK_Combine_Analysis path, forced
	-data_cfg	<file>	data_cfg	
	-detail_cfg	<file>	input detail cfg
	-type		<str>	combine
	-o		<file>	output xml file
	-cloud		<null>	generate xml for cloud
        -h      Help

Example: perl $0 -combine BMK_7_Combine -data_cfg data.cfg -detail_cfg detail.cfg -type combine -o HTML/combine.xml 

USAGE
        print $usage;
        exit;
}

