结构分析结果文件:
|-- BMK_1_SNP_Analysis	SNP分析结果
	|-- AllSample.SNP_density.png	SNP 类型密度图.注:横轴为基因上平均每1000bp 序列中分布的 SNP 数目,纵轴为基因数。
	|-- AllSample.SNP_density.stat	SNP 类型密度统计结果文件
		|-- #Sample:样品名称
		|-- Interval:SNP 密度区间(个/ Kbp ) 
		|-- GeneNum:在相应(第二列)密度区间的基因数
	|-- AllSample.snp.stat	SNP 位点统计表
		|-- #BMK-ID:百迈客对样品的统一编号
		|-- SNP Number: SNP 位点总数
		|-- Genic SNP:基因区 SNP 位点总数
		|-- Intergenic SNP:基因间区 SNP 位点总数
		|-- Transition:转换类型的 SNP 位点数目在总 SNP 位点数目中所占的百分比
		|-- Transversion:颠换类型的 SNP 位点数目在总 SNP 位点数目中所占的百分比
		|-- Heterozygosity:杂合型 SNP 位点数目在总 SNP 位点数目中所占的百分比
	|-- BMK_1_InDel_anno	Indel 注释结果文件夹
		|-- all.indel.anno.stat.png	Indel 注释分类图
		|-- final_Indel.anno.stat	Indel 注释分类统计文件
		|-- *.indel.anno.stat.png	单个样本Indel 注释分类图
	|-- BMK_2_SNP_anno	SNP 注释结果文件夹
		|-- all.snp.anno.stat.png	SNP 注释分类图
		|-- final_SNP.anno.stat		SNP 注释分类统计文件,常见 SNP 类型分类见结题报告2.6.1
		|-- *.snp.anno.stat.png		单个样本 SNP 注释分类图
	|-- BMK_3_SNP_type	SNP 类型结果文件夹
		|-- All.snp.type.png	SNP 突变类型分布图,横轴为 SNP 突变类型，纵轴为相应的 SNP 数目
		|-- All.snp_type.stat	SNP 突变类型统计文件
		|-- *.snp.type.png	单个样本 SNP 突变类型分布图
	|-- final.*.anno.gatk.all.list	SNP/Indel位点信息文件
		|-- #Chr:SNP/InDel 位点所在染色体编号
		|-- Pos:SNP/InDel 位点在染色体上的位置
		|-- Gene_id:SNP/InDel 位点所在的基因或原来未注释的基因区(表中用Intergenic表示)
		|-- Ref:所选参考基因组中的 SNP/InDel 等位
		|-- Alt:测序样品中识别到的 SNP/InDel
		|-- **:样品 ** 中 SNP/InDel 位点的分型, SNP 类型见附录
		|-- Depth:样品 ** 中 SNP/InDel 位点的测序深度
		|-- AlleDp:样品 ** 中 SNP/InDel 位点各等位点的支持深度
		|-- Effect:SNP/InDel 所在区域或类型, Effect 具体说明详见:http://snpeff.sourceforge.net/SnpEff_manual.html
		|-- Codon_change:编码改变方式,未改变用点表示
########
其中Alt列中的 SNP 位点的碱基缩写表见下:
Nucleic_Acid_Code	Meaning	Mnemonic
A	A	Adenine
C	C	Cytosine
G	G	Guanine
T	T	Thymine
U	U	Uracil
R	A or G	puRine
Y	C,T or U	pYrimidines
K	G,T or U 	bases which are Ketones
M	A or C	bases with aMino groups
S	C or G	Strong interaction
W	A,T or U	Weak interaction
B	not A (i.e. C,G,T or U)	B comes after A
D	not C (i.e. A,G,T or U)	D comes after C
H	not G (i.e. A,C,T or U)	H comes after G
V	neither T nor U (i.e. A, C or G)	V comes after U
N	A C G T U	Nucleic acid
########
|-- BMK_2_miRNA_Target	(重要文件夹)miRNA靶向关系结果.植物用TargetFinder软件进行靶基因预测;动物用miRanda和targetscan进行靶基因预测,最终取两款软件的交集作为最终靶向关系对
	|-- circRNA/gene/lncRNA.mir2target.list	miRNA与circRNA/gene/lncRNA靶向关系结果
		|--#ID	miRNA ID
		|--Target	与miRNA有靶向关系的RNA的ID,多个RNA用";"分隔
	|-- circRNA/gene/lncRNA.mir2target.stat	统计结果,第一行是有靶向关系的miRNA的个数,第二行是有靶向关系的RNA的个数
	|-- circRNA/gene/lncRNA.miRNA_target_info.xls
		|-- #miRNA	miRNA ID
		|-- target	靶向RNA ID
		|-- miRanda	0表示该靶向关系对在该软件中并未预测到，1表示该靶向关系对在该软件中预测到
		|-- targetscan	0表示该靶向关系对在该软件中并未预测到，1表示该靶向关系对在该软件中预测到
|-- BMK_3_Gene_Structure_Optimize
	|-- Human.geneStructure.optimize.xls	基因结构优化结果
		|-- GeneID:基因ID
		|-- Locus:基因座，格式为“染色体编号:起点坐标-终点坐标”
		|-- Strand:正负链
		|-- Site:优化的位置,3'或5'UTR
		|-- OriginalRegion:原来注释的第一个或最后一个外显子的起止坐标
		|-- OptimizedRegion:延伸之后的第一个或最后一个外显子的起止坐标
|-- BMK_4_Alt_splice	ASprofile软件对StringTie预测出的基因模型对每个样品的12类可变剪切事件分别进行分类和数量统计
	|-- As_event_stat*.png	可变剪切类型统计图.图中横坐标表示属于一种可变剪接的转录本数，纵坐标表示12种可变剪切类型。
	|-- *.AS.list
		|-- event_id	可变剪切事件ID
		|-- event_type	可变剪切事件类型,具体见下行
		|-- gene_id	基因ID
		|-- chrom	所在染色体
		|-- event_start	事件开始位置
		|-- event_end	时间结束位置
		|-- event_pattern	时间模式,即剪切位点的信息
		|-- strand	方向
########
12类可变剪切事件可分为：
(1) TSS: Alternative 5' first exon (transcription start site)	第一个外显子可变剪切
(2) TTS: Alternative 3' last exon (transcription terminal site)	最后一个外显子可变剪切
(3) SKIP: Skipped exon(SKIP_ON,SKIP_OFF pair)	单外显子跳跃
(4) XSKIP: Approximate SKIP (XSKIP_ON,XSKIP_OFF pair)	单外显子跳跃（模糊边界）
(5) MSKIP: Multi-exon SKIP (MSKIP_ON,MSKIP_OFF pair)	多外显子跳跃
(6) XMSKIP: Approximate MSKIP (XMSKIP_ON,XMSKIP_OFF pair)	多外显子跳跃（模糊边界)
(7) IR: Intron retention (IR_ON, IR_OFF pair)	单内含子滞留
(8) XIR: Approximate IR (XIR_ON,XIR_OFF pair)	单内含子滞留（模糊边界）
(9) MIR: Multi-IR (MIR_ON, MIR_OFF pair)	多内含子滞留
(10) XMIR: Approximate MIR (XMIR_ON, XMIR_OFF pair)	多内含子滞留（模糊边界)
(11) AE: Alternative exon ends (5', 3', or both)	可变 5'或3'端剪切
(12) XAE: Approximate AE 可变 5'或3'端剪切（模糊边界）
########


|-- BMK_5_Gene_Fusion	基因融合事件.使用Fusionmap在转录组中研究基因融合事件.Fusionmap首先通过比对到基因组和转录本中双末端(pairend)关系的序列寻找候选的基因融合，然后采用通过与nt等数据库比较，过滤掉假阳性结果。
	|-- *_circos.png	基因融合circos图.红色的线代表同一染色体上发生的融合事件，绿色的线代表不同染色体上发生的融合事件。
	|-- *_FusionReport.xls	基因融合事件统计
		|-- FusionID	融合的唯一标识ID
		|-- *.HISAT_aln.sorted.bam.UniqueMappingPosition1	表示支持该融合的具有唯一断点的Reads数；（注意：是Reads的断点，不是fusion的断点）
		|-- *.HISAT_aln.sorted.bam.UniqueMappingPosition2	
		|-- *.HISAT_aln.sorted.bam.Count	
		|-- Gene1	
		|-- Strand1	
		|-- Chromosome1	
		|-- Start1	 
		|-- End1	
		|-- Gene2	
		|-- Strand2	
		|-- Chromosome2
		|-- Start2
		|-- End2
		|-- Filter
