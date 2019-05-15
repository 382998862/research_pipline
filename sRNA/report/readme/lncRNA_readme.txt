.
|-- BMK_1_Assembly_Result
	|-- L*.StringTie.transcripts.gtf	StringTie软件组装转录本结果
		|-- 第一列:Chr:染色体序号
		|-- 第二列:Source:来源
		|-- 第三列:Type:转录本或者其外显子
		|-- 第四列:Start Site:起始坐标
		|-- 第五列: End Site:终止坐标
		|-- 第六列:Score:对应位置得分
		|-- 第七列:Strand:链的信息
		|-- 第八列:frame:frame 的信息
		|-- 第九 列:Description:描述信息
|-- BMK_2_LncRNA_Prediction
	|-- BMK_1_Software_Result	候选 lncRNA 的四种软件预测结果
		|-- CNCI.xls
		|-- CPAT.xls
		|-- CPC.xls
		|-- Pfam.xls
		|-- Software_veen.xls
		|-- venn.png
	|-- LncRNA_circos.png	lncRNA 在染色体上分布环状图,描述不同种类lncRNA 在染色体上的分布
注:图片的最外层是基因组的染色体环,由外到内依次是正义链 lncRNA 环(绿色),基因间区 lncRNA 环(红色),反义 lncRNA 环(灰色),内含子 lncRNA 环(蓝色).
	|-- LncRNA_classification.png	lncRNA 分类柱状图,横坐标为 4种不同类型的 LncRNA(其中:lincRNA:基因间区长链非编码 RNA,Antisense-lncRNA:反义长链非编码 RNA,Intronic-lncRNA:内含子长链非编码 RNA,sense_lncRNA:正义长链非编码 RNA), 纵坐标为对应的 LncRNA 数量
	|-- LncRNA_classification.xls	lncRNA 对应的class code 类型
	|-- LncRNA.fa	lncRNA 预测结果 fasta 文件(重要文件)
注:FASTA 格式每一个序列单元以“>”开头,直到出现下一个“>”之前为止。“>”开头的行为序列 ID 行,后面紧接着基因 ID;下面一行或多行为该基因的碱基序列。
	|-- LncRNA.gff	lncRNA 预测结果 gff 文件(重要文件)
		|-- 第一列:#Seq_ID:染色体号
		|-- 第二列:Source:注释信息的来源
		|-- 第三列:Type:注释特征(Feature)类型
		|-- 第四,五列:Start/End:特 征序列的起止位置
		|-- 第六列:Score:得分,数字,注释信息可能性的说明,“.”表示缺失值
		|-- 第七列:Strand:特征序列所在的正负链
		|-- 第八列:Phase:仅对注释类型为 CDS 有效,表示起始编码的位置,有效值为0、1、2,“.”表示缺失值
		|-- 第九列:Attributes:以多个键值对组成的注释信息描述
	|-- LncRNA.gtf	lncRNA 预测结果 gtf 文件
		|-- 第一列:Chr:染色体序号
		|-- 第二列:Source:来源
		|-- 第三列:Type:转录本或者其外显子
		|-- 第四列:Start Site:起始坐标
		|-- 第五列:End Site:终止坐标
		|-- 第六列:Score:对应位置得分
		|-- 第七列:Strand:链的信息
		|-- 第八列:frame:frame 的信息
		|-- 第九列:Description:描述信息
|-- BMK_3_LncRNA_Expression
	|-- lncRNA_counts.xls	所有 lncRNA 的表达量(比对片段数)矩阵文件(重要文件)
	|-- lncRNA_expression.xls	所有 lncRNA 的表达量(fpkm)矩阵文件(重要文件)
	|-- PNG
		|-- all.fpkm_box.png	各样品表达量箱线图,图中横坐标代表不同的样品，纵坐标表示样品表达量FPKM的对数值。该图从表达量的总体离散角度来衡量各样品表达水平。
		|-- all.fpkm_density.png	各样品表达量密度分布图,图中不同颜色的曲线代表不同的样品，横坐标表示对应样品FPKM的对数值，纵坐标表示概率密度。
		|-- *.fpkm_density.png	单个样品表达量密度分布图
		|-- sample_coefficient.txt	样本间相关性系数表
		|-- sample_corelation.pdf	样本间相关性系数图
		|-- sample_corelation.png
		|-- sample_pvalue.txt
|-- BMK_4_LncRNA_Target
	|-- Cis_target_gene.xls	lncRNA 顺式靶基因预测结果文件(重要文件)
	|-- Trans_target_gene.xls	(当样品数目不少于5时)lncRNA 反式靶基因预测结果文件,通过样本间lncRNA 与 mRNA 的表达量相关性分析方法来预测lncRNA的反式靶基因。此方法需要一定样品数量保证其预测的准确性，因此当样本量小于5时不进行此项分析。采用Pearson 相关系数法分析样本间 lncRNA 与 mRNA 的相关性，取相关性绝对值大于0.9且显著性p值小于0.01的基因作为lncRNA的反式靶基因.(重要文件)
	|--WGCNA_Result	共表达分析结果目录,样品数少于5个时没有共表达分析结果,只有PCA分析部分结果。结果文件说明见文件夹内部的readme文件
|-- BMK_5_DEG_Analysis
	|-- BMK_1_All_DEG
		|-- All.DEG_final_Cis_anno.xls	所有差异LncRNA/表达量/顺式靶基因/顺式靶基因注释整合结果(整合文件,信息全面,使用方便,比较重要)
			|-- #Gene:顺式靶基因ID
			|-- ID:差异LncRNA ID
			|-- 样本ID:表达量列
			|-- *_vs_*.Pvalue/FDR
			|-- *_vs_*.log2FC
			|-- *_vs_*.regulated
		|-- All.DEG_final.pdf
		|-- All.DEG_final.png	所有差异LncRNA表达量聚类热图,横坐标代表样品名称及样品的聚类结果，纵坐标代表的差异lncRNA 及lncRNA的聚类结果。图中不同的列代表不同的样品，不同的行代表不同的基因。颜色代表了基因在样品中的表达量水平log10（lncRNA+0.000001）
		|-- All.DEG_final_Trans_anno.xls	所有差异LncRNA/表达量/反式靶基因/反式靶基因注释整合结果(样本数量不少于5)
		|-- All.DEG_final.xls	所有差异LncRNA表达量/Pvalue或者FDR/regulated
		|-- Veen.pdf
		|-- Veen.png	(9>差异组合数目>=2)
	|-- BMK_2_*_vs_*	条件 1 样品和条件 2 样品的差异表达分析结果目录(*表示比较的样本编号,如 BMK_2_L01_L02_vs_L03_L04 表示L01,L02 与 L03,L04 做差异比较)。在注释统计表中不存在的内容表示没有对该数据库做注释,不存在的文件夹表示没有对该项内容做分析
		|-- Cis_Anno_enrichment	差异LncRNA顺式靶基因注释
	      	|-- anno
        		|-- *_vs_*.annotation.xls	差异基因顺式靶基因注释信息和差异LncRNA信息整合文件
		        |-- *_vs_*.GO_enrichment.stat.xls	GO二级节点注释分类统计结果文件
		        |-- *_vs_*.GO.pdf	GO二级节点注释分类统计图
		        |-- *_vs_*.GO.png
		|-- enrich	差异表达lncRNA顺式靶基因GO富集分析/KEGG富集分析分析.此部分分析是采用R包clusterProfiler对lncRNA顺式靶基因分别进行生物学过程，分子功能和细胞组分的富集分析以及KEGG通路富集分析。富集分析采用超几何检验方法来寻找与整个基因组背景相比显著富集的GO条目/KEGG通路。此目录下的BP/CC/MF分别是Biological_Process/Cellular_Component/Molecular_Function的简写
			|-- BMK_1_GO_enrich
				|-- *_vs_*_(BP/CC/MF)_enrich_barplot.pdf	
				|-- *_vs_*_(BP/CC/MF)_enrich_barplot.png	差异表达lncRNA顺式靶基因GO富集条形图
				|-- *_vs_*_(BP/CC/MF)_enrich_dotplot.pdf
				|-- *_vs_*_(BP/CC/MF)_enrich_dotplot.png	差异表达lncRNA顺式靶基因GO富集图(BP分支),横坐标为GeneRatio即注释在该条目中的感兴趣基因占所有差异lncRNA顺式靶基因数的比例, 纵坐标每一个BP条目。点的大小代表该通路中注释的差异表达lncRNA顺式靶基因数，点的颜色代表超几何检验的校正后的p值。
				|-- *_vs_*_Biological_Process_enrich.list	富集结果文件
					|--ID      GO节点
					|--Description   GO节点名称  
					|--GeneRatio	富集到此节点上的差异靶基因的数目/所有差异靶基因的数目
					|--BgRatio	富集到此节点上的所有靶基因的数目/所有靶基因的数目
					|--enrich_factor 富集因子  
					|--pvalue  Pvalue
					|--qvalue  Qvalue
					|--geneID  富集到此通路上的差异靶基因的ID
					|--gene_symbol	对应的Name
			|-- BMK_2_KEGG_enrich
				|-- BMK_1_Keggmap	KEGG通路
					|-- ko*.html
					|-- ko*.png 相对于对照组来说，红色框标记的酶与上调基因有关，绿色框标记的酶与下调基因有关。蓝色框标记的酶与上调和下调基因均有关，框内的数字代表酶的编号（EC number），而整个通路由多种酶催化的复杂生化反应构成，此通路图中与差异表达基因相关的酶均用不同的颜色标出，研究人员可以根据自己的研究对象间的差异，重点研究某些代谢通路相关基因的差异表达情况，通过通路解释表型差异的根源。
			        |-- *_vs_*_KEGG_pathway_enrich_barplot.pdf
			        |-- *_vs_*_KEGG_pathway_enrich_barplot.png	差异表达lncRNA顺式靶基因KEGG富集条形图,注：横坐标代表差异表达lncRNA顺式靶基因注释在通路中的基因数，纵坐标代表通路，柱的颜色代表校正后的p值。
			        |-- *_vs_*_KEGG_pathway_enrich_dotplot.pdf
		        	|-- *_vs_*_KEGG_pathway_enrich_dotplot.png	差异表达lncRNA顺式靶基因KEGG富集图
			        |-- *_vs_*_KEGG_pathway_enrich.list	同BP
				|-- 	KEGG通路
		|-- *_vs_*.all	该组合中所有LncRNA的表达量/log2FC/Pvalue/FDR
		|-- *_vs_*.DEG_final.Cis_Target.xls	差异LncRNA增加顺式靶基因ID
		|-- *_vs_*.DEG_final.Trans_Target.xls	差异LncRNA增加反式靶基因ID
		|-- *_vs_*.DEG_final.xls	差异LncRNA的信息
		|-- *_vs_*.heatmap.pdf
		|-- *_vs_*.heatmap.png	差异表达lncRNA聚类热图.
					注：横坐标代表样品名称及样品的聚类结果，纵坐标代表的差异lncRNA及lncRNA的聚类结果。图中不同的列代表不同的样品，不同的行代表不同的基因。颜色代表了基因在样品中的表达量水平log10（lncRNA+0.000001）。
		|-- *_vs_*_MA.pdf
		|-- *_vs_*_MA.png	差异表达MA图.
					注：差异表达基因MA图中每一个点代表一个lncRNA。横坐标为A值：log2(FPKM)，即两样品中表达量均值的对数值；纵坐标为M值：log2(FC)，即两样品间lncRNA表达量差异倍数的对数值，用于衡量表达量差异的大小。图中绿色的点代表显著下调的lncRNA,红色的点代表显著上调的lncRNA, 黑色的点代表表达差异不显著的lncRNA。
		|-- *_vs_*_Volcano.pdf
		|-- *_vs_*_Volcano.png	差异表达火山图.
					注：差异表达火山图中的每一个点表示一个lncRNA,横坐标表示某一个基因在两样品中表达量差异倍数的对数值;纵坐标表示pvalue的负对数值。横坐标绝对值越大，说明表达量在两样品间的表达量倍数差异越大;纵坐标值越大，表明差异表达越显著，筛选得到的差异表达基因越可靠。图中绿色的点代表下调差异表达lncRNA,红色的点代表上调差异表达lncRNA, 黑色代表非差异表达lncRNA。
	|-- Trans_Anno_enrichment	(样品数目不少于5个时有此部分结果)反式靶基因注释信息,同顺式靶基因
		|-- anno
			|-- *_vs_*.annotation.xls
			|-- *_vs_*.GO_enrichment.stat.xls
			|-- *_vs_*.GO.pdf
			|-- *_vs_*.GO.png
		|-- enrich
			|-- BMK_1_GO_enrich
				|-- *_vs_*_(BP/CC/MF)_enrich_barplot.pdf
				|-- *_vs_*_(BP/CC/MF)_enrich_barplot.png
				|-- *_vs_*_(BP/CC/MF)_enrich_dotplot.pdf
				|-- *_vs_*_(BP/CC/MF)_enrich_dotplot.png	
				|-- *_vs_*_(BP/CC/MF)_enrich.list
			|-- BMK_2_KEGG_enrich
                	        |-- BMK_1_Keggmap
        	                        |-- ko*.html
	                                |-- ko*.png
				|-- *_vs_*_KEGG_pathway_enrich_barplot.pdf
				|-- *_vs_*_KEGG_pathway_enrich_barplot.png
				|-- *_vs_*_KEGG_pathway_enrich_dotplot.pdf
				|-- *_vs_*_KEGG_pathway_enrich_dotplot.png
				|-- *_vs_*_KEGG_pathway_enrich.list
	|-- DEG_Cis_anno.stat	差异LncRNA的顺式靶基因注释统计结果
	|-- DEG.stat	差异LncRNA统计结果
	|-- DEG_Trans_anno.stat	差异LncRNA的反式靶基因注释统计结果
|-- BMK_6_Lnc_disease
	|-- LncRNA_disease.xls	LncRNA疾病注释结果
		|-- #LncRNA_id,项目中lncRNA ID
		|-- gene_symbol:lncRNA gene symbol
		|-- DB_lncRNA，数据库中比对上的lncRNA
		|-- Source，lncRNA来源
		|-- disease：相关的疾病
		|-- regulate：调控方式
		|-- description：疾病相关描述
		|-- chr：染色体
		|-- start：起始
		|-- end：终止
		|-- strand：正负链
		|-- pmid：文章
|-- BMK_7_Lnc_precursor
	|-- lncRNA_precursor.xls	LncRNA前体分析结果
		|-- LncRNA_ID：LncRNA ID
		|-- gene_symbol:lncRNA gene symbol
		|-- miRNA_hp：miRNA前体ID
		|-- LncRNA_start：LncRNA比对起始位置
		|-- LncRNA_end：LncRNA比对终止位置
