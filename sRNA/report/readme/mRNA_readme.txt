mRNA分析结果目录
|-- BMK_1_NewGene	新基因预测结果
	|-- BMK_1_NewGene_Anno	新基因注释结果
		|-- All_Database_annotation.xls	新基因整合各个数据库注释结果
		|-- BMK_1_Annotation	NewGene 功能注释分析结果目录
			|-- *.Cog/Kog_class.txt	COG/KOG 数据库注释结果文件
				|-- Gene name:基因名称
				|-- Portein name in COG:COG 中注释的相应蛋白质的名字
				|-- E_value:E 值
				|-- Identify:比对上的碱基数占总比对碱基长度的百分比
				|-- Score:打分值
				|-- Organism:物种名
				|-- COG(KOG) id:COG (KOG)编号
				|-- COG(KOG) class defination:COG(KOG) 分类名称
				|-- Function code:COG(KOG) 功能分类编号
				|-- Functional categories:COG(KOG) 功能类别名称
				|-- Function class defination:COG(KOG) 功能类别定义
			|-- *.eggNOG_class.txt	eggNOG 功能注释文件
				|-- Query:基因名称
				|-- Match:基因比对上的 eggNOG 数据库中的蛋白序列
				|-- eggNOG:eggNOG 数据库的 ID 编号
				|-- score:基因和序列比对的打分
				|-- Functional Category:eggNOG 数据库对应的功能分类编号
				|-- Description:比对上的序列功能描述
				|-- Function class defination:eggNOG 数据库分类类别的功能描述
			|-- *.GO.anno.txt	各基因注释到 GO	数据库的信息统计
				|-- Gene:基因名称
				|-- Number:某基因注释到的 GO 类别数目
				|-- GO_Anno:注释到的 GO 类别
			|-- *.GO.list.txt	每个基因注释到 GO 数据库上的编号
				|-- 基因 ID
				|-- 相应基因注释到的 GO 编号
			|-- *.Kegg.ko	每个基因的 kegg 注释信息
				|-- Gene_id:基因编号
				|-- KO|e_value|Database_Genes|Anno:四列值,分别表示基因的注释到的 K 编号,blastx 比对上基因 E 值,比对上的基因 ID 以及它的功能注释
			|-- *.Kegg.pathway KEGG 上注释到的基因信息
				|-- pathway:pathway 名称
				|-- pathway_id:pathway 编号(一个 ko 编号表示一个通路)
				|-- Gene_number:注释到某 pathway 上的基因数目
				|-- Gene_id:注释到某 pathway 上的基因名称
				|-- KOs:KO 分类信息(一个代谢通路上有哪些 KO 分类,它是蛋白质(酶)的一个分类体系,序列高度相似,并且在同一条通路上有相似功能的蛋白质被归为一组,然后打上 K 标签)
			|-- *.Pfam.anno.txt	所有的基因在 Pfam 数据库中的注释结果
				|-- Gene_ID:基因名称
				|-- Pfam_IDs :比对上的 Pfam 数据库的编号
				|-- Pfam_Description:比对上的 Pfam 数据库的描述
			|-- *.NR/Swissprot.anno.txt	基因在 NR/swissprot 上的注释结果
				|-- SwissprotGeneID(NrGeneID):基因名称
				|-- Database_ID:注释到数据库的基因编号
				|-- E_value:E 值
				|-- IdeNTity:比对上的碱基数占总比对碱基长度的百分比
				|-- Score:打分值
				|-- Annotation:功能注释
		|-- BMK_2_Statistic	
			|-- *.eggNOG.cluster.png	eggNOG注释分类统计图.横坐标为eggNOG各分类内容，纵坐标为基因数目。在不同的功能类中，基因所占多少反映对应时期和环境下代谢或者生理偏向等内容，可以结合研究对象在各个功能类的分布做出科学的解释。
			|-- *.eggNOG.cluster.stat	eggNOG注释分类统计表
				|-- #ID     class ID
				|-- Class_Name    class Name  
				|-- Numbers	基因数目
			|-- *.GO.png	GO 二级节点注释统计图(png格式),所有基因 GO 分类统计结果图,此图展示的是在全部基因背景下 GO 各二级功能的基因比例和个数。横坐标:GO 分类；纵坐标:左边为基因数目所占百分比,右边为基因数目
			|-- *.nr.lib.png	NR 注释中比对上不同物种分布统计,展示了注释到某物种上的基因数目以及占总检测到的基因的比例
			|-- *.nr.lib.stat	NR 注释中比对上不同物种分布统计表
				|-- #Species_Name	物种名称,只展示最多的前9个物种+Other
				|-- Homologous_Number	数目
				|-- Ratio	比例
		|-- BMK_3_KEGG_map	KEGG通路注释图
			|-- ko*.png	
		|-- Function_Annotation.stat.xls	新基因各个数据库的统计结果
	|-- Human.newGene_final.filtered.gff	新基因预测gff文件
		|-- #Seq_ID:染色体号
		|-- Source:注释信息的来源,Cufflinks 软件
		|-- Type:注释特征(Feature)类 型
		|-- Start/End:特征序列的起止位置
		|-- Score:得分,数字,注释信息可能性的说明,“.”表示缺失值
		|-- Strand:特征序列所在的正负 链
		|-- Phase:仅对注释类型为 CDS 有效,表示起始编码的位置,有效值为0、1、2,“.”表示缺失值
		|-- Attributes:以多个键值对组成的注释 信息描述
	|-- Human.newGene.longest_transcript.fa	新基因最长转录本fa文件
|-- BMK_2_geneExpression
	|-- gene_counts.xls	所有基因的表达量(比对片段数)矩阵文件(重要文件)
	|-- gene_expression.xls	所有基因的表达量(fpkm)矩阵文件(重要文件)
	|-- PNG
		|-- all.fpkm_box.png	 各样本的FPKM箱线图
		|-- all.fpkm_density.png	各样本FPKM密度分布对比图
		|-- *.fpkm_density.png	各个样本FPKM密度分布图
		|-- sample_coefficient.txt	相关性系数表
		|-- sample_corelation.pdf
		|-- sample_corelation.png	相关性系数图
		|-- sample_pvalue.txt
|-- BMK_3_DEG_Analysis
	|-- BMK_1_All_DEG	所有差异表达基因表达量和聚类分析,统计结果
		|-- All.DEG_final_anno.xls	所有差异基因表达信息+差异信息+注释信息(整合文件,信息全面,使用方便,比较重要)
			|-- #ID:差异基因ID
			|-- L*:样品表达量
			|-- *_vs_*.Pvalue/FDR:差异Pvalue/FDR值,所有组合都有列出,若无,则用"--"代替
			|-- *_vs_*.log2FC:差异log2FC值,所有组合都有列出,若无,则用"--"代替
			|-- *_vs_*.regulated:差异上下调(up/down),若无,则用"--"代替
			|-- COG/KOG/GO/KEGG/NR/eggNOG/Swissport 各个数据库注释结果
		|-- All.DEG_final.pdf
		|-- All.DEG_final.png	所有差异基因聚类热图
		|-- All.DEG_final.xls	所有差异基因表达信息+差异信息整合
		|-- Veen.pdf
		|-- Veen.png	Veen图(9>差异组合数目>=2)
	|-- BMK_*_*_vs_*	条件 1 样品和条件 2 样品的差异表达分析结果目录注:(*表示比较的样本编号,如 BMK_4_L01_L02 _vs_L03_L04 表示 L01,L02 与 L03,L04 做差异比较)
		|-- Anno_enrichment
			|-- anno	注释统计
				|-- *_vs_*.annotation.xls	差异基因Count+FPKM+Pvalue/FDR+log2FC+regulated+各个数据库注释结果
				|-- *_vs_*.GO_enrichment.stat.xls	差异基因GO二级节点注释统计
					|-- #GO_classify1	GO二级节点ID
					|-- GO_classify2	GO三级节点名称
					|-- All	注释到该节点的所有基因的数目
					|-- DE	注释到该节点的差异基因的数目
				|-- *_vs_*.GO.pdf
				|-- *_vs_*.GO.png	差异表达基因GO注释分类统计图
							注:横坐标为GO分类，纵坐标左边为基因数目所占百分比，右边为基因数目。此图展示的是在差异表达基因背景和全部基因背景下GO各二级功能的基因富集情况，体现两个背景下各二级功能的地位，具有明显比例差异的二级功能说明差异表达基因与全部基因的富集趋势不同，可以重点分析此功能是否与差异相关。
			|-- enrich	富集分析,采用R包clusterProfiler分别进行GO/KEGG富集分析,此目录下的BP/CC/MF分别是Biological_Process/Cellular_Component/Molecular_Function的简写
				|-- BMK_1_GO_enrich
					|-- *_vs_*_(BP/CC/MF)_enrich_barplot.pdf
	                                |-- *_vs_*_(BP/CC/MF)_enrich_barplot.png	差异表达基因GO富集条形图
					|-- *_vs_*_(BP/CC/MF)_enrich_dotplot.pdf
					|-- *_vs_*_(BP/CC/MF)_enrich_dotplot.png	差异表达基因GO富集图(BP分支),横坐标为GeneRatio即注释在该条目中的感兴趣基因占所有差异基因数的比例,纵坐标每一个BP条目。点的大小代表该通路中注释的差异表达基因数，点的颜色代表超几何检验的校正后的p值。
					|-- *_vs_*_(BP/CC/MF)_enrich.list富集结果文件
						|-- ID GO节点
	                                        |-- Description   GO节点名称
        	                                |-- GeneRatio 富集到此节点上的差异靶基因的数目/所有差异靶基因的数目
	                                        |-- BgRatio 富集到此节点上的所有靶基因的数目/所有靶基因的数目
        	                                |-- enrich_factor 富集因子
                	                        |-- pvalue Pvalue
                        	                |-- qvalue Qvalue
                                	        |-- geneID 富集到此通路上的差异靶基因的ID
	                                        |-- gene_symbol  对应的Name
				|-- BMK_2_KEGG_enrich
					|-- BMK_1_Keggmap KEGG通路注释图
						|-- ko*.html
						|-- ko*.png 相对于对照组来说，红色框标记的酶与上调基因有关，绿色框标记的酶与下调基因有关。蓝色框标记的酶与上调和下调基因均有关，框内的数字代表酶的编号（EC number），而整个通路由多种酶催化的复杂生化反应构成，此通路图中与差异表达基因相关的酶均用不同的颜色标出，研究人员可以根据自己的研究对象间的差异，重点研究某些代谢通路相关基因的差异表达情况，通过通路解释表型差异的根源。
					|-- *_vs_*_KEGG_pathway_enrich_barplot.pdf
	                                |-- *_vs_*_KEGG_pathway_enrich_barplot.png 同GO
	                                |-- *_vs_*_KEGG_pathway_enrich_dotplot.pdf
        	                        |-- *_vs_*_KEGG_pathway_enrich_dotplot.png 同GO
	                                |-- *_vs_*_KEGG_pathway_enrich.list     同GO
        	                        |-- *_vs_*_KEGG_pathway_GSEA.xls        同GO
                	                |-- *_vs_*.Phase.png 差异表达基因KEGG通路富集散点图,图中每一个图形表示一个KEGG通路，通路名称见右侧图例。横坐标为富集因子（Enrichment Factor），表示差异基因中注释到该通路的基因比例与所有基因中注释到该通路的基因比例的比值。富集因子越大，表示差异表达基因在该通路中的富集水平越显著。纵坐标为-log10(Q-value)，其中Q-value为多重假设检验校正之后的P-value。因此，纵坐标越大，表示差异表达基因在该通路中的富集显著性越可靠。
	                                |-- *_vs_*.Phase.tiff
				|-- BMK_3_GSEA	GSEA的分析结果，选择富集的前5个GO节点/pathway做富集图
					|-- BMK_1_Plot
						|-- *_vs_*_(BP/CC/MF)_GO_*.preranked.pdf
						|-- *_vs_*_(BP/CC/MF)_GO_*.preranked.png	GSEA preranked图,图中横坐标代表排序后的基因集的位置信息，纵坐标代表打分信息。每一条竖线代表感兴趣基因集合中的基因，竖线的高低代表打分的高低。
						|-- *_vs_*_(BP/CC/MF)_GO_*.runningScore.pdf
						|-- *_vs_*_(BP/CC/MF)_GO_*.runningScore.png	GSEA runningScore图,图中横坐标代表排序后的基因集的位置信息，纵坐标为动态的富集打分。图形顶部的竖线代表感兴趣基因集合中的基因。绿色的曲线代表每个位置上基因集合的富集打分。
						|-- *_vs_*_KEGG_pathway_ko*.preranked.pdf
						|-- *_vs_*_KEGG_pathway_ko*.preranked.png	同GO
						|-- *_vs_*_KEGG_pathway_ko*.runningScore.pdf
						|-- *_vs_*_KEGG_pathway_ko*.runningScore.png	同GO
					|-- *_vs_*_(BP/CC/MF)_GSEA.xls	GSEA分析结果,此文件已经按照qvalue列升序排序
						|-- ID	GO 节点ID
						|-- Description  节点描述   
						|-- setSize	基因个数
						|-- enrichmentScore 富集分数
						|-- NES     
						|-- pvalue  
						|-- qvalue  
						|-- rank
						|-- core_enrichment	基因ID
		|-- *_vs_*.all	该差异组合中所有基因表达量信息+Pvalue/FDR+log2FC+regulated信息(重要文件)
                |-- *_vs_*.DEG_cosmic.xls	癌症基因功能注释(重要文件)
			|-- ENSG_ID	基因ID
			|-- log2FC  差异log2FC
			|-- regulated	上下调(up/down)
			|-- Description	Cosmic中的描述信息
			|-- Tumour_Types(Somatic)	来自上代的突变
			|-- Tumour_Types(Germline)	体细胞突变
		|-- *_vs_*.DEG_TF.xls	转录因子注释(重要文件)
			|-- ENSG_ID
			|-- log2FC 
			|-- regulated
			|-- Family 转录因子家族
                |-- *_vs_*.DEG_final.xls	差异基因信息(重要文件)
                |-- *_vs_*.heatmap.pdf		差异mRNA 热图.pdf
                |-- *_vs_*.heatmap.png		差异mRNA 热图.png
                |-- *_vs_*_MA.pdf		差异miRNA MA图.pdf
                |-- *_vs_*_MA.png		差异miRNA MA图.png
                |-- *_vs_*_Volcano.pdf		差异mRNA 火山图.pdf
                |-- *_vs_*_Volcano.png		差异mRNA 火山图.pdf

	|-- BMK_PPI	差异基因蛋白互作结果
		|-- *_vs_*.DEG.detail.txt
			|-- #Query_id1	
			|-- Type    
			|-- Query_id2
			|-- Subject_id1     
			|-- Subject_id2     
			|-- Mode    
			|-- Score
		|-- *_vs_*.ppi.cytoscapeInput.sif	可直接导入cytoscape中的文件
	|-- DEG_anno.stat	差异基因注释数据库统计结果
	|-- DEG.stat	差异基因数目(up/down)统计结果
