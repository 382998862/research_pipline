.
|-- BMK_1_miRNA_Prediction
	|-- BMK_1_miR_Seq
		|-- All_miRNA.expressed.fa	预测得到的所有 miRNA 的 mature 序列(重要文件)
		|-- All_miRNA_Pre.expressed.fa	预测得到的所有 miRNA 的前体序列(重要文件)
		|-- Known_mature_expressed.fa	已知 miRNA 的 mature 序列
		|-- Known_Pre_expressed.fa	已知 miRNA 的前体序列
		|-- Novel_mature_expressed.fa	新预测的 miRNA 的 mature 序列
		|-- Novel_Pre_expressed.fa	新预测的 miRNA 的前体序列
	|-- BMK_2_Len_Distribution
		|-- All_miRNA.length.pdf
		|-- All_miRNA.length.png	所有 miRNA 长度分布图
		|-- All_miRNA.stat		所有 miRNA 长度分布统计文件
			|-- Length 长度
			|-- Number 此长度的miRNA的数目
		|-- Known_miRNA.length.pdf
		|-- Known_miRNA.length.png	已知 miRNA 长度分布图
		|-- Known_miRNA.stat		已知 miRNA 长度分布统计文件
		|-- Novel_miRNA.length.pdf
		|-- Novel_miRNA.length.png	新预测 miRNA 长度分布图
		|-- Novel_miRNA.stat		新预测 miRNA 长度分布统计文件
		|-- *.known_miRNA_len.stat	各样本已知 miRNA 长度分布统计文件
		|-- *.predicted_miRNA_len.stat	各样本新预测 miRNA 长度分布统计文件
	|-- BMK_3_Base_Distribution
		|-- All.base_bias.pdf
		|-- All.base_bias.png		miRNA 各位点碱基分布图,图中横坐标表示序列的各位点；纵坐标表示miRNA各位点碱基所占百分比
		|-- All.base_bias.stat		miRNA 各位点碱基分布统计文件,文件第一列为pos,之后每一列是该pos是分别是A/U/G/C的miRNA的数目
		|-- All.first_base.pdf
		|-- All.first_base.png		miRNA 首位碱基的分布图
		|-- All.first_base.stat		miRNA 首位碱基的统计文件
		|-- Known.base_bias.pdf
		|-- Known.base_bias.png		已知miRNA各位点碱基分布图
		|-- Known.base_bias.stat	已知miRNA各位点碱基分布统计文件
		|-- Known.first_base.pdf
		|-- Known.first_base.png	已知miRNA首位碱基的分布图
		|-- Known.first_base.stat	已知miRNA首位碱基的统计文件
		|-- Novel.base_bias.pdf
		|-- Novel.base_bias.png		新预测的miRNA各位点碱基分布图
		|-- Novel.base_bias.stat	新预测的miRNA各位点碱基分布统计文件
		|-- Novel.first_base.pdf
		|-- Novel.first_base.png	新预测的miRNA首位碱基的分布图
		|-- Novel.first_base.stat	新预测的miRNA首位碱基的统计文件
	|-- BMK_4_Summary_Stat	(重要文件夹)
		|-- Summary_Known_miRNA.txt	已知 miRNA 的总结文件,鉴定miRNA时,我们将比对上参考基因组的reads序列与已知miRNA数据库miRBase（v21）中的成熟miRNA序列进行比对。序列与已知miRNA完全相同的reads被认为是本项目中鉴定到的已知miRNA,因此此文件中的部分额列为"--"
			|-- miRNA	miRNA ID
			|-- score_total	"--"
			|-- mature_sequence	从miRbase数据库中提取的已知miRNA的成熟序列
			|-- star_sequence	"--"
			|-- GenomeID	"--"
			|-- strand	"--"
			|-- start	"--"
			|-- end	"--"
			|-- pre_seq	"--"
			|-- hairpin_stru	"--"
			|-- hairpin_energy	"--"
			|-- S*	miRNA在各个样品中的表达量信息
		|-- Summary_Novel_miRNA.txt	新预测 miRNA的总结文件,新miRNA用miRDeep2软件包进行预测
			|-- miRNA       	miRNA ID
                        |-- score_total		软件打分
                        |-- mature_sequence	预测得到的miRNA 的 mature 序列
                        |-- star_sequence       Star 序列
                        |-- GenomeID		序列所在染色体
                        |-- strand      	序列所在链的方向
                        |-- start       	前体序列的开始位置
                        |-- end			前体序列的结束位置
                        |-- pre_seq     	前体序列
                        |-- hairpin_stru        发夹结构
                        |-- hairpin_energy      发夹结构能量值
                        |-- S*  miRNA在各个样品中的表达量信息
		|-- Total_miRNA.stat	所有miRNA数据统计文件,统计了每个样品中预测得到的已知和Novel miRNA的数量及Total miRNA的数量
	|-- BMK_5_PDF	(重要文件夹)
		|-- *.pdf	二级结果预测
|-- BMK_2_miRNA_Expression
	|-- miRNA_counts.xls	所有 miRNA 的 count 数文件(重要文件)
	|-- miRNA_expression.xls	所有 miRNA 的表达量文件(重要文件)
	|-- PNG
		|-- all.fpkm_box.png	各样品的TPM箱线图
		|-- all.fpkm_density.png	各样品TPM密度分布对比图
		|-- *.fpkm_density.png	单个样品TPM密度分布图
		|-- sample_coefficient.txt	相关性系数表
		|-- sample_corelation.pdf	相关性系数图
		|-- sample_corelation.png
		|-- sample_pvalue.txt
|-- BMK_3_DEG_Analysis	差异分析结果
	|-- BMK_1_All_DEG	所有差异表达miRNA表达量和聚类分析,统计结果
		|-- All.DEG_final_anno.xls	所有差异miRNA表达量+差异信息+靶基因各个数据库注释结果(整合文件,信息全面,比较重要)
			|-- #Gene:差异miRNA靶向gene ID
			|-- ID:差异miRNA ID
                        |-- S*:样品表达量
                        |-- *_vs_*.Pvalue/FDR:差异Pvalue/FDR值,所有组合都有列出,若无,则用"--"代替
                        |-- *_vs_*.log2FC:差异log2FC值,所有组合都有列出,若无,则用"--"代替
                        |-- *_vs_*.regulated:差异上下调(up/down),若无,则用"--"代替
                        |-- 靶基因在COG/KOG/GO/KEGG/NR/eggNOG/Swissport各个数据库注释结果
		|-- All.DEG_final.pdf
		|-- All.DEG_final.png	所有差异miRNA聚类热图
		|-- All.DEG_final_target.xls	所有差异miRNA表达信息+差异信息+靶基因整合
		|-- All.DEG_final.xls	所有差异miRNA表达量+差异信息
			|-- #ID	差异miRNA ID
			|-- S*:样品表达量
			|-- *_vs_*.Pvalue/FDR:差异Pvalue/FDR值,所有组合都有列出,若无,则用"--"代替
			|-- *_vs_*.log2FC:差异log2FC值,所有组合都有列出,若无,则用"--"代替
			|-- *_vs_*.regulated:差异上下调(up/down),若无,则用"--"代替
		|-- Veen.pdf
		|-- Veen.png	Veen图(9>差异组合数目>=2)
	|-- BMK_*_*_vs_*
		|-- Anno_enrichment
			|-- anno	注释统计
				|-- *_vs_*.annotation.xls
				|-- *_vs_*.GO_enrichment.stat.xls
				|-- *_vs_*.GO.pdf
				|-- *_vs_*.GO.png
			|-- enrich	靶基因富集分析,采用R包clusterProfiler分别进行GO/KEGG富集分析,此目录下的BP/CC/MF分别是Biological_Process/Cellular_Component/Molecular_Function的简写
				|-- BMK_1_GO_enrich
					|-- *_vs_*_(BP/CC/MF)_enrich_barplot.pdf
					|-- *_vs_*_(BP/CC/MF)_enrich_barplot.png	差异表达基因GO富集条形
					|-- *_vs_*_(BP/CC/MF)_enrich_dotplot.pdf
					|-- *_vs_*_(BP/CC/MF)_enrich_dotplot.png	差异表达基因GO富集图(BP分支),横坐标为GeneRatio即注释在该条目中的感兴趣基因占所有差异lncRNA顺式靶基因数的比例,纵坐标每一个BP条目。点的大小代表该通路中注释的差异表达lncRNA顺式靶基因数，点的颜色代表超几何检验的校正后的p值。
					|-- *_vs_*_(BP/CC/MF)_enrich.list	富集结果文件
						|--ID	GO节点
	                                        |--Description	GO节点名称
        	                                |--GeneRatio	富集到此节点上的差异靶基因的数目/所有差异靶基因的数目
                	                        |--BgRatio	富集到此节点上的所有靶基因的数目/所有靶基因的数目
                        	                |--enrich_factor	富集因子
                                	        |--pvalue	Pvalue
                                        	|--qvalue	Qvalue
	                                        |--geneID	富集到此通路上的差异靶基因的ID
        	                                |--gene_symbol  对应的Name
				|-- BMK_2_KEGG_enrich
					|-- BMK_1_Keggmap    KEGG通路注释图
						|-- ko*.html
						|-- ko*.png 相对于对照组来说，红色框标记的酶与上调基因有关，绿色框标记的酶与下调基因有关。蓝色框标记的酶与上调和下调基因均有关，框内的数字代表酶的编号（EC number），而整个通路由多种酶催化的复杂生化反应构成，此通路图中与差异表达基因相关的酶均用不同的颜色标出，研究人员可以根据自己的研究对象间的差异，重点研究某些代谢通路相关基因的差异表达情况，通过通路解释表型差异的根源。
					|-- *_vs_*_KEGG_pathway_enrich_barplot.pdf
					|-- *_vs_*_KEGG_pathway_enrich_barplot.png	同GO
					|-- *_vs_*_KEGG_pathway_enrich_dotplot.pdf
					|-- *_vs_*_KEGG_pathway_enrich_dotplot.png	同GO
					|-- *_vs_*_KEGG_pathway_enrich.list	同GO
		|-- *_vs_*.all	该差异组合中所有基因表达量信息+Pvalue/FDR+log2FC+regulated信息(重要文件)
		|-- *_vs_*.DEG_final.Target.xls	该差异组合中所有基因表达量信息+Pvalue/FDR+log2FC+regulated+靶基因信息(重要文件)
		|-- *_vs_*.DEG_final.xls	差异基因信息(重要文件)
		|-- *_vs_*.heatmap.pdf	差异miRNA 热图.pdf
		|-- *_vs_*.heatmap.png	差异miRNA 热图.png
		|-- *_vs_*_MA.pdf	差异miRNA MA图.pdf
		|-- *_vs_*_MA.png	差异miRNA MA图.png
		|-- *_vs_*_Volcano.pdf	差异miRNA 火山图.pdf
		|-- *_vs_*_Volcano.png	差异miRNA 火山图.png
	|-- DEG_anno.stat	差异基因注释数据库统计结果
	|-- DEG.stat	差异基因数目(up/down)统计结果
|-- BMK_4_miRNA_Family	miRNA家族分析结果(重要文件夹)
	|-- Family_In_All.xls	miRBase 数据库中Chromalveolata|Metazoa|Mycetozoa|Viridiplantae|Viruses 某一大类下各物种 miRNA 在各家族中分类及数量
	|--第一列:表示家族名称
	|--第二列:表示物种在该家族中的数目
	|-- Family_In_miR.xls	所有 miRNA 的家族分类
	|--第一列:表示 miRNA 名称
	|--第二列:表示该 miRNA 的家族分类
|-- BMK_5_miRNA_Edit
	|-- miRNAEdit_cutoff_20.xls	count 数大于 20 作为阈值,过滤得到的碱基编辑文件
	|--第一列:表示 sRNA 和 pre-miRNA 的序列及匹配情况
	|--第二列:表示 round 数和碱基编辑类型,第一轮:完全匹配r0-r0;第二轮:允许一个错配,r1-M3,r1-M5,r1-MM;第三轮:第二轮过滤后的,去除 3 ‘端或5’端的一个碱基,允许一个错配, r2-M5, r2-M3
	|--第三列:表示 sRNA 的长度
	|--第四列:表示相对于前体的碱基变异情况
	|--第五列:表示 sRNA 在前体的位置
	|--第六列:到最后一列表示 sRNA 在各样本中的 count 数
|-- BMK_6_HMDD	microRNA人疾病相关数据库注释(重要文件夹)
	|-- HMDD_related_miRNA.txt	人疾病相关数据库注释结果
		|-- #Pre-miRNA	前体miRNA ID
		|-- related_miRNA	mature miRNA ID
		|-- Disease_name	疾病名称
		|-- PMID	可在NCBI中直接搜索
		|-- Description	描述
	|-- HMDD_related_miRNA_with_target_gene.txt	文件第三列增加了miRNA对应的靶基因的symbol
