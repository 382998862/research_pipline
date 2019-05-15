.
|-- BMK_1_circRNA_Prediction
	|-- BMK_1_CIRC_info	(此部分结果是circRNA预测的结果,比较重要)
		|-- BMK_1_CIRI	CIRI预测结果(如果使用了此款软件进行circRNA预测)
			|-- *.CIRI_Predict.xls	CIRI 软件预测环状 RNA 结果
		|-- BMK_2_find_circ	find_circ软件预测结果(如果使用了此款软件进行circRNA预测)
			|-- *.find_circ_Predict.xls	find_circ 软件预测环状 RNA 结果
				|-- chrom：染色体名称
				|-- start：左剪切位点（以0开始）
				|-- end：右剪切位点（以0开始）
				|-- name：环状RNA名称
				|-- n_reads：number of reads supporting the junction (BED'score')
				|-- strand：+/-；n_uniq：支持junction的唯一reads数目
				|-- uniq_bridges：锚点比对唯一的reads数目
				|-- best_qual_left：支持左剪切junction最好的锚点比对分数
				|-- best_qual_right：支持右剪切junction最好的锚点比对分数
				|-- tissues：样品名称，以逗号分割
				|-- tiss_counts：每个样品对应的reads数，以逗号分割
				|-- edits：锚点延伸过程的错配数
				|-- anchor_overlap：断点在一个锚点内的核苷酸数目
				|-- breakpoints：以侧翼GT/AG方式打断reads的数目
				|-- signal：侧翼核苷酸剪切信号(GT/AG)
				|-- strandmatch：'MATCH','MISMATCH' or 'NA'（非链特异性分析）
				|-- category：描述junction的关键字
		|-- *.circ.intersect.xls	CIRI和find_circ取交集后整理的结果（如果同时使用这两款软件进行circRNA预测；否则是经过整理后的某一款软件的结果）
			|-- chr		circRNA所在染色体
			|-- start	circRNA开始位置
			|-- end		circRNA结束位置
			|-- Mean_junction	平均junction reads数目,即同一个circRNA,求两款软件的平均junction reads数目
			|-- circRNA_type
			|-- gene_id	circRNA 的来源基因ID
			|-- strand	circRNA 所在链
	|-- BMK_2_statistics
		|-- all.chromosome_distrbution.pdf
		|-- all.chromosome_distrbution.png	所有环状 RNA 在染色体上的分布图,图中横坐标表示染色体，纵坐标表示对应染色体上的环状RNA junction reads数量
		|-- All_circRNA.circRNA.length.distribution.pdf
		|-- All_circRNA.circRNA.length.distribution.png	环状 RNA 长度分布图,图中横坐标代表circRNA的长度区间；纵坐标代表该长度区间内的环状RNA数目。
		|-- circRNA_circos.png	circos 图(每个样品表达量以及环状 RNA 分布).
					注：CircRNA分布在基因组的统计以及各个样品表达量的circos图。其中，柱状图表示各个位置的CircRNA的个数，热图表示每个样品的表达量分布，0~1：蓝色，1~10：绿色，10~100：橘色，100～：红色，热图由外向内依次表示T01、T02、T03等。
		|-- circRNA_circos.svg
		|-- *.chromosome_distrbution.pdf
		|-- *.chromosome_distrbution.png	每个样品中环状RNA在染色体上的分布图
		|-- *.circRNA.type.png	环状 RNA 来源分布图,即circRNA来源于基因组上的位置(intron/exon/intergenic_region)
	|-- BMK_3_Known	(如果老师关注预测得到的circRNA是否是已知circRNA,可以查看这部分结果)
		|-- circRNA_known.xls	已知circRNA分析结果
			|-- 预测得到的circRNA ID
			|-- 对应的已知circRNA 的ID
			|-- blast identity 值
		|-- Known_Unknown_pie.png	已知circRNA与新预测circRNA数目统计饼图
	|-- circRNA.fa	预测的环状 RNA 的序列文件(重要文件)
	|-- circRNA_newname.xls	环状 RNA 重新命名的文件
		|-- New_ID:circRNA 新名字,命名规则是物种简写_来源基因Symbol_排序号
		|-- circRNA_ID:预测ID
		|-- L*:junction reads数目
		|-- geneLength:circRNA长度,这里的长度是只成环之后的长度,而不是直接用预测的Start-End
		|-- gene_id:来源基因ID
		|-- gene_symbol:来源基因symbol
	|-- overlap_alitisplice.xls	可变剪切文件
|-- BMK_2_circRNA_Expression
		|-- circRNA_counts.xls	环状 RNA junction reads count结果文件(重要文件)
		|-- circRNA_expression.xls	环状 RNA 表达量结果文件(重要文件)
		|-- PNG
			|-- all.fpkm_box.png	各样品的SRPBM/TPM箱线图
			|-- all.fpkm_density.png	 各样品SRPBM密度分布对比图
			|-- *.fpkm_density.png	每个样品SRPBM密度分布图
			|-- sample_coefficient.txt	相关性系数表
			|-- sample_corelation.pdf
			|-- sample_corelation.png	样品相关性关系图
			|-- sample_pvalue.txt
|-- BMK_3_DEG_Analysis
	|-- BMK_1_All_DEG	所有差异circRNA统计结果
		|-- All.DEG_final_anno.xls	所有差异circRNA表达量+差异信息+靶基因各个数据库注释结果(所有差异circRNA整合文件,信息较全面,方便查看,比较重要)
			|-- #Gene:差异circRNA靶向gene ID
                        |-- ID:差异circRNA ID
                        |-- L*:样品表达量
                        |-- *_vs_*.Pvalue/FDR:差异Pvalue/FDR值,所有组合都有列出,若无,则用"--"代替
                        |-- *_vs_*.log2FC:差异log2FC值,所有组合都有列出,若无,则用"--"代替
                        |-- *_vs_*.regulated:差异上下调(up/down),若无,则用"--"代替
                        |-- 靶基因在COG/KOG/GO/KEGG/NR/eggNOG/Swissport各个数据库注释结果
		|-- All.DEG_final.pdf
		|-- All.DEG_final.png	聚类热图
		|-- All.DEG_final_target.xls	所有差异circRNA表达信息+差异信息+靶基因整合
		|-- All.DEG_final.xls	所有差异miRNA表达量+差异信息
			|-- #ID 差异miRNA ID
                        |-- L*:样品表达量
                        |-- *_vs_*.Pvalue/FDR:差异Pvalue/FDR值,所有组合都有列出,若无,则用"--"代替
                        |-- *_vs_*.log2FC:差异log2FC值,所有组合都有列出,若无,则用"--"代替
                        |-- *_vs_*.regulated:差异上下调(up/down),若无,则用"--"代替
		|-- Veen.pdf
		|-- Veen.png	Veen图(9>差异组合数目>=2)
	|-- BMK_*_*_vs_*
		|-- Anno_enrichment
			|-- anno
			|-- *_vs_*.annotation.xls	差异基因Count+FPKM+Pvalue/FDR+log2FC+regulated+各个数据库注释结果
			|-- *_vs_*.GO_enrichment.stat.xls	差异基因GO二级节点注释统计
			|-- *_vs_*.GO.pdf
			|-- *_vs_*.GO.png
		|-- enrich	(重要文件夹)差异表达circRNA靶基因GO富集分析/KEGG富集分析分析.此部分分析是采用R包clusterProfiler对circRNA靶基因分别进行生物学过程，分子功能和细胞组分的富集分析以及KEGG通路富集分析。富集分析采用超几何检验方法来寻找与整个基因组背景相比显著富集的GO条目/KEGG通路。此目录下的BP/CC/MF分别是Biological_Process/Cellular_Component/Molecular_Function的简写。
			|-- *_vs_*_(BP/CC/MF)_enrich_barplot.pdf
			|-- *_vs_*_(BP/CC/MF)_enrich_barplot.png	差异表达circRNA靶基因GO富集条形图
			|-- *_vs_*_(BP/CC/MF)_enrich_dotplot.pdf
			|-- *_vs_*_(BP/CC/MF)_enrich_dotplot.png	差异表达circRNA靶基因GO富集图,横坐标为GeneRatio即注释在该条目中的感兴趣基因占所有差异circRNA靶基因数的比例,纵坐标每一个BP/CC/MF条目。点的大小代表该通路中注释的差异表达circRNA靶基因数，点的颜色代表超几何检验的校正后的p值。
			|-- *_vs_*_(BP/CC/MF)_enrich.list	富集结果文件
				|--ID:GO节点
                                |--Description:GO节点名称
                                |--GeneRatio:富集到此节点上的差异靶基因的数目/所有差异靶基因的数目
                                |--BgRatio:富集到此节点上的所有靶基因的数目/所有靶基因的数目
                                |--enrich_factor:富集因子
                                |--pvalue:Pvalue
                                |--qvalue:Qvalue
                                |--geneID:富集到此通路上的差异靶基因的ID
                                |--gene_symbol:对应的Name
			|-- *_vs_*_KEGG_pathway_enrich_barplot.pdf
			|-- *_vs_*_KEGG_pathway_enrich_barplot.png	差异表达circRNA靶基因KEGG富集条形图,注：横坐标代表差异表达circRNA靶基因注释在通路中的基因数，纵坐标代表通路，柱的颜色代表校正后的p值。
			|-- *_vs_*_KEGG_pathway_enrich_dotplot.pdf
			|-- *_vs_*_KEGG_pathway_enrich_dotplot.png	差异表达lncRNA顺式靶基因KEGG富集图
		 	|-- *_vs_*_KEGG_pathway_enrich.list	同BP
		|-- pathway
			|-- kegg_map
		        	|-- *.html
		        	|-- *.png 相对于对照组来说，红色框标记的酶与上调基因有关，绿色框标记的酶与下调基因有关。蓝色框标记的酶与上调和下调基因均有关，框内的数字代表酶的编号（EC number），而整个通路由多种酶催化的复杂生化反应构成，此通路图中与差异表达基因相关的酶均用不同的颜色标出，研究人员可以根据自己的研究对象间的差异，重点研究某些代谢通路相关基因的差异表达情况，通过通路解释表型差异的根源。
	|-- *_vs_*.all	该组合中所有LncRNA的表达量/log2FC/Pvalue/FDR(重要文件)
	|-- *_vs_*.DEG_final.Target.xls	差异circRNA文件增加靶基因ID(重要文件)
	|-- *_vs_*.DEG_final.xls	差异circRNA文件(重要文件)
	|-- *_vs_*.heatmap.pdf	差异circRNA 热图.pdf
	|-- *_vs_*.heatmap.png	差异circRNA 热图.png
	|-- *_vs_*_MA.pdf	差异circRNA MA图.pdf
	|-- *_vs_*_MA.png	差异circRNA MA图.png
	|-- *_vs_*_Volcano.pdf	差异circRNA 火山图.pdf
	|-- *_vs_*_Volcano.png	差异circRNA 火山图.png
	|-- DEG_anno.stat	各组差异circRNA上下调数目统计
	|-- DEG.stat		各组差异circRNA在不同数据库注释数目统计
