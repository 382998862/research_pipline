.
|-- BMK_1_Global	样品表达量全局展示结果
	|-- BMK_1_All
		|-- Total.pie.png	四种RNA的数量统计图
		|-- Total.stat.xls
		|-- W*	每个样本预测得到的所有RNA表达量全局展示结果
			|-- circRNA.txt	W*样品中circRNA的位置信息(Chr Start End Count)
			|-- lncRNA.txt	W*样品中lncRNA的位置信息(Chr Start End Count)
			|-- miRNA.txt	W*样品中miRNA前体的位置信息(Chr Start End Count)
			|-- mRNA.txt	W*样品中mRNA的位置信息(Chr Start End Count)
			|-- W*.barplot.txt	exp_barplot,作图文件
			|-- W*.circos.png	表达量全局展示Circos图，图中最外圈为染色体信息，之后依次是mRNA、lncRNA、circRNA、miRNA的表达水平。这里，表达量均取对数，miRNA位置采用其前体RNA的位置信息。
			|-- W*.circos.svg
			|-- W*.exp_barplot.pdf
			|-- W*.exp_barplot.png	RNA在不同染色体的表达分布，图中不同颜色代表不同的RNA类型，横坐标代表不同染色体，纵坐标代表RNA在该染色体上表达量的比例。同一种RNA在不同染色体上的总分布和为1。这里，表达量均取对数。
	|-- BMK_2_DEG	差异RNA表达量全局展示结果
		|-- *_vs_*
			|-- circRNA.txt	W*样品中差异circRNA的位置信息(Chr Start End Count)
			|-- gene.txt	W*样品中差异lncRNA的位置信息(Chr Start End Count)
			|-- lncRNA.txt	W*样品中差异miRNA前体的位置信息(Chr Start End Count)
			|-- miRNA.txt	W*样品中差异mRNA的位置信息(Chr Start End Count)
			|-- *_vs_*.Circos.png	表达量全局展示,Circos图,图中最外圈为染色体信息，之后依次是mRNA、lncRNA、circRNA、miRNA。每个差异分组的圈图中红色代表上调，蓝色代表下调，高度代表显著性（-log10(FDR)）
			|-- *_vs_*.Circos.svg
			|-- *_vs_*.pie.png	差异表达的RNA数目分布图
			|-- *_vs_*.stat.xls	差异RNA数目统计
|-- BMK_2_conservation	保守性分析结果。LncRNA 的序列保守性相对 mRNA 要低，使用phastCons软件 (http://compgen.bscb.cornell.edu/phast/) 分别对 mRNA 和lncRNA 以及circRNA进行保守性打分，得到 lncRNA 和 mRNA 以及circRNA的保守性分值累积分布图
	|-- All.phastcons.boxplot.pdf
	|-- All.phastcons.boxplot.png	不同类型RNA保守性分值柱状图
	|-- All.phastcons.cdf.pdf
	|-- All.phastcons.cdf.png	不同类型RNA保守性分值累积分布图
	|-- All.phastcons.density.pdf
	|-- All.phastcons.density.png	不同类型RNA保守性分值概率密度曲线
	|-- all_type.score	所有RNA分值统计文件，第一列是分值，第二列是RNA类型
	|-- circRNA.phastcons.cdf.pdf
	|-- circRNA.phastcons.cdf.png	circRNA保守性分值累积分布图
	|-- circRNA.PhastCons.score
	|-- gene.phastcons.cdf.pdf
	|-- gene.phastcons.cdf.png	gene保守性分值累积分布图
	|-- gene.PhastCons.score
	|-- lncRNA.phastcons.cdf.pdf
	|-- lncRNA.phastcons.cdf.png	lncRNA保守性分值累积分布图
	|-- lncRNA.PhastCons.score
|-- BMK_3_ceRNA	ceRNA分析结果
	|-- ceRNA_pair_adjust_p_Sig_diff_node.xls	差异ceRNA关系对对应的节点文件
	|-- ceRNA_pair_adjust_p_Sig_diff.xls	差异ceRNA关系对
	|-- ceRNA_pair_adjust_p_Sig.xls	ceRNA关系对
	|-- Coexpression_ceRNA_pair_adjust_p_Sig_diff_node.xls	若样本数目不少于5个，则将共表达的结果考虑到ceRNA中，即共表达差异ceRNA关系对对应的节点文件
	|-- Coexpression_ceRNA_pair_adjust_p_Sig_diff.xls	若样本数目不少于5个，则将共表达的结果考虑到ceRNA中，即共表达差异ceRNA关系对
	|-- Coexpression_ceRNA_pair_adjust_p_Sig.xls	若样本数目不少于5个，则将共表达的结果考虑到ceRNA中，即ceRNA关系对
	|-- random	关键基因通路整合分析
		|-- KEGG
			|-- Key_KEGG_pathway_enrich_barplot.pdf
			|-- Key_KEGG_pathway_enrich_barplot.png
			|-- Key_KEGG_pathway_enrich_dotplot.pdf
			|-- Key_KEGG_pathway_enrich_dotplot.png
			|-- Key_KEGG_pathway_enrich.list
			|-- network.pdf
			|-- network.png
			|-- network.txt
			|-- node.txt
		|-- Key_gene.txt
		|-- Key_RNA_top0.10.txt
		|-- Net_score.txt
|-- BMK_4_cytoscape	lncRNA/gene/circRNA/miRNA靶向关系结果

########
说明：此部分结果是给予差异RNA进行的分析，差异RNA的数量会影响结果的有无。若发现结果为空或者数量极少，则建议老师根据自己感兴趣的RNA个性化分析这部分结果。
重要：此部分结果是分别以lncRNA/circRNA/miRNA/mRNA为中心，整合与其有靶向关系的其他RNA，并构建Cytoscape输入文件，利用此部分结果可以绘制以差异RNA为主体的精美的网络图
########
文件说明：m代表mRNA，l代表lncRNA，c代表circRNA，s代表miRNA，绘制网络图会用到此部分结果，建议筛选出关键RNA之后再绘制。
*Attribution.list是cytoscape输入属性文件
*Interaction.list是cytoscape输入互作文件
*.degset.xls	每个Class包含的RNA的数目和ID
*.DEG.xls	以差异RNA为中心，整合与其有靶向关系的其他差异RNA的信息
*.venn.png	Venn图
*.vennset.xls	Venn图每个部分包含的RNA的数目和ID
*.ceRNA_pair_adjust_p_Sig_diff.xls	基于ceRNA的分析结果，提取三者均差异的ceRNA关系对
########

	|-- circRNA
		|-- *_vs_*.cm.Attribution.list
		|-- *_vs_*.cm.Interaction.list
		|-- *_vs_*.cs.Attribution.list
		|-- *_vs_*.cs.Interaction.list
		|-- *_vs_*.degset.xls	每个Class包含的RNA的数目和ID
		|-- *_vs_*.DEG.xls	以差异circRNA为中心进行靶向关系联合分析，将差异circRNA及与其Host gene和与其具有靶向关系的差异miRNA，整理成表格
		|-- *_vs_*.venn.png	Venn图
		|-- *_vs_*.vennset.xls	Venn图每个部分包含的RNA的数目和ID
	|-- circRNA-miRNA-mRNA
		|-- *_vs_*.circRNA-miRNA-mRNA.Attribution.list
		|-- *_vs_*.circRNA-miRNA-mRNA.ceRNA_pair_adjust_p_Sig_diff.xls
		|-- *_vs_*.circRNA-miRNA-mRNA.Interaction.list
	|-- gene
		|-- *_vs_*.Cis.degset.xls
		|-- *_vs_*.Cis.venn.png
		|-- *_vs_*.Cis.vennset.xls
		|-- *_vs_*.DEG.xls
		|-- *_vs_*.mc.Attribution.list
		|-- *_vs_*.mc.Interaction.list
		|-- *_vs_*.ml.Cis.Attribution.list
		|-- *_vs_*.ml.Cis.Interaction.list
		|-- *_vs_*.ml.Trans.Attribution.list
		|-- *_vs_*.ml.Trans.Interaction.list
		|-- *_vs_*.ms.Attribution.list
		|-- *_vs_*.ms.Interaction.list
		|-- *_vs_*.Trans.degset.xls
		|-- *_vs_*.Trans.venn.png
		|-- *_vs_*.Trans.vennset.xls
	|-- lncRNA
		|-- *_vs_*.Cis.degset.xls
		|-- *_vs_*.Cis.venn.png
		|-- *_vs_*.Cis.vennset.xls
		|-- *_vs_*.DEG.xls
		|-- *_vs_*.lm.Cis.Attribution.list
		|-- *_vs_*.lm.Cis.Interaction.list
		|-- *_vs_*.lm.Trans.Attribution.list
		|-- *_vs_*.lm.Trans.Interaction.list
		|-- *_vs_*.ls.Attribution.list
		|-- *_vs_*.ls.Interaction.list
		|-- *_vs_*.Trans.degset.xls
		|-- *_vs_*.Trans.venn.png
		|-- *_vs_*.Trans.vennset.xls
	|-- lncRNA-miRNA-mRNA
		|-- *_vs_*.lncRNA-miRNA-mRNA.Attribution.list
		|-- *_vs_*.lncRNA-miRNA-mRNA.ceRNA_pair_adjust_p_Sig_diff.xls
		|-- *_vs_*.lncRNA-miRNA-mRNA.Interaction.list
	|-- sRNA
		|-- *_vs_*.degset.xls
		|-- *_vs_*.DEG.xls
		|-- *_vs_*.sc.Attribution.list
		|-- *_vs_*.sc.Interaction.list
		|-- *_vs_*.sl.Attribution.list
		|-- *_vs_*.sl.Interaction.list
		|-- *_vs_*.sm.Attribution.list
		|-- *_vs_*.sm.Interaction.list
		|-- *_vs_*.venn.png
		|-- *_vs_*.vennset.xls
|-- BMK_5_coexp	共表达分析结果

########
说明：此目录下的每个文件均为四列，RNA1/RNA2分别是用于共表达分析的RNA,coefficient为共表达相关性系数，pvalue为显著性Pvalue)
########
	|-- Diff.Sig.circRNA_miRNA.xls	circRNA和差异miRNA共表达关系对，circRNA和miRNA至少一个是差异RNA
	|-- Diff.Sig.lncRNA_mRNA.xls	lncRNA和差异miRNA共表达关系对，lncRNA和miRNA至少一个是差异RNA
	|-- Sig.circRNA_miRNA.xls	circRNA和miRNA共表达关系对
	|-- Sig.lncRNA_mRNA.xls		lncRNA和miRNA共表达关系对
