+--------------------------------------------------------------------------------------------------+
|                        有参考基因组转录组项目信息分析结果说明文档                                |
+--------------------------------------------------------------------------------------------------+

目录结构及文件说明：
****************************************************************************************************
BMKXXXXXX-XXX_Transcriptome_final/
|-- rawdata        #测序数据评估目录
|   |-- AllSample_GC_Q.stat  #测序数据评估统计表
|   `-- PNG                  #以Cycle为单位对测序数据进行评估作图的目录
|       |-- Sample1.acgtn.png          #样品1测序数据碱基含量分布图
|       |-- Sample1.cycleQ.png         #样品1测序数据平均质量值分布图
|       |-- Sample1.quality.png        #样品1测序数据碱基质量值分布图
|       |-- Raw_data_ratio_Sample1.png #样品1原始数据组成图
|       |-- ... ...
|       `-- Sample_data_relation.xls   #分析样品名与数据对应列表
|-- AllGene   全部基因的注释结果
|   |-- AllGene_Anno    基因注释部分
|       |-- All_Database_annotation.xls                 #功能注释整合文件,xls格式
|       |-- Function_Annotation.stat.xls                #注释数目统计表
|       |-- Species_Unigene.Cog_class.stat.xls          #COG注释分类统计表
|       |-- Species_Unigene.Cog_class.txt               #COG数据库注释结果文件
|       |-- Species_Unigene.Cog.cluster.png             #COG注释分类统计图（PNG格式）
|       |-- Species_Unigene.GO.anno.txt                 #GO数据库注释结果文件
|       |-- Species_Unigene.GO_enrichment.stat.xls      #GO二级节点注释统计表
|       |-- Species_Unigene.GO.list.txt                 #基因编号与注释到的GO节点编号对应表
|       |-- Species_Unigene.GO.pdf                      #GO二级节点注释统计图（PDF格式）
|       |-- Species_Unigene.GO.png                      #GO二级节点注释统计图（PNG格式）
|       |-- Species_Unigene.GO_tree.stat.xls            #GO节点与注释到该节点的基因编号对应表
|       |-- Species_Unigene.Kegg.globe                  #KEGG全局通路与注释到该全局通路的基因（Unigene）编号对应表
|       |-- Species_Unigene.Kegg.ko                     #KEGG数据库注释结果文件
|       |-- Species_Unigene.Kegg.pathway                #KEGG通路与注释到该通路的基因编号对应表
|       |-- Species_Unigene.nr.anno.evalue.png          #nr数据库BLAST比对E-value分布饼图（PNG格式）
|       |-- Species_Unigene.nr.anno.evalue.stat         #nr数据库BLAST比对E-value分布统计表
|       |-- Species_Unigene.nr.anno.evalue.svg          #nr数据库BLAST比对E-value分布饼图（SVG格式）
|       |-- Species_Unigene.nr.anno.identity.png        #nr数据库BLAST比对identity分布饼图（PNG格式）
|       |-- Species_Unigene.nr.anno.identity.stat       #nr数据库BLAST比对identity分布统计表
|       |-- Species_Unigene.nr.anno.identity.svg        #nr数据库BLAST比对identity分布饼图（SVG格式）
|       |-- Species_Unigene.nr.anno.txt                 #nr数据库注释结果文件
|       |-- Species_Unigene.nr.lib.png                  #nr数据库中注释到的物种分布饼图（PNG格式）
|       |-- Species_Unigene.nr.lib.stat                 #nr数据库中注释到的物种分布统计表
|       |-- Species_Unigene.nr.lib.svg                  #nr数据库中注释到的物种分布饼图（SVG格式）
|       |-- Species_Unigene.Swissprot.anno.txt          #Swiss-Prot数据库注释结果文件
|       |-- Integrated_Function.annotation.xls          #功能注释整合文件
|       |   `-- Kegg_map                                    #注释到的KEGG通路图目录
|       |       |-- ko00010.png
|       |-- All_gene.fa                                     #所有的基因序列
|       `-- All_Gene.longest_transcript.fa                  #所有基因序列的最长转录本序列
|-- NewGene
|       |-- Species.newGene_final.filtered.gff              #新基因的预测结果gff
|       |-- Species.newGene.longest_transcript.fa           #新基因最长转录本序列
|       |-- NewGene_Anno
|       |-- All_Database_annotation.xls                                           #功能注释整合文件,xls格式                                 
|       |-- Function_Annotation.stat.xls					      #注释数目统计表                                           
|       |-- Species.newGene.longest_transcript.fa.Cog_class.txt		      #COG数据库注释结果文件                                       
|       |-- Species.newGene.longest_transcript.fa.Cog.cluster.png		      #COG注释分类统计图（PNG格式）                                    
|       |-- Species.newGene.longest_transcript.fa.Cog.cluster.stat		      #COG注释分类统计                             
|       |-- Species.newGene.longest_transcript.fa.GO.anno.txt		      #GO数据库注释结果文件                                     
|       |-- Species.newGene.longest_transcript.fa.GO_enrichment.stat.xls	      #GO二级节点注释统计表                                     
|       |-- Species.newGene.longest_transcript.fa.GO.list.txt		      #新基因编号与注释到的GO节点编号对应表                       
|       |-- Species.newGene.longest_transcript.fa.GO.pdf			      #GO二级节点注释统计图（PDF格式）                          
|       |-- Species.newGene.longest_transcript.fa.GO.png			      #GO二级节点注释统计图（PNG格式）                          
|       |-- Species.newGene.longest_transcript.fa.GO_tree.stat.xls		      #GO节点与注释到该节点的新基因编号对应表                     
|       |-- Species.newGene.longest_transcript.fa.Kegg.globe		      #KEGG全局通路与注释到该全局通路的基因（Unigene）编号对应表
|       |-- Species.newGene.longest_transcript.fa.Kegg.ko			      #KEGG数据库注释结果文件                                   
|       |-- Species.newGene.longest_transcript.fa.Kegg.pathway		      #KEGG通路与注释到该通路的新基因编号对应表                   
|       |-- Species.newGene.longest_transcript.fa.nr.anno.txt		      #nr数据库注释结果文件              
|       |-- Species.newGene.longest_transcript.fa.nr.lib.png		      #nr数据库中注释到的物种分布饼图（PNG格式）                       
|       |-- Species.newGene.longest_transcript.fa.nr.lib.stat		      #nr数据库中注释到的物种分布统计表               
|       |-- Species.newGene.longest_transcript.fa.nr.lib.svg		      #nr数据库中注释到的物种分布饼图（SVG格式）             
|       |-- Species.newGene.longest_transcript.fa.Swissprot.anno.txt	      #Swiss-Prot数据库注释结果文件                                  
|       |-- Integrated_Function.annotation.xls				      #功能注释整合文件                               ）             
|       |   `-- Kegg_map							      #注释到的KEGG通路图目录                                        
|       |       |-- ko00062.png							                     
|       |       |-- ko00190.png							                     
|       `-- New_gene.fa								      #新基因基因序列               
|-- geneExpression #基因（Unigene）表达量分析结果目录
|   |-- Sample1.fpkm_density.png       #样品1的FPKM密度分布图
|   |-- Sample1.geneExpression.xls     #样品1基因表达量分析结果文件
|   |-- Sample1.insertSize.png         #样品1插入片段长度模拟分布图（png）
|   |-- Sample1.mappedStat.xls         #样品1测序数据与参考序列的比对统计表
|   |-- Saturation_Sample1.png         #样品1测序数据饱和度模拟图
|   |-- Sample1.coverage.png           #样品1Mapped Reads在参考基因组上的位置及覆盖深度分布图
|   |-- Sample1.randcheck.png          #样品1测序Reads在转录本上的位置分布图
|   |-- Sample1.type.png               #样品1测序Reads在不同区域分布饼图
|   |-- Sample1.isoExpression.xls      #样品1转录本表达量分析结果文件
|   |-- ... ...
|   |-- Total.randchecks.png            #各样品测序Reads在转录本上的位置分布整合图
|   |-- all.fpkm_box.png               #各样品FPKM箱线图
|   |-- all.fpkm_density.png           #各样品FPKM密度分布对比图
|   |-- AllSample.isoforms_expression.xls #各样品转录本表达整合表
|   |-- cor_plot                       #两两样品的表达量相关性散点图目录
|       |-- Sample1_vs_Sample2.cor.png     #样品1和样品2的表达量相关性散点图
|   |   `-- ... ...
|   |-- free_com.cor                   #两两样品的表达量相关性（皮尔逊相关系数的平方）统计表
|   `-- sample_cluster.png             #两两样品的表达量相关性热图
|-- DEG_Analysis   #差异表达分析结果目录
|   |-- All_gene_counts.list           #所有基因的表达量（比对片段数）矩阵文件
|   |-- All_gene_fpkm.list             #所有基因的表达量（fpkm）矩阵文件
|   |-- Condition1_vs_Condition2       #条件1样品和条件2样品的差异表达分析目录[1][2]
|       |-- Condition1_vs_Condition2.cor.png          #样品相关性点图
|       |-- Condition1_vs_Condition2.all              #未经log2fc和FDR过滤的差异分析结果
|       |-- Condition1_vs_Condition2.DEG.CytoscapeInput.txt       #差异基因蛋白互作分析结果
|       |-- Condition1_vs_Condition2.FC_count.png     #差异表达基因MA图
|       |-- Condition1_vs_Condition2.FC_FDR.png       #差异表达基因火山图
|       |-- Condition1_vs_Condition2.DEG.final.xls    #差异表达基因集，即差异表达分析结果文件
|       |-- Condition1_vs_Condition2.annotation.xls   #差异表达基因功能注释整合文件
|       |-- Condition1_vs_Condition2.nr.anno.txt      #差异表达基因nr库功能注释结果文件
|       |-- Condition1_vs_Condition2.Swissprot.anno.txt#差异表达基因Swiss-Prot库功能注释结果文件
|       |-- Cog_Anno                                  #差异表达基因COG注释目录
|       |   |-- Condition1_vs_Condition2.Cog_class.txt          #差异表达基因COG注释结果文件
|       |   |-- Condition1_vs_Condition2.Cog.classfy.png        #差异表达基因COG注释分类统计图
|       |   `-- Condition1_vs_Condition2.Cog_class.txt.stat     #差异表达基因COG注释分类统计表
|       |-- DEG_Cluster                               #差异表达基因表达模式聚类分析目录
|       |    |-- Condition1_vs_Condition2hclust_eudience_complete.txt    #差异表达基因表达模式聚类树型结果文件
|       |    |-- Condition1_vs_Condition2.pdf         #差异表达基因表达模式聚类图（pdf格式）
|       |    |-- Condition1_vs_Condition2.png         #差异表达基因表达模式聚类图（png格式）
|       |-- go_enrichment                             #差异表达分析GO注释富集分析结果目录
|       |   |-- Condition1_vs_Condition2.GO.anno.txt            #差异表达基因GO功能注释结果文件
|       |   |-- Condition1_vs_Condition2.GO.Biological.stat     #差异表达基因Fisher精确检验GO富集结果文件（生物学过程）
|       |   |-- Condition1_vs_Condition2.GO.Biological.xls      #差异表达基因Fisher精确检验GO富集统计文件（生物学过程）
|       |   |-- Condition1_vs_Condition2.GO.Cellular.stat       #差异表达基因Fisher精确检验GO富集结果文件（细胞组分）
|       |   |-- Condition1_vs_Condition2.GO.Cellular.xls        #差异表达基因Fisher精确检验GO富集统计文件（细胞组分）
|       |   |-- Condition1_vs_Condition2.GO.Molecular.stat      #差异表达基因Fisher精确检验GO富集结果文件（分子功能）
|       |   |-- Condition1_vs_Condition2.GO.Molecular.xls       #差异表达基因Fisher精确检验GO富集统计文件（分子功能）
|       |   |-- Condition1_vs_Condition2.GO_enrichment.stat.xls #差异表达基因GO二级节点注释统计表
|       |   |-- Condition1_vs_Condition2.GO.list.txt            #差异表达基因编号与注释到的GO节点编号对应表
|       |   |-- Condition1_vs_Condition2.GO.png                 #差异表达基因GO二级节点注释统计图（png格式）
|       |-- Graph                                     #差异表达基因功能富集统计作图目录
|       |   |-- Condition1_vs_Condition2.KEGG.list              #差异表达基因的KEGG富集结果文件
|       |   |-- Condition1_vs_Condition2.KEGG.Phase.png         #差异表达基因KEGG通路富集散点图
|       |   |-- Condition1_vs_Condition2.topGO_BP.pdf           #差异表达基因topGO富集有向无环图（pdf格式，生物学过程）
|       |   |-- Condition1_vs_Condition2.topGO_BP.png           #差异表达基因topGO富集有向无环图（png格式，生物学过程）
|       |   |-- Condition1_vs_Condition2.topGO_BP.xls           #差异表达基因topGO富集结果文件（生物学过程）
|       |   |-- Condition1_vs_Condition2.topGO_CC.pdf           #差异表达基因topGO富集有向无环图（pdf格式，细胞组分）
|       |   |-- Condition1_vs_Condition2.topGO_CC.png           #差异表达基因topGO富集有向无环图（png格式，细胞组分）
|       |   |-- Condition1_vs_Condition2.topGO_CC.xls           #差异表达基因topGO富集结果文件（细胞组分）
|       |   |-- Condition1_vs_Condition2.topGO_MF.pdf           #差异表达基因topGO富集有向无环图（pdf格式，分子功能）
|       |   |-- Condition1_vs_Condition2.topGO_MF.png           #差异表达基因topGO富集有向无环图（png格式，分子功能）
|       |   `-- Condition1_vs_Condition2.topGO_MF.xls           #差异表达基因topGO富集结果文件（分子功能）
|   |   `-- pathway                                   #差异表达基因KEGG功能富集分析结果目录
|   |       |-- kegg_enrichment                                 #差异表达基因KEGG功能富集分析结果目录
|   |       |   |-- Condition1_vs_Condition2.KEGG.png           #差异表达基因KEGG分类图（png格式）
|   |       |   |-- Condition1_vs_Condition2.KEGG.stat          #差异表达基因KEGG分类统计表
|   |       |   |-- Condition1_vs_Condition2.KEGG.svg           #差异表达基因KEGG分类图（svg格式）
|   |       |   |-- Condition1_vs_Condition2.KEGG.tree.stat     #差异表达基因KEGG分类统计详表
|   |       |   |-- Condition1_vs_Condition2.KEGG.xls           #KEGG通路注释结果统计文件
|   |       |   `-- Condition1_vs_Condition2.Kegg.ko            #差异表达基因KEGG注释结果文件
|   |       `-- kegg_map                                        #差异表达基因注释到的KEGG通路图目录
|   |           |-- ko00010.html                                     #差异表达基因注释到的KEGG通路图（html格式）
|   |           |-- ko00010.png                                      #差异表达基因注释到的KEGG通路图（png格式）
|   |           `-- ... ...
|   |-- ... ...
|   |-- All_DEG_veen.genes             #各差异表达基因集元素列表
|   |-- All_DEG_veen.png               #各差异表达基因集韦恩图
|   |-- All_DEG_veen.stat              #各差异表达基因集数目统计表
|   |-- DEG.anno.stat                  #各差异表达基因集的功能注释统计表
|   `-- DEG.stat                       #各差异表达基因集数目统计表
|-- Alt_splice     #可变剪接事件分析结果目录
|   |-- As_event_stat*.png             #各样品可变剪接事件数目统计图
|   |-- As_Result.stat                 #各样品预测的可变剪接事件结果统计
|   |-- Sample1.W9.fpkm                #样品1预测的可变剪接事件表格
|   `-- ... ...

|-- SNP_Analysis
|   |-- indel_anno
|   |   |-- all.anno.stat.png
|   |   |-- Dorsal.anno.stat.png
|   |   |-- final.indel.anno.gatk.Dorsal.list
|   |   |-- final.indel.anno.gatk.Spine.list
|   |   |-- final.indel.anno.gatk.Ventral.list
|   |   |-- final_Indel.anno.stat
|   |   |-- Spine.anno.stat.png
|   |   `-- Ventral.anno.stat.png
|   |-- Pairwised_SNP
|   |   |-- AlleDp.AlleDp.parwised_snp.hete_hete.list
|   |   |-- AlleDp.AlleDp.parwised_snp.hete_homo.list
|   |   |-- AlleDp.AlleDp.parwised_snp.homo_hete.list
|   |   |-- AlleDp.AlleDp.parwised_snp.homo_homo.list
|   |   |-- AlleDp.AlleDp.parwised_snp.list
|   |   |-- Alt.AlleDp.parwised_snp.hete_hete.list
|   |   |-- Alt.AlleDp.parwised_snp.hete_homo.list
|   |   |-- Alt.AlleDp.parwised_snp.homo_hete.list
|   |   |-- Alt.AlleDp.parwised_snp.homo_homo.list
|   |   |-- Alt.AlleDp.parwised_snp.list
|   |   `-- parwised_snp.stat
|   |-- snp_anno
|   |   |-- all.anno.stat.png
|   |   |-- Dorsal.anno.stat.png
|   |   |-- final.snp.anno.gatk.Dorsal.list
|   |   |-- final.snp.anno.gatk.Spine.list
|   |   |-- final.snp.anno.gatk.Ventral.list
|   |   |-- final_SNP.anno.stat
|   |   |-- Spine.anno.stat.png
|   |   `-- Ventral.anno.stat.png
|   |-- Spine.snp.type.png
|   `-- Ventral.snp.type.png


|-- SNP_Analysis   #SNP分析结果目录
|   |-- AllSample.SNP_density.png      #SNP密度分布图（png格式）
|   |-- AllSample.SNP_density.stat     #SNP密度分布统计表
|   |-- AllSample.snp.stat             #SNP数量统计表
|   |-- All.snp.type.png               #SNP突变类型分布图
|   |-- All.snp_type.stat              #SNP突变类型统计
|   |-- Sample1.snp.type.png           #样品1的SNP突变类型分布图
|   |-- Sample1.snp.stat.xls           #样品1的各染色体SNP位点统计表
|   |-- final.indel.anno.gatk.all.list #各样品inDel分析位点及注释信息统计表
|   |-- final.indel.anno.gatk.vcf      #各样品InDel分析位点及注释信息（vcf格式）
|   |-- final.snp.anno.gatk.all.list   #各样品SNP分析位点及注释信息统计表
|   |-- final.snp.anno.gatk.vcf        #各样品SNP分析位点及注释信息（vcf格式）
|   |-- indel_anno           #InDel注释目录
|   |   |-- all.anno.stat.png          #所有的InDel注释分类图
|   |   |-- final_Indel.anno.stat      #所有的InDel注释统计表
|   |   |-- Sample1.anno.stat.png      #样品1的InDel注释分类图
|   |   |-- final.indel.anno.gatk.Sample1.list #样品1的InDel注释统计表
|   |-- snp_anno           #InDel注释目录
|   |   |-- all.anno.stat.png          #所有的SNP注释分类图
|   |   |-- final_snp.anno.stat        #所有的SNP注释统计表
|   |   |-- Sample1.anno.stat.png      #样品1的SNP注释分类图
|   |   |-- final.snp.anno.gatk.Sample1.list #样品1的SNP注释统计表
|   |-- Pairwised_SNP
|   |   |-- AlleDp.AlleDp.parwised_snp.hete_hete.list
|   |   |-- AlleDp.AlleDp.parwised_snp.hete_homo.list
|   |   |-- AlleDp.AlleDp.parwised_snp.homo_hete.list
|   |   |-- AlleDp.AlleDp.parwised_snp.homo_homo.list
|   |   |-- AlleDp.AlleDp.parwised_snp.list
|   |   |-- Alt.AlleDp.parwised_snp.hete_hete.list
|   |   |-- Alt.AlleDp.parwised_snp.hete_homo.list
|   |   |-- Alt.AlleDp.parwised_snp.homo_hete.list
|   |   |-- Alt.AlleDp.parwised_snp.homo_homo.list
|   |   |-- Alt.AlleDp.parwised_snp.list
|   |   `-- parwised_snp.stat
|   |-- ... ...
|   |-- sam.merge.snp                  #各样品的SNP位点信息整合文件
|   |-- snp.total.stat                 #各样品SNP位点统计表
|   |-- AllSample.SNP_density.png      #各样品SNP位点密度分布图
|   `-- AllSample.SNP_density.stat     #各样品SNP位点密度分布统计表
|   |-- final.indel.vcf                #INDEL提取文件（vcf格式）
|   |-- final.snp.list                 #SNP文件，最主要的结果文件
|   |-- final.snp.vcf                  #SNP提取文件（vcf格式）
|-- Gene_Structure_Optimize  #基因结构优化分析目录
|   `-- Species.geneStructure.optimize.xls       #基因结构优化结果文件
`-- readme.txt     #说明文档

注：
[1] 对于设立生物学重复（Biological Replicates）的实验，同一条件会有多个样品（通常3个及以上），我们采用DESeq软件进行两条件样品组间的差异表达分析；而对于没有设立生物学重复的实验，同一条件只有1个样品，我们采用EBSeq软件进行两个样品间的差异表达分析。
[2] 在分析结果中，使用“A_vs_B”的方式命名差异表达基因集，如T1_vs_T2或T1_T2_vs_T3_T4等。通常情况下，对于两个样品之间的差异表达基因集，A表示对照样品、野生型样品或前一个时间点样品；而B表示对应的处理样品、突变型样品或后一个时间点样品。相应地，对于两个条件（即两组样品）之间的差异表达基因集，A表达含有多个重复样品（Duplicates）的对照组、野生型组或前一个时间点样品组；B表示对应的处理组、突变型组、后一个时间点样品组。根据两（组）样品之间表达水平的相对高低，差异表达基因可以划分为上调基因（Up-regulated Gene）和下调基因（Down-regulated Gene）。样品（组）B中的表达水平高于样品（组）A中的基因称之为上调基因；反之为下调基因。因此，上调和下调是相对的，由所给A和B的顺序决定，若更换A和B的顺序会完全反过来，但这不会对分析结果产生实质性的影响。



基本文件格式说明：
****************************************************************************************************
1. FASTQ格式：
===========================================================================
FASTQ格式是一种用于存放核酸序列及其对应的质量值的文本格式，是存储高通量测序数据的标准格式。我们的测序数据均以FASTQ格式存储，通常以“.fastq” 或 “.fq”为文件后缀（以“.fq.gz”为文件后缀的是压缩过的FASTQ格式文件）。格式如下：

@HWI-7001455:116:H8E71ADXX:2:2212:8408:20433 1:Y:0:ATCACG
CTTCAACCAGGTCACCGGCATCAACGTCATCAACTTCTACGCGCCGTTCATGTTCCGGACCATCGGGCTCAAGGAGAGCGCGTCCCTCATGTCGGCCGTGG
+
;<5;=<??#############################################################################################
@HWI-7001455:116:H8E71ADXX:2:2107:14756:56485 1:Y:0:ATCACG
CAGGTGCTGACAGCAATCGGTAACTTCAGCATCTGCTCCATTGGGGTCGGCATCCTCGTCGAGATCATCGTCATGTTCCCAATCCAGCACCGGAAGTACCG
+
;5;?#################################################################################################

FASTQ文件中通常每4行对应一个序列单元：
第一行 以@开头，后面接着序列标识（ID）以及其他可选的描述信息；
第二行 序列行；
第三行 以“+”开头，后面接着可选的描述信息；
第四行 第二行每个字母对应的质量值编码，长度必须和第二行的序列长度相同。


2. FASTA格式：
===========================================================================
FASTA格式是一种用于表示核苷酸序列或多肽序列的标准文本格式，其中核苷酸或氨基酸用单个字母表示，序列之前放置序列名字和注释信息。我们使用的参考基因组序列、组装得到的新基因或新转录本序列都以FASTA格式存储。FASTA格式文件通常以“.fasta” 或 “.fa”结尾。以CDS文件为例：

>c1000.graph_c0|orf1 type=3prime_partial len=462nt loc=c1000.graph_c0:62-523:-
TTGGACAGGTCTTGGATGGTTTTGGAGACTGCAGGGGGTCACGGCCGCAACCCCAATGGC
ATCCTCGGCCTTCTATGCCGGGAACACTTCCCTGGGCTTGTCGAGTACGCCGGAGTGACG
AGCCCAGCGTACACCTTCGACCACTACGCCGTCGCCCCCGATGCAGTAGATCGGGACGGC
AGACAATTCAACAACAAGGCGGAGCGGGTCAAGCAAGAGCTGTGGGTAAGTCTTCCTCGC
ACTATATTGTTTAATAAGTCGCATTTGTTGCATATTCTTGAAATAATGTATGGATACATC
GTCTTTGTATGTAGGATTTCTTCAGGTGCGATGCTGGATACGAGGCCAGGGCGGATGTGG
TGTCTACGACGTGCTGTAAGAAGCTCGTCGTGGACATGCACTACGAGGCGCGCATCCAGG
CCATCGTCACTTACCACGGCTCCGTCCTTGGGGAGAAGGTGA
>c1020.graph_c0|orf1 type=complete len=705nt loc=c1020.graph_c0:550-1254:+
ATGCTGAATCTAAGTAGCTCTGCTGTATCTGCGCCATCAAGAACTCATGTGGATCATGCC
GCTCTTACCGGTGCCTCTCATCCAGCTTCTACGGTCAAGACACGCATGTGCACCAAGTAC
AACACTACAGAAGGCTGCAAGTTCGGTGATAAGTGCCATTTCGCTCATAGCGAAAGAGAG
CTTGCGAAGCCAGCCTACATGTCTCAAGAAGGACCTCCTATGGGTGGTCGATATGGACGA
GCTGAACCTATGCAACCAGCTGCCATGGGCCCTCCAGCAGGAAACTTTGGTGCCTCGGCG
ACTGCCAAGATCAGTGTGGACGCCTCTCTGGCCGGTGGCATAATCGGCAAGGGTGGGGTC
AACACTAAGCAGATAAGTAGAATTACAGGCGTCAAGCTCTCCATCCGCGACCACGAGTCT
AACCCCAACCTAAAGAACATCGAGCTGGAAGGCAATTTTGACCAGATCAAGCAAGCCAGC
GACTTGGTGCGTGATCTCATCGCAACCATCAGCGCAAGCATGCCAGCGAAAAATCCATCT
GCTGCCGCGGCACCAGCAGGAGGAGGCCGAGGTGGTGGTCCAGGGGGCAAGAGCAGCAAC
TACAAGACGAAGCTGTGCGAGAACTTCTTGAAGGGTGCCTGTACTTTTGGTGACCGGTGC
CACTTCGCCCATGGCGAGACGGAGCAGCGGAGAGGTGCTGCATGA

FASTA格式每一个序列单元由描述信息和序列数据两部分构成，每一个序列单元从“>”（大于号）开始到下一个“>”之前结束。
描述信息与“>”放于第一行，“>”后面紧接序列标识（ID，如c1000.graph_c0|orf1），二者之间不能有空格，之后是可选的描述；
序列数据从第二行开始，可放于单行或多行，通常每一行为60或100个字符，直到出现“>”之前为止。
“>”标志着一个新的序列单元的开始。

3. GFF格式：
===========================================================================
通用特征格式（General Feature Format、Generic Feature Format，GFF）是桑格研究所定义的用于描述基因和其它DNA、RNA或蛋白质序列特征的一种简单、方便的文本格式，又称为基因发现格式（Gene-Finding Format）。目前已经成为序列注释的通用格式，比如基因组的特征注释，许多软件都支持输入或者输出GFF格式，比如GBrowse、IGB等等。我们分析结果中新基因（Unigene）基因组位置注释信息等都以GFF格式存储，文件一般以“gff”、“gtf”或“gff3”等作为扩展名。格式示意如下：

##gff-version 3
chr01	BMK	gene	1000	9000	.	+	.	ID=gene00001;Name=DEMO
chr01	BMK	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=DEMO.1
chr01	BMK	exon	1050	1500	.	+	.	Parent=mRNA00001
chr01	BMK	five_prime_UTR	1050	1200	.	+	.	Parent=mRNA00001
chr01	BMK	CDS	1201	1500	.	+	0	Parent=mRNA00001;Name=demoprotein.1
chr01	BMK	exon	3000	3902	.	+	.	Parent=mRNA00001
chr01	BMK	CDS	3000	3902	.	+	0	Parent=mRNA00001;Name=demoprotein.1
chr01	BMK	exon	5000	5500	.	+	.	Parent=mRNA00001
chr01	BMK	CDS	5000	5500	.	+	0	Parent=mRNA00001;Name=demoprotein.1
chr01	BMK	exon	7000	9000	.	+	.	Parent=mRNA00001
chr01	BMK	CDS	7000	7600	.	+	0	Parent=mRNA00001;Name=demoprotein.1
chr01	BMK	three_prime_UTR	7601	9000	.	+	0	Parent=mRNA00001;Name=demoprotein.1

GFF文本文件由表头行和特征行构成。表头行可选，一般以“#”开头，用于存放文件描述信息，放于文件开头。
特征行每行共9列，由TAB键分割。特征行各列说明如下表：
第一列：用于建立该特征坐标体系的序列标识；
第二列：产生特征的算法或数据库来源，比如”GeneScan”、”GenBank” 等，缺失值用“.”表示；
第三列：特征的类型，如Gene、cDNA、mRNA等等；
第四列：（基于1的坐标体系下）特征的开始位置，开始位置一般要小于或等于结束位置；
第五列：（基于1的坐标体系下）特征的结束位置，结束位置不能大于序列的长度；
第六列：特征打分，浮点值，是特征可能性的说明，可以是序列比对时的E值或者基因预测时的P值。缺失值为“.”；
第七列：特征所在的链。“+”表示正链，“-”表示负链，“?”表示未知；
第八列：仅对“CDS”特征类型有效，表示起始编码的位置，有效值为0、1、2，缺失值为“.”；
第九列：以多对标签-值构成的特征属性描述。

GTF格式是GFF2的改进版本，有时候也被当做GFF2.5。


----------------------------------------------------------------------------------------------------
tree v1.5.0 1996-2009 by Steve Baker
copyright (c) BMK 2014 by Simon Young

