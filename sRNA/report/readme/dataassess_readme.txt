测序数据评估及比对统计结果文件:
|-- BMK_1_Data_Assess	数据评估结果
	|-- miRNA.All_sample_filter.stat.xls	miRNA数据过滤统计表,过滤掉原始序列中含有接头的序列和低质量序列等一系列不符合要求的序列,得到符合要求的高质量序列
		|-- Samples:miRNA样本信息单样品名称
		|-- BMK-ID:百迈客编号
		|-- Raw_reads:测序原始数据
		|-- Low_quality:质量值低于 30 的碱基所占比例超过 20%的 reads
		|-- Containing'N'reads: 含至少 10%的未知碱基 N 的 reads
		|-- Length<18:去掉接头后,长度小于 18 个核苷酸的 reads 数
		|-- Length>30:去掉接头后,长度小于 18 个核苷酸的 reads 数
		|-- Clean_reads:质量值大于或等于 30 的碱基的 Reads 数
		|-- Q30(%):碱基质量值大于 30 的比例
	|-- miRNA.AllSample_GC_Q.stat.xls
		|-- SampleID:样本 ID
		|-- ReadSum:碱基数目
		|-- GC(%):GC 含量
		|-- N(%):含 N 的 reads 所占百分比
		|-- Q20(%):碱基质量值大于 20 的比例
		|-- Q30(%):碱基质量值大于 30 的比例
	|-- miRNA_PNG
		|-- *.error.png	miRNA样品测序错误率分布图,png格式;横坐标是碱基的位置,纵坐标是错误率,此图反映每个测序反应的测序质量,通常测序序列5'端前几个碱基的错误率相对较高
	|-- nonRib.AllSample_GC_Q.stat.xls	LncRNA样本的测序数据质量评估统计
		|-- SampleID:客户的样本编号
		|-- BMK-ID:百迈客样品分析编号,为公司内部编号,与客户样品的编号对应
		|-- ReadSum:Clean Data 中 pair-end Reads 总数
		|-- BaseSum:Clean Data 中的总碱基数
		|-- GC(%):Clean Data 中的 GC 含量,即Clean Data中 G 和 C 两种碱基占总碱基的百分比
		|-- N(%):Clean Data 中含 N 的碱基占总碱基的百分比
		|-- Q30(%):Clean Data 质量值大于或等于 30 的碱基所占的百分比
	|-- nonRib_PNG
		|-- *.acgtn.png	GC分布图
		|-- *.quality.png	碱基错误率分布图
|-- BMK_2_Mapped_Statistics	样本比对结果统计
	|-- circRNA	
		|-- All.mappedStat.xls	所有样本比对效率统计结果
			|-- BMK-ID：百迈客对样品的统一编号
			|-- Total Reads：Clean Reads数目，按单端计
			|-- Mapped Reads：比对到参考基因组上的Reads数目及在Clean Reads中占的百分比
			|-- Uniq Mapped Reads：比对到参考基因组唯一位置的Reads数目及在Clean Reads中占的百分比
			|-- Multiple Mapped Reads：比对到参考基因组多个位置的数目及在Clean Reads中占的百分比
			|-- Reads Map to '+'：比对到参考基因组正链的Reads数目及在Clean Reads中占的百分比
			|-- Reads Map to '-'：比对到参考基因组负链的Reads数目及在Clean Reads中占的百分比
		|-- *.type.png	样本测序 Reads 在不同区域分布饼图。图中将基因组分为外显子区、基因间区、内含子区,区域大小为 Map 到相应区域的 Reads 在所有 Mapped Reads 中所占的百分比。理论上,来自成熟 mRNA 的 Reads 应比对到外显子区。Reads 比对到内含子是由于 mRNA 前体和发生可变剪切的内含子保留;Reads比对到基因间区是由于基因组注释不完善。
   	|-- lncRNA
		|-- All.mappedStat.xls	同circRNA部分的统计
		|-- *.mappedStat.xls	单个样本测序数据与参考序列的比对统计表
			|-- Total Reads: Clean Reads 数目及其百分比(100%)
			|-- 比对到参考基因组上的 Reads 数目及在 Clean Reads 中占的百分比
			|-- 比对到参考基因组唯一位置的 Reads 数目及在 Clean Reads 中占的百分比
			|-- 比对到参考基因组多处位置的 Reads 数目及在 Clean Reads 中占的百分比	
			|-- 成对的(Paired)Reads 均比对到参考基因组上的 Reads 数目及其在 Clean Reads 中占的百分比
			|-- Single Map: 成对的(Paired)Reads中只有一条比对到参考基因组的 Reads 数目及在 Clean Reads 中占的百分比
			|-- Only Map Plus Strand:比对到参考基因组正链的 Reads 数目及在 Clean Reads 中占的百分比
			|-- Only Map Minus Strand 比对到参考基因组负链的 Reads 数目及在 Clean Reads 中 占的百分比
		|-- *.map.png	单个样本Mapped Reads在参考基因组上的位置及覆盖深度分布图。图中横坐标为染色体位置,纵坐标为覆盖深度以 2 为底的对数值,以 10kb 作为区间单位长度,划分染色体成多个小窗口(Window),统计落在各个窗口内的 Mapped Reads作为其覆盖深度。蓝色为正链,绿色为负链。
		|-- *.type.png	 同circRNA部分
	|-- miRNA
		|-- All_ncRNA_mapped.stat	sRNA所有样本分类注释统计结果
			|-- #ID:BMK-ID
			|-- Total:Clean reads 数目
			|-- rRNA:过滤掉的核糖体RNA的reads数目及比例
			|-- snRNA:细胞质小RNA
			|-- scRNA:核内小RNA
			|-- snoRNA:核仁小RNA
			|-- tRNA:转运RNA
			|-- Repbase:重复序列
			|-- Unannotated:包含miRNA的Unannotated reads
		|-- All_sample_map.stat	与参考基因组比对统计结果
			|-- BMK-ID：百迈客编号
			|-- Total_Reads：未被注释的用于与参考基因组比对的reads数目
			|-- Mapped_Reads：比对到参考基因组的Clean Reads
			|-- Mapped_reads(+)：比对正链上的Clean Reads数目
			|-- Mapped_reads(-)：比对到负链上的Clean Reads数目
		|-- *	
			|-- *.chro_distribution.png	样本的 Reads 在参考基因组上的位置及覆盖深度分布图。横坐标为染色体位置;纵坐标为染色体上对应位置的覆盖深度取以 2 为底的对数值
			|-- *.clean.fa	样品中 Clean Reads 的 fasta 文件
			|-- *.Clean_reads.length.png	样本 Clean Reads 长度分布统计图。横坐标为 small RNA 长度,纵坐标为特定长度的 small RNA 数目
			|-- *.Clean_reads.length.stat	样本 Clean Reads 长度分布统计表
				|-- length:表示 small RNA 长度
				|-- reads:表示该长度的 small RNA 数量
			|-- *.Genome_mapped_reads.length.png 样本比对到参考基因组上的 reads 长度分布统计图。横坐标为 small RNA 长度,纵坐标为特定长度的 small RNA 数目
			|-- *.Genome_mapped_reads.length.stat	样本比对到参考基因组上的 reads 长度分布统计表
				|-- length:表示 small RNA 长度
				|-- reads:表示该长度的 small RNA 数量
			|-- *.ncRNA_mapped.stat	单个样本分类注释统计结果
|-- BMK_3_Library_Assessment	文库质量评估
	|-- circRNA
		|-- *.insertSize.png	circRNA文库插入片段分布长度模拟分布图(png)
注:通过插入片段两端的 Reads 在参考基因组上的比对起止点之间的距离计算插入片段长度。图中横坐标为双端 Reads 在参考基因组上的比对起止点之间的距离,范围为 0 到 800bp;纵坐标为比对起止点之间不同距离的双端 Reads的数目。插入片段长度检验插入片段长度的离散程度能直接反映出文库制备过程中磁珠纯化的效果。
	|-- lncRNA
		|-- *.insertSize.png	同circRNA文库
		|-- *.randcheck.png	测序 Reads 在转录本上的位置分布图
注:横坐标为标准化后的 mRNA 位置,纵坐标为对应位置区间内 Reads 在总 Mapped Reads 中所占百分比。由于参考的 mRNA 长度不同,作图时把每个 mRNA 按照长度划分成 100 个区间,进而统计每一区间内的 Mapped Reads 数目及所占的比例,图中反映的是所有 mRNA 各个区间内的Mapped Reads 比例的汇总。
		|-- *.Saturation.png	测序数据的饱和度分布图
注:使用各样品的 Mapped Data 对检测到的不同表达情况的基因数目饱和情况进行模拟,绘制曲线图,可查看随着测序数据量的增加,检测到的不同表达量的基因数目是否趋于饱和。
图中是随机抽取 10%、20%、30%......90%的总体测序数据单独进行基因定量分析的结果;横坐标代表抽取数据定位到基因组上的 Reads 数占总定位的 reads 数的百分比,纵坐标代表所有抽样结果中表达量差距小于 15%的 Gene 在各个 FPKM范围的百分比。不同颜色代表不同的基因表达量范围。
		|-- Total.randcheck.png	所有样本测序 Reads 在转录本上的位置分布图
注:图中不同颜色的曲线代表不同的样品,横坐标为标准化后的 mRNA 位置,纵坐标为对应位置区间内 Reads 在总 Mapped Reads 中所占百分比。
