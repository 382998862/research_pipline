use Getopt::Long;
use Getopt::Std;
use Config::General;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);

my ($in,$od,$data_cfg,$detail_cfg);
GetOptions(
	"h|?"   	        =>\&USAGE,
	"od:s"     		=>\$od,
	"in:s"			=>\$in,
	"step:s"		=>\$step,
	"cfg1:s"		=>\$data_cfg,
	"cfg2:s"		=>\$detail_cfg,
)or &USAGE;
&USAGE unless ($in);


$step||= join(",",1..9);
my @steps=split(/,/,$step);
my %steps=();
foreach my $s(@steps){$steps{$s}++;}
$od||="$in/BMK_Result";

`mkdir $od`	unless(-d $od);
`cp $Bin/readme/1.Whole.pdf $od/Readme.pdf`;
$od=abs_path($od);
$in=abs_path($in);

my ($odir,$dir,%config,@diff);
&read_Config("$detail_cfg");
if(exists $steps{1}){
print "Begin BMK_1_rawData\n";
$odir="$od/BMK_1_rawData";
`rm -r $odir`	if(-d "$odir");
`mkdir $odir`	unless(-d "$odir");
`mkdir $odir/BMK_1_Data_Assess`  unless(-d "$odir/BMK_1_Data_Assess");
`cp $Bin/readme/BMK_1.1_Data_Assess.pdf $odir/BMK_1_Data_Assess/Readme.pdf`;

`cp -r $in/Data_Assess/lncRNA/PNG $odir/BMK_1_Data_Assess/nonRib_PNG`;
`cp -r $in/Data_Assess/miRNA/PNG $odir/BMK_1_Data_Assess/miRNA_PNG`;
`rm $odir/BMK_1_Data_Assess/nonRib_PNG/*.cycleQ.png`;
`rm $odir/BMK_1_Data_Assess/miRNA_PNG/*.cycleQ.png`;
`rm $odir/BMK_1_Data_Assess/miRNA_PNG/*.acgtn.png`;
`rm $odir/BMK_1_Data_Assess/miRNA_PNG/*.quality.png`;
`cp $in/Data_Assess/miRNA/All_sample_filter.stat $odir/BMK_1_Data_Assess/miRNA.All_sample_filter.stat.xls`;
`cp $in/Data_Assess/miRNA/QC/AllSample_GC_Q.stat $odir/BMK_1_Data_Assess/miRNA.AllSample_GC_Q.stat.xls`;
`cp $in/Data_Assess/lncRNA/AllSample_GC_Q.stat $odir/BMK_1_Data_Assess/nonRib.AllSample_GC_Q.stat.xls`;

$dir="$in/Basic_Analysis";

`mkdir $odir/BMK_2_Mapped_Statistics`		unless(-d "$odir/BMK_2_Mapped_Statistics");
`mkdir $odir/BMK_2_Mapped_Statistics/lncRNA`	unless(-d "$odir/BMK_2_Mapped_Statistics/lncRNA");
`cp $Bin/readme/BMK_1.2_lncRNA_map.pdf $odir/BMK_2_Mapped_Statistics/lncRNA/Readme.pdf`;
`cp $dir/Hisat_Stringtie/Mapped/*map.png $odir/BMK_2_Mapped_Statistics/lncRNA`;
`cp $dir/Hisat_Stringtie/Mapped/*type.png $odir/BMK_2_Mapped_Statistics/lncRNA`;
`cp $dir/Hisat_Stringtie/Mapped/*.mappedStat.xls $odir/BMK_2_Mapped_Statistics/lncRNA`;

`mkdir $odir/BMK_2_Mapped_Statistics/circRNA`	unless(-d "$odir/BMK_2_Mapped_Statistics/circRNA");
`cp $Bin/readme/BMK_1.3_circRNA_map.pdf $odir/BMK_2_Mapped_Statistics/circRNA/Readme.pdf`;
`cp $dir/circRNA_Bwa/Map_Stat/*readDensity.png $odir/BMK_2_Mapped_Statistics/circRNA`;
`cp $dir/circRNA_Bwa/Map_Stat/*.type.png $odir/BMK_2_Mapped_Statistics/circRNA`;
`cp $dir/circRNA_Bwa/Map_Stat/Total_Type.png $odir/BMK_2_Mapped_Statistics/circRNA`;
`cp $dir/circRNA_Bwa/Map_Stat/*.mappedStat.xls $odir/BMK_2_Mapped_Statistics/circRNA`;

`mkdir $odir/BMK_2_Mapped_Statistics/miRNA`	unless(-d "$odir/BMK_2_Mapped_Statistics/miRNA");
`cp $Bin/readme/BMK_1.4_sRNA_map.pdf $odir/BMK_2_Mapped_Statistics/miRNA/Readme.pdf`;
`cp $dir/sRNA_Alignment/All_sample_map.stat $odir/BMK_2_Mapped_Statistics/miRNA/All_sample_map.stat.xls`;
my @sRNAs=glob("$dir/sRNA_Alignment/*/*.clean.fa");
my @samples=map{(split(/\./,basename $_))[0]} @sRNAs;
my %stat=();
my @types=("Total","rRNA","scRNA","snRNA","snoRNA","tRNA","Repbase","Unannotated");
open(OUT,">$odir/BMK_2_Mapped_Statistics/miRNA/All_ncRNA_mapped.stat.xls")||die $!;
print OUT "#ID\t",join("\t",@types),"\n";
foreach my $s(@samples){
	`mkdir $odir/BMK_2_Mapped_Statistics/miRNA/$s` unless(-d "$odir/BMK_2_Mapped_Statistics/miRNA/$s");
	`cp $dir/sRNA_Alignment/$s/$s.clean.fa $odir/BMK_2_Mapped_Statistics/miRNA/$s`;
	`cp $dir/sRNA_Alignment/$s/Len_stat/$s.clean_len.length.png $odir/BMK_2_Mapped_Statistics/miRNA/$s/$s.Clean_reads.length.png`;
	`cp $dir/sRNA_Alignment/$s/Len_stat/$s.clean_len.length.pdf $odir/BMK_2_Mapped_Statistics/miRNA/$s/$s.Clean_reads.length.pdf`;
	`cp $dir/sRNA_Alignment/$s/Len_stat/$s.map_genome.length.png $odir/BMK_2_Mapped_Statistics/miRNA/$s/$s.Genome_mapped_reads.length.png`;
	`cp $dir/sRNA_Alignment/$s/Len_stat/$s.map_genome.length.pdf $odir/BMK_2_Mapped_Statistics/miRNA/$s/$s.Genome_mapped_reads.length.pdf`;
	`cp $dir/sRNA_Alignment/$s/genome_map/$s.png $odir/BMK_2_Mapped_Statistics/miRNA/$s/$s.chro_distribution.png`;
	`cp $dir/sRNA_Alignment/$s/Len_stat/$s.clean_len.stat $odir/BMK_2_Mapped_Statistics/miRNA/$s/$s.Clean_reads.length.stat.xls`;
	`cp $dir/sRNA_Alignment/$s/Len_stat/$s.map_genome.stat $odir/BMK_2_Mapped_Statistics/miRNA/$s/$s.Genome_mapped_reads.length.stat.xls`;
	`cp $dir/sRNA_Alignment/$s/Data.stat	$odir/BMK_2_Mapped_Statistics/miRNA/$s/$s.ncRNA_mapped.stat.xls`;
	open(STAT,"$dir/sRNA_Alignment/$s/Data.stat")||die $!;
	while(<STAT>){
		chomp;my @tmp=split;
		$stat{$s}{$tmp[0]}="$tmp[1] ($tmp[2])";	
	}
	close(STAT);
	print OUT "$s";
	foreach my $t(@types){	print OUT "\t$stat{$s}{$t}";}
	print OUT "\n";
}
close(OUT);

`mkdir $odir/BMK_3_Library_Assessment`		unless(-d "$odir/BMK_3_Library_Assessment");
`cp $Bin/readme/BMK_1.5_Library_Assessment.pdf	$odir/BMK_3_Library_Assessment/Readme.pdf`;
`mkdir $odir/BMK_3_Library_Assessment/lncRNA`	unless(-d "$odir/BMK_3_Library_Assessment/lncRNA");
`cp $dir/Hisat_Stringtie/Lib_type_stat/fq_map.xls $odir/BMK_3_Library_Assessment/lncRNA/Lib_assess.xls`;
`cp $dir/Hisat_Stringtie/Mapped/*.randcheck.png $odir/BMK_3_Library_Assessment/lncRNA`;
`cp $dir/Hisat_Stringtie/Mapped/*.insertSize.png $odir/BMK_3_Library_Assessment/lncRNA`;
`cp $dir/Hisat_Stringtie/Mapped/*.Saturation.png $odir/BMK_3_Library_Assessment/lncRNA`;

`mkdir $odir/BMK_3_Library_Assessment/circRNA`	unless(-d "$odir/BMK_3_Library_Assessment/circRNA");
`cp $dir/circRNA_Bwa/Map_Stat/*.insertSize.png $odir/BMK_3_Library_Assessment/circRNA`;
}
if(exists $steps{2}){	###lncRNA
print "Begin BMK_2_LncRNA\n";

$odir="$od/BMK_2_LncRNA";
$dir="$in/Basic_Analysis";
`mkdir $odir`	unless(-d $odir);
my @gtfs=glob("$dir/Hisat_Stringtie/StringTie/*/StringTie_asm.gtf");
my @samples=map{(split(/\//,$_))[-2]} @gtfs;
`mkdir $odir/BMK_1_Assembly_Result` unless(-d "$odir/BMK_1_Assembly_Result");
`cp $Bin/readme/BMK_2.1_Assembly_Result.pdf $odir/BMK_1_Assembly_Result/Readme.pdf`;
foreach my $s(@samples){
	`cp $dir/Hisat_Stringtie/StringTie/$s/StringTie_asm.gtf $odir/BMK_1_Assembly_Result/$s.StringTie.transcripts.gtf`;
	`sed -i '1d' $odir/BMK_1_Assembly_Result/$s.StringTie.transcripts.gtf && sed -i '1i\#Chr\tSource\tType\tStart Site\tEnd Site\tScore\tStrand\tframe\tDescription' $odir/BMK_1_Assembly_Result/$s.StringTie.transcripts.gtf`;
}
#############
`mkdir $odir/BMK_2_LncRNA_Prediction`	unless(-d "$odir/BMK_2_LncRNA_Prediction");
`cp $Bin/readme/BMK_2.2_LncRNA_Prediction.pdf $odir/BMK_2_LncRNA_Prediction/Readme.pdf`;
`mkdir $odir/BMK_2_LncRNA_Prediction/BMK_1_Software_Result`	unless(-d "$odir/BMK_2_LncRNA_Prediction/BMK_1_Software_Result");
`cp $Bin/readme/BMK_2.3_four_softwares.pdf $odir/BMK_2_LncRNA_Prediction/BMK_1_Software_Result/Readme.pdf`;
`cp $dir/LncRNA_Analysis/lnc_predict/code_filter/CNCI/CNCI.index $odir/BMK_2_LncRNA_Prediction/BMK_1_Software_Result/CNCI.xls`;
`cp $dir/LncRNA_Analysis/lnc_predict/code_filter/CPAT/cpat.txt $odir/BMK_2_LncRNA_Prediction/BMK_1_Software_Result/CPAT.xls`;
`cp $dir/LncRNA_Analysis/lnc_predict/code_filter/CPC/lnc_code_filter.result.txt $odir/BMK_2_LncRNA_Prediction/BMK_1_Software_Result/CPC.xls`;
`cp $dir/LncRNA_Analysis/lnc_predict/code_filter/Pfam/Pfam_result.txt $odir/BMK_2_LncRNA_Prediction/BMK_1_Software_Result/Pfam.xls`;
`cp $dir/LncRNA_Analysis/lnc_predict/code_filter/list.txt $odir/BMK_2_LncRNA_Prediction/BMK_1_Software_Result/Software_veen.xls`;
`cp $dir/LncRNA_Analysis/lnc_predict/code_filter/venn.png $odir/BMK_2_LncRNA_Prediction/BMK_1_Software_Result/venn.png`;

`cp $dir/LncRNA_Analysis/lnc_predict/All_filter_final.fa $odir/BMK_2_LncRNA_Prediction/LncRNA.fa`;
`cp $dir/LncRNA_Analysis/lnc_predict/All_filter_final.gff $odir/BMK_2_LncRNA_Prediction/LncRNA.gff`;
`cp $dir/LncRNA_Analysis/lnc_predict/All_filter_final.gtf $odir/BMK_2_LncRNA_Prediction/LncRNA.gtf`;
`cp $dir/LncRNA_Analysis/lnc_predict/circos/circos.png $odir/BMK_2_LncRNA_Prediction/LncRNA_circos.png`;
`cp $dir/LncRNA_Analysis/lnc_predict/circos/circos.svg $odir/BMK_2_LncRNA_Prediction/LncRNA_circos.svg`;
`cp $dir/LncRNA_Analysis/lnc_predict/Lnc_filter/lncRNA_filter_pie.stat.png $odir/BMK_2_LncRNA_Prediction/LncRNA_classification.png`;
`cp $dir/LncRNA_Analysis/lnc_predict/Lnc_filter/lnc_class_code.xls $odir/BMK_2_LncRNA_Prediction/LncRNA_classification.xls`;
`cp $dir/LncRNA_Analysis/lnc_predict/Lnc_filter/lncRNA_class_id.list $odir/BMK_2_LncRNA_Prediction/lncRNA_classification_name.xls`;
############
`mkdir $odir/BMK_3_LncRNA_Expression`		unless(-d "$odir/BMK_3_LncRNA_Expression");
`cp $Bin/readme/BMK_2.4_LncRNA_Expression.pdf $odir/BMK_3_LncRNA_Expression/Readme.pdf`;
`cp $dir/Hisat_Stringtie/prepDE/All_lncRNA_counts.list $odir/BMK_3_LncRNA_Expression/lncRNA_counts.xls`;
`cp $dir/Hisat_Stringtie/prepDE/All_lncRNA_fpkm.list $odir/BMK_3_LncRNA_Expression/lncRNA_expression.xls`;

`mkdir $odir/BMK_3_LncRNA_Expression/PNG`	unless(-d "$odir/BMK_3_LncRNA_Expression/PNG");
`cp $Bin/readme/BMK_2.5_corr.pdf $odir/BMK_3_LncRNA_Expression/PNG/Readme.pdf`;
&getDensity("$in/DEG_Analysis/lncRNA/density","$odir/BMK_3_LncRNA_Expression/PNG");

############
`mkdir $odir/BMK_4_LncRNA_Target`	unless(-d "$odir/BMK_4_LncRNA_Target");
`cp $Bin/readme/BMK_2.6_LncRNA_Target.pdf $odir/BMK_4_LncRNA_Target/Readme.pdf`;
`cp $dir/LncRNA_Analysis/Lnc_target_predict/Cis_target_gene.xls $odir/BMK_4_LncRNA_Target`;
`cp $dir/LncRNA_Analysis/Lnc_target_predict/Trans/Trans_target_gene.xls $odir/BMK_4_LncRNA_Target`	if(-e "$dir/LncRNA_Analysis/Lnc_target_predict/Trans/Trans_target_gene.xls");
`rm -r $odir/BMK_4_LncRNA_Target/WGCNA_Result`	if(-d "$odir/BMK_4_LncRNA_Target/WGCNA_Result");
`cp -r $dir/Tophat_Cufflinks/Lnc_target_predict/Trans_target $odir/BMK_4_LncRNA_Target/WGCNA_Result`	if(-d "$dir/Tophat_Cufflinks/Lnc_target_predict/Trans_target");

############
&getDEG("$in/DEG_Analysis/lncRNA","$odir/BMK_5_DEG_Analysis");
`cp $Bin/readme/BMK_2.7_DEG_Analysis.pdf $odir/BMK_5_DEG_Analysis/Readme.pdf`;
&DEG_anno_stat("$in/DEG_Analysis/lncRNA","Cis_Anno_enrichment","$odir/BMK_5_DEG_Analysis/DEG_Cis_anno.stat.xls");
&DEG_anno_stat("$in/DEG_Analysis/lncRNA","Trans_Anno_enrichment","$odir/BMK_5_DEG_Analysis/DEG_Trans_anno.stat.xls");


if(-d "$dir/LncRNA_Analysis/Lnc_disease"){
	`mkdir $odir/BMK_6_Lnc_disease`	unless(-d "$odir/BMK_6_Lnc_disease");
	`cp $Bin/readme/BMK_2.8_Lnc_disease.pdf $odir/BMK_6_Lnc_disease/Readme.pdf`;
	`cp $dir/LncRNA_Analysis/Lnc_disease/LncRNA_disease.txt $odir/BMK_6_Lnc_disease/LncRNA_disease.xls`;

}
`mkdir $odir/BMK_7_Lnc_precursor` unless(-d "$odir/BMK_7_Lnc_precursor");
`cp $Bin/readme/BMK_2.9_Lnc_precursor.pdf $odir/BMK_7_Lnc_precursor/Readme.pdf`;
`cp $in/Combine/lncRNA_precursor/lncRNA_precursor.txt $odir/BMK_7_Lnc_precursor/lncRNA_precursor.xls`;
}
if(exists $steps{3}){###mRNA
print "Begin BMK_3_mRNA\n";
$odir="$od/BMK_3_mRNA";
`mkdir $odir`   unless(-d $odir);
$dir="$in/Anno_Integrate/New_Anno/Result";
`mkdir $odir/BMK_1_NewGene`		unless(-d "$odir/BMK_1_NewGene");
`cp $Bin/readme/BMK_3.1_NewGene.pdf $odir/BMK_1_NewGene/Readme.pdf`;
`cp $in/Basic_Analysis/Hisat_Stringtie/genePredict/final_track/*.newGene_final.filtered.gff $odir/BMK_1_NewGene`;
`cp $in/Basic_Analysis/Hisat_Stringtie/genePredict/final_track/*.newGene.longest_transcript.fa $odir/BMK_1_NewGene`;

`mkdir -p $odir/BMK_1_NewGene/BMK_1_NewGene_Anno/BMK_1_Annotation`      unless(-d "$odir/BMK_1_NewGene/BMK_1_NewGene_Anno/BMK_1_Annotation");
`cp $dir/All_Database_annotation.xls $odir/BMK_1_NewGene/BMK_1_NewGene_Anno`;
`cp $dir/Integrated_Function.annotation.xls $odir/BMK_1_NewGene/BMK_1_NewGene_Anno`;
`cp $dir/Function_Annotation.stat.xls $odir/BMK_1_NewGene/BMK_1_NewGene_Anno`;

`cp $Bin/readme/BMK_3.2_NewGene_anno.pdf $odir/BMK_1_NewGene/BMK_1_NewGene_Anno/Readme.pdf`;
my @anno_txt = glob "$dir/*txt";
foreach my $file(@anno_txt){
	my $File = basename $file;
	$File =~s/\.txt/\.xls/;
	`cp $file $odir/BMK_1_NewGene/BMK_1_NewGene_Anno/BMK_1_Annotation/$File`;
}
my @kegg = glob "$dir/*ko";
push @kegg,glob "$dir/*pathway";
foreach my $kfile (@kegg){
	my $KFile = basename $kfile;
	`cp $kfile $odir/BMK_1_NewGene/BMK_1_NewGene_Anno/BMK_1_Annotation/$KFile.xls`;
}
`cp $dir/*.GO_tree.stat.xls $odir/BMK_1_NewGene/BMK_1_NewGene_Anno/BMK_1_Annotation`;

`mkdir $odir/BMK_1_NewGene/BMK_1_NewGene_Anno/BMK_2_Statistic`      unless(-d "$odir/BMK_1_NewGene/BMK_1_NewGene_Anno/BMK_2_Statistic");
`cp $dir/*png $odir/BMK_1_NewGene/BMK_1_NewGene_Anno/BMK_2_Statistic`;
my @stat = glob "$dir/*stat";
foreach my $stat(@stat){
	my $Stat = basename $stat;
	`cp $stat $odir/BMK_1_NewGene/BMK_1_NewGene_Anno/BMK_2_Statistic/$Stat.xls`;
}

`rm -r $odir/BMK_1_NewGene/BMK_1_NewGene_Anno/BMK_3_KEGG_map`	if(-d "$odir/BMK_1_NewGene/BMK_1_NewGene_Anno/BMK_3_KEGG_map");
`cp -r $dir/Kegg_map $odir/BMK_1_NewGene/BMK_1_NewGene_Anno/BMK_3_KEGG_map`;


$dir="$in/Basic_Analysis";
`mkdir $odir/BMK_2_geneExpression`	unless(-d "$odir/BMK_2_geneExpression");
`cp $Bin/readme/BMK_3.3_geneExpression.pdf $odir/BMK_2_geneExpression/Readme.pdf`;
`cp $dir/Hisat_Stringtie/prepDE/All_gene_counts.list $odir/BMK_2_geneExpression/gene_counts.xls`;
`cp $dir/Hisat_Stringtie/prepDE/All_gene_fpkm.list $odir/BMK_2_geneExpression/gene_expression.xls`;

`mkdir $odir/BMK_2_geneExpression/PNG`      unless(-d "$odir/BMK_2_geneExpression/PNG");
`cp $Bin/readme/BMK_3.4_gene_corr.pdf $odir/BMK_2_geneExpression/PNG/Readme.pdf`;
&getDensity("$in/DEG_Analysis/gene/density","$odir/BMK_2_geneExpression/PNG");

if(-d "$in/DEG_Analysis/gene/kmeans"){
	`mkdir $odir/BMK_2_geneExpression/Kmeans_Cluster` unless (-d "$odir/BMK_2_geneExpression/Kmeans_Cluster");
	`cp $in/DEG_Analysis/gene/kmeans/k-means.png $odir/BMK_2_geneExpression/Kmeans_Cluster`;
	`cp $in/DEG_Analysis/gene/kmeans/kmeans_cluster.txt $odir/BMK_2_geneExpression/Kmeans_Cluster/kmeans_cluster.xls`;
	`cp $in/DEG_Analysis/gene/kmeans/cluster_gene.xls $odir/BMK_2_geneExpression/Kmeans_Cluster/`;
}

&getDEG("$in/DEG_Analysis/gene","$odir/BMK_3_DEG_Analysis","gene");
`cp $Bin/readme/BMK_3.5_DEG_Analysis.pdf $odir/BMK_3_DEG_Analysis/Readme.pdf`;
`cp $Bin/readme/BMK_3.5_Kmeans.pdf  $odir/BMK_2_geneExpression/Kmeans_Cluster/Readme.pdf`;
&DEG_anno_stat("$in/DEG_Analysis/gene","Anno_enrichment","$odir/BMK_3_DEG_Analysis/DEG_anno.stat.xls");

}

if(exists $steps{4}){###circRNA
print "Begin BMK_4_circRNA\n";
$odir="$od/BMK_4_circRNA";
`mkdir $odir`   unless(-d $odir);
$dir="$in/Basic_Analysis/circRNA_analysis";
`mkdir $odir/BMK_1_circRNA_Prediction`   unless(-d "$odir/BMK_1_circRNA_Prediction");
`cp $Bin/readme/BMK_4.1_circRNA_Prediction.pdf $odir/BMK_1_circRNA_Prediction/Readme.pdf`;
`cp $dir/circRNA.fa $odir/BMK_1_circRNA_Prediction`;
`cp $dir/new_name/circRNA_newname.xls $odir/BMK_1_circRNA_Prediction`;
`cp $dir/overlap_alitisplice/overlap_alitisplice.xls $odir/BMK_1_circRNA_Prediction`;

`mkdir $odir/BMK_1_circRNA_Prediction/BMK_1_CIRC_info`   unless(-d "$odir/BMK_1_circRNA_Prediction/BMK_1_CIRC_info");
`cp $Bin/readme/BMK_4.2_circRNA_single_sample.pdf $odir/BMK_1_circRNA_Prediction/BMK_1_CIRC_info/Readme.pdf`;
my $circ_software=0;
my @CIRI=glob("$dir/circRNA_identify/*.circ");
if(@CIRI>=1){
	$circ_software++;
	`mkdir $odir/BMK_1_circRNA_Prediction/BMK_1_CIRC_info/BMK_$circ_software\_CIRI` unless (-d "$odir/BMK_1_circRNA_Prediction/BMK_1_CIRC_info/BMK_$circ_software\_CIRI");
	foreach my $c(@CIRI){
		my $ss=(split /\./,(split(/\//,$c))[-1])[0];
		`cp $dir/circRNA_identify/$ss.circ $odir/BMK_1_circRNA_Prediction/BMK_1_CIRC_info/BMK_$circ_software\_CIRI/$ss.CIRI_Predict.xls`;
	}
}
my @intersect=glob("$dir/circRNA_identify/*.circ.intersect");
foreach my $cc(@intersect){
	my $ff=(split/\./,(split(/\//,$cc))[-1])[0];
	`cp $cc $odir/BMK_1_circRNA_Prediction/BMK_1_CIRC_info/$ff.circ.intersect.xls`;
}

if(-d "$dir/circRNA_identify/find_circ/"){
	$circ_software++;
	`mkdir $odir/BMK_1_circRNA_Prediction/BMK_1_CIRC_info/BMK_$circ_software\_find_circ` unless (-d "$odir/BMK_1_circRNA_Prediction/BMK_1_CIRC_info/BMK_$circ_software\_find_circ");
	my @find_circ=glob("$dir/circRNA_identify/find_circ/*/circ_candidates.bed");
	foreach my $f(@find_circ){
		my $s=(split(/\//,$f))[-2];
		`cp $f $odir/BMK_1_circRNA_Prediction/BMK_1_CIRC_info/BMK_$circ_software\_find_circ/$s.find_circ_Predict.xls`;
	}
}
`mkdir $odir/BMK_1_circRNA_Prediction/BMK_2_statistics`   unless(-d "$odir/BMK_1_circRNA_Prediction/BMK_2_statistics");
`cp $Bin/readme/BMK_4.3_circRNA_statistic.pdf $odir/BMK_1_circRNA_Prediction/BMK_2_statistics/Readme.pdf`;
`cp $dir/statistics/*png $odir/BMK_1_circRNA_Prediction/BMK_2_statistics`;
`cp $dir/statistics/*pdf $odir/BMK_1_circRNA_Prediction/BMK_2_statistics`;
`cp $dir/statistics/statistic.xls $odir/BMK_1_circRNA_Prediction/BMK_2_statistics`;
if(-e "$dir/circos/circosResult/circos.png"){
	`cp $dir/circos/circosResult/circos.png $odir/BMK_1_circRNA_Prediction/BMK_2_statistics/circRNA_circos.png`;
	`cp $dir/circos/circosResult/circos.svg $odir/BMK_1_circRNA_Prediction/BMK_2_statistics/circRNA_circos.svg`;
}
if(-e "$dir/statistics/Known/ratio.out"){
	`mkdir $odir/BMK_1_circRNA_Prediction/BMK_3_Known`	unless(-d "$odir/BMK_1_circRNA_Prediction/BMK_3_Known");
	`cp $Bin/readme/BMK_4.4_known_circRNA.pdf $odir/BMK_1_circRNA_Prediction/BMK_3_Known/Readme.pdf`;
	`cp $dir/statistics/Known/ratio.out $odir/BMK_1_circRNA_Prediction/BMK_3_Known/circRNA_known.xls`;
	`cp $dir/statistics/Known/circ_pred.pie.png $odir/BMK_1_circRNA_Prediction/BMK_3_Known/Known_Unknown_pie.png`;
	`cp $dir/statistics/Known/forpie.list $odir/BMK_1_circRNA_Prediction/BMK_3_Known/Known_Unknown_pie.xls`;
	`cp $dir/statistics/Known/new.list $odir/BMK_1_circRNA_Prediction/BMK_3_Known/new.list.xls`;
	if(exists $config{medical} && $config{medical}=~/GRCh37|GRCh38/){
		`cp $dir/statistics/Known/circRNA_Circ2disease_annotation.xls $odir/BMK_1_circRNA_Prediction/BMK_3_Known/`;
		`cp $dir/statistics/Known/circRNA_circBank_annotation.xls $odir/BMK_1_circRNA_Prediction/BMK_3_Known/`;
	}
}


`mkdir $odir/BMK_2_circRNA_Expression`		unless(-d "$odir/BMK_2_circRNA_Expression");
`cp $Bin/readme/BMK_4.5_circRNA_Expression.pdf $odir/BMK_2_circRNA_Expression/Readme.pdf`;
`mkdir $odir/BMK_2_circRNA_Expression/PNG`	unless(-d "$odir/BMK_2_circRNA_Expression/PNG");
`cp $Bin/readme/BMK_4.6_circRNA_corr.pdf $odir/BMK_2_circRNA_Expression/PNG/Readme.pdf`;
`cp $dir/expression/All_gene_counts.list $odir/BMK_2_circRNA_Expression/circRNA_counts.xls`;
`cp $dir/expression/All_gene_expression.list $odir/BMK_2_circRNA_Expression/circRNA_expression.xls`;
&getDensity("$in/DEG_Analysis/circRNA/density","$odir/BMK_2_circRNA_Expression/PNG");

&getDEG("$in/DEG_Analysis/circRNA","$odir/BMK_3_DEG_Analysis");
`cp $Bin/readme/BMK_4.7_DEG_Analysis.pdf $odir/BMK_3_DEG_Analysis/Readme.pdf`;
&DEG_anno_stat("$in/DEG_Analysis/circRNA","Anno_enrichment","$odir/BMK_3_DEG_Analysis/DEG_anno.stat.xls");

}
if(exists $steps{5}){###miRNA
print "Begin BMK_5_miRNA\n";
$odir="$od/BMK_5_miRNA";
`mkdir $odir`   unless(-d $odir);
$dir="$in/Basic_Analysis/sRNA_Analysis/miRDeep2/miRNA_Quantify";
`mkdir $odir/BMK_1_miRNA_Prediction`   unless(-d "$odir/BMK_1_miRNA_Prediction");
`cp $Bin/readme/BMK_5.1_miRNA_Prediction.pdf $odir/BMK_1_miRNA_Prediction/Readme.pdf`;
`mkdir $odir/BMK_1_miRNA_Prediction/BMK_1_miR_Seq`	unless(-d "$odir/BMK_1_miRNA_Prediction/BMK_1_miR_Seq");
`cp $dir/*fa $odir/BMK_1_miRNA_Prediction/BMK_1_miR_Seq`;

`mkdir $odir/BMK_1_miRNA_Prediction/BMK_2_Len_Distribution`	unless(-d "$odir/BMK_1_miRNA_Prediction/BMK_2_Len_Distribution");
`cp $dir/*length.png $odir/BMK_1_miRNA_Prediction/BMK_2_Len_Distribution`;
`cp $dir/*length.pdf $odir/BMK_1_miRNA_Prediction/BMK_2_Len_Distribution`;

my @S = glob "$dir/*.len.stat";
push @S,glob "$dir/*miRNA.stat";
foreach my $s (@S){
	my $ss = basename $s;
	`cp $s $odir/BMK_1_miRNA_Prediction/BMK_2_Len_Distribution/$ss.xls`;
}
`rm $odir/BMK_1_miRNA_Prediction/BMK_2_Len_Distribution/Total_miRNA.stat.xls`;

`mkdir $odir/BMK_1_miRNA_Prediction/BMK_3_Base_Distribution`	unless(-d "$odir/BMK_1_miRNA_Prediction/BMK_3_Base_Distribution");
my @base = glob "$dir/*bias.stat";
push @base,glob "$dir/*base.stat";
foreach my $c (@base){
	my $cc = basename $c;
	`cp $c $odir/BMK_1_miRNA_Prediction/BMK_3_Base_Distribution/$cc.xls`;
}
`cp $dir/*bias.png $odir/BMK_1_miRNA_Prediction/BMK_3_Base_Distribution`;
`cp $dir/*bias.pdf $odir/BMK_1_miRNA_Prediction/BMK_3_Base_Distribution`;
`cp $dir/*base.pdf* $odir/BMK_1_miRNA_Prediction/BMK_3_Base_Distribution`;
`cp $dir/*base.png* $odir/BMK_1_miRNA_Prediction/BMK_3_Base_Distribution`;

`mkdir $odir/BMK_1_miRNA_Prediction/BMK_4_Summary_Stat`	unless(-d "$odir/BMK_1_miRNA_Prediction/BMK_4_Summary_Stat");
`cp $dir/sum_known_miRNA.txt $odir/BMK_1_miRNA_Prediction/BMK_4_Summary_Stat/Summary_Known_miRNA.xls`;
`cp $dir/sum_predicted_miRNA.txt $odir/BMK_1_miRNA_Prediction/BMK_4_Summary_Stat/Summary_Novel_miRNA.xls`;
`cp $dir/Total_miRNA.stat $odir/BMK_1_miRNA_Prediction/BMK_4_Summary_Stat/Total_miRNA.stat.xls`;

`mkdir $odir/BMK_1_miRNA_Prediction/BMK_5_PDF`		unless(-d "$odir/BMK_1_miRNA_Prediction/BMK_5_PDF");
`cp $dir/*/*pdf $odir/BMK_1_miRNA_Prediction/BMK_5_PDF`;

###################
`mkdir $odir/BMK_2_miRNA_Expression`      unless(-d "$odir/BMK_2_miRNA_Expression");
`cp $Bin/readme/BMK_5.2_miRNA_Expression.pdf $odir/BMK_2_miRNA_Expression/Readme.pdf`;
`cp $dir/All_miRNA_expression.list $odir/BMK_2_miRNA_Expression/miRNA_expression.xls`;
`cp $dir/All_miRNA.count.list $odir/BMK_2_miRNA_Expression/miRNA_counts.xls`;

`mkdir $odir/BMK_2_miRNA_Expression/PNG`      unless(-d "$odir/BMK_2_miRNA_Expression/PNG");
`cp $Bin/readme/BMK_5.3_miRNA_corr.pdf $odir/BMK_2_miRNA_Expression/PNG/Readme.pdf`;
&getDensity("$in/DEG_Analysis/sRNA/density","$odir/BMK_2_miRNA_Expression/PNG");

###################
&getDEG("$in/DEG_Analysis/sRNA","$odir/BMK_3_DEG_Analysis");
`cp $Bin/readme/BMK_5.4_DEG_Analysis.pdf $odir/BMK_3_DEG_Analysis/Readme.pdf`;
&DEG_anno_stat("$in/DEG_Analysis/sRNA","Anno_enrichment","$odir/BMK_3_DEG_Analysis/DEG_anno.stat.xls");

###################
`mkdir $odir/BMK_4_miRNA_Family`	unless(-d "$odir/BMK_4_miRNA_Family");
`cp $Bin/readme/BMK_5.5_miRNA_Family.pdf $odir/BMK_4_miRNA_Family/Readme.pdf`;
`cp $in/miRNA_Family/family.xls $odir/BMK_4_miRNA_Family/Family_In_All.xls`;
`cp $in/miRNA_Family/result.txt $odir/BMK_4_miRNA_Family/Family_In_miR.xls`;

###################
`mkdir $odir/BMK_5_miRNA_Edit`		unless(-d "$odir/BMK_5_miRNA_Edit");
`cp $Bin/readme/BMK_5.6_miRNA_Edit.pdf $odir/BMK_5_miRNA_Edit/Readme.pdf`;
`cp $in/miRNA_edit/Results/MapResults/cutoff_20_MapResults.txt $odir/BMK_5_miRNA_Edit/miRNAEdit_cutoff_20.xls`;

###################
if(-d "$in/Basic_Analysis/sRNA_Analysis/HMDD"){
	`mkdir $odir/BMK_6_HMDD`	unless(-d "$odir/BMK_6_HMDD");
	`cp $Bin/readme/BMK_5.7_HMDD.pdf $odir/BMK_6_HMDD/Readme.pdf`;
	`cp $in/Basic_Analysis/sRNA_Analysis/HMDD/HMDD_related_miRNA.txt $odir/BMK_6_HMDD/HMDD_related_miRNA.xls`;
	`cp $in/Basic_Analysis/sRNA_Analysis/HMDD/HMDD_related_miRNA_with_target_gene.txt $odir/BMK_6_HMDD/HMDD_related_miRNA_with_target_gene.xls`;
}
`cp $in/Basic_Analysis/sRNA_Analysis/miRNA_loc/miRNA_pos.list $od/`;
`cp $in/Basic_Analysis/sRNA_Analysis/miRNA_loc/miRNA_mature.fa $od/`;
`cp $in/Basic_Analysis/sRNA_Analysis/miRNA_loc/miRNA_pre_seq.fa $od/`;
}
if(exists $steps{6}){
print "Begin BMK_6_Structure\n";
$odir="$od/BMK_6_Structure";
`mkdir $odir`   unless(-d $odir);
###############
$dir="$in/SNP_Analysis";
`mkdir $odir/BMK_1_SNP_Analysis`			unless(-d "$odir/BMK_1_SNP_Analysis");
`cp $Bin/readme/BMK_6.1_SNP_Analysis.pdf $odir/BMK_1_SNP_Analysis/Readme.pdf`;
`cp $dir/stat/AllSample.snp.stat $odir/BMK_1_SNP_Analysis/AllSample.snp.stat.xls`;
`cp $dir/stat/AllSample.SNP_density.stat $odir/BMK_1_SNP_Analysis/AllSample.SNP_density.stat.xls`;
`cp $dir/stat/AllSample*png $odir/BMK_1_SNP_Analysis/`;
`cp $dir/stat/final.snp.anno.gatk.all.list $odir/BMK_1_SNP_Analysis/final.snp.anno.gatk.all.list.xls`;
`cp $dir/stat/final.indel.anno.gatk.all.list $odir/BMK_1_SNP_Analysis/final.indel.anno.gatk.all.list.xls`;
##
`mkdir $odir/BMK_1_SNP_Analysis/BMK_1_InDel_anno`	unless(-d "$odir/BMK_1_SNP_Analysis/BMK_1_InDel_anno");
`cp $dir/stat/indel_anno/*png $odir/BMK_1_SNP_Analysis/BMK_1_InDel_anno`;
`cp $dir/stat/indel_anno/final_Indel.anno.stat $odir/BMK_1_SNP_Analysis/BMK_1_InDel_anno/final_Indel.anno.stat.xls`;
##
`mkdir $odir/BMK_1_SNP_Analysis/BMK_2_SNP_anno`		unless(-d "$odir/BMK_1_SNP_Analysis/BMK_2_SNP_anno");
`cp $dir/stat/snp_anno/*png $odir/BMK_1_SNP_Analysis/BMK_2_SNP_anno`;
`cp $dir/stat/snp_anno/final_SNP.anno.stat $odir/BMK_1_SNP_Analysis/BMK_2_SNP_anno/final_SNP.anno.stat.xls`;
##
`mkdir $odir/BMK_1_SNP_Analysis/BMK_3_SNP_type`		unless(-d "$odir/BMK_1_SNP_Analysis/BMK_3_SNP_type");
`cp $dir/stat/*.snp.type.png $odir/BMK_1_SNP_Analysis/BMK_3_SNP_type`;
`cp $dir/stat/All.snp_type.stat $odir/BMK_1_SNP_Analysis/BMK_3_SNP_type/All.snp_type.stat.xls`;

################
$dir="$in/miRNA_Target";
`mkdir $odir/BMK_2_miRNA_Target`		unless(-d "$odir/BMK_2_miRNA_Target");
`cp $Bin/readme/BMK_6.2_miRNA_Target.pdf $odir/BMK_2_miRNA_Target/Readme.pdf`;
`cp $dir/circRNA/mir2target.list $odir/BMK_2_miRNA_Target/circRNA.mir2target.list.xls`;
`cp $dir/lncRNA/mir2target.list $odir/BMK_2_miRNA_Target/lncRNA.mir2target.list.xls`;
`cp $dir/gene/mir2target.list $odir/BMK_2_miRNA_Target/gene.mir2target.list.xls`;
`cp $dir/circRNA/mir2target.stat $odir/BMK_2_miRNA_Target/circRNA.mir2target.stat.xls`;
`cp $dir/lncRNA/mir2target.stat $odir/BMK_2_miRNA_Target/lncRNA.mir2target.stat.xls`;
`cp $dir/gene/mir2target.stat $odir/BMK_2_miRNA_Target/gene.mir2target.stat.xls`;
`cp $dir/circRNA/miRNA_target_info.xls $odir/BMK_2_miRNA_Target/circRNA.miRNA_target_info.xls`	if(-e "$dir/circRNA/miRNA_target_info.xls");
`cp $dir/lncRNA/miRNA_target_info.xls $odir/BMK_2_miRNA_Target/lncRNA.miRNA_target_info.xls`	if(-e "$dir/lncRNA/miRNA_target_info.xls");
`cp $dir/gene/miRNA_target_info.xls $odir/BMK_2_miRNA_Target/gene.miRNA_target_info.xls`	if(-e "$dir/gene/miRNA_target_info.xls");

################
`mkdir $odir/BMK_3_Gene_Structure_Optimize`	unless(-d "$odir/BMK_3_Gene_Structure_Optimize");
`cp $Bin/readme/BMK_6.3_Gene_Structure_Optimize.pdf $odir/BMK_3_Gene_Structure_Optimize/Readme.pdf`;
`cp $in/Gene_Structure_Optimize/*.geneStructure.optimize.xls $odir/BMK_3_Gene_Structure_Optimize`;

################
`mkdir $odir/BMK_4_Alt_splice`			unless(-d "$odir/BMK_4_Alt_splice");
`cp $Bin/readme/BMK_6.4_Alt_splice.pdf $odir/BMK_4_Alt_splice/Readme.pdf`;
$dir="$in/Alitsplice_Analysis";
`cp $dir/*png $odir/BMK_4_Alt_splice/`;
my @as=glob("$dir/*fpkm");
foreach my $a(@as){
        my $base=(split(/\./,basename $a))[0];
	print "$base\n";
        `awk 'OFS="\t"{\$NF="";print}' $dir/$base.fpkm| sed 's/\t\$//g' >$odir/BMK_4_Alt_splice/$base.AS.list.xls `;
}       

################
if(-d "$in/Gene_Fusion"){
	`mkdir $odir/BMK_5_Gene_Fusion`		unless(-d "$odir/BMK_5_Gene_Fusion");
	`cp $Bin/readme/BMK_6.5_Gene_Fusion.pdf $odir/BMK_5_Gene_Fusion/Readme.pdf`;
	`cp $in/Gene_Fusion/result/png/* $odir/BMK_5_Gene_Fusion`;
	`cp $in/Gene_Fusion/result/report/* $odir/BMK_5_Gene_Fusion`;
}

}
if(exists $steps{7}){
	print "Begin BMK_7_Combine\n";
	$odir="$od/BMK_7_Combine";	
	`mkdir $odir` unless(-d $odir);
	$dir="$in/Combine/Global";
	`mkdir $odir/BMK_1_Global` unless(-d "$odir/BMK_1_Global");
	`mkdir $odir/BMK_1_Global/BMK_1_All` unless(-d "$odir/BMK_1_Global/BMK_1_All");
	`cp $Bin/readme/BMK_7.1_Global.pdf $odir/BMK_1_Global/BMK_1_All/Readme.pdf`;
	opendir(DIR,"$dir/barplot")||die $!;
	my @samples=grep{/^[^\.]/ && -d "$dir/barplot/$_"} readdir(DIR);
	closedir(DIR);
	foreach my $s(@samples){
		`mkdir $odir/BMK_1_Global/BMK_1_All/$s`	unless(-d "$odir/BMK_1_Global/BMK_1_All/$s");
		`cp $dir/barplot/$s/$s* $odir/BMK_1_Global/BMK_1_All/$s && cp $dir/barplot/$s/barplot_t.txt $odir/BMK_1_Global/BMK_1_All/$s/$s.barplot.xls`;
		`cp -r $dir/circos/$s/circos.png $odir/BMK_1_Global/BMK_1_All/$s/$s.circos.png`;
		`cp -r $dir/circos/$s/circos.svg $odir/BMK_1_Global/BMK_1_All/$s/$s.circos.svg`;
		my @a = glob "$dir/circos/$s/*.txt";
		foreach my $a(@a){
			my $aa = basename $a;
			my $AA = (split /\./,$aa)[0];
			`cp -r $dir/circos/$s/$AA.txt $odir/BMK_1_Global/BMK_1_All/$s/$AA.xls`;
		}
	}
	`cp $dir/pie/Total.* $odir/BMK_1_Global/BMK_1_All/`;

	$dir="$in/Combine/Diff_stat";
	`mkdir $odir/BMK_1_Global/BMK_2_DEG` unless(-d "$odir/BMK_1_Global/BMK_2_DEG");
	`cp $Bin/readme/BMK_7.2_Global_diff.pdf $odir/BMK_1_Global/BMK_2_DEG/Readme.pdf`;
	opendir(DIR,$dir)||die $!;
	my @vs=grep{/_vs_/ && -d "$dir/$_"} readdir(DIR);
	closedir(DIR);
	foreach my $v(@vs){
		`mkdir $odir/BMK_1_Global/BMK_2_DEG/$v`	unless(-d "$odir/BMK_1_Global/BMK_2_DEG/$v");
		`cp $dir/$v/circos.png $odir/BMK_1_Global/BMK_2_DEG/$v/$v.Circos.png`;
		`cp $dir/$v/circos.svg $odir/BMK_1_Global/BMK_2_DEG/$v/$v.Circos.svg`;
		`cp $dir/$v/$v.pie.png $odir/BMK_1_Global/BMK_2_DEG/$v/$v.pie.png`;
		my @b = glob "$dir/$v/*txt";
		foreach my $b(@b){
			my $bb = basename $b;
			my $BB = (split /\./,$bb)[0];
			`cp $dir/$v/$BB.txt $odir/BMK_1_Global/BMK_2_DEG/$v/$BB.xls`;
		}
		`cp $dir/$v/*xls $odir/BMK_1_Global/BMK_2_DEG/$v`;
	}

$dir="$in/Combine/conservation";
`mkdir $odir/BMK_2_conservation`	unless(-d "$odir/BMK_2_conservation");
`cp $Bin/readme/BMK_7.3_conservation.pdf $odir/BMK_2_conservation/Readme.pdf`;
`cp $dir/*png $odir/BMK_2_conservation`;
`cp $dir/*pdf $odir/BMK_2_conservation`;
my @c = glob "$dir/*score";
foreach my $c(@c){
	my $cc = basename $c;
	`cp $c $odir/BMK_2_conservation/$cc.xls`;
}
$dir="$in/Combine/ceRNA";
`mkdir $odir/BMK_3_ceRNA_random` unless(-d "$odir/BMK_3_ceRNA_random");
`cp $Bin/readme/BMK_7.4_ceRNA.pdf $odir/BMK_3_ceRNA_random/Readme.pdf`;
my $ceRNA_file = (glob "$dir/*ceRNA_pair_adjust_p_Sig.txt")[0];
my $ceRNA_xls = basename $ceRNA_file;
$ceRNA_xls =~ s/\.txt$/.xls/;
`cp $ceRNA_file $odir/BMK_3_ceRNA_random/$ceRNA_xls`;

$dir = "$in/Combine/random";
my @edge_color_xls = glob "$in/Combine/random/*/*_edge.color.xls";
my @Key_gene_xls = glob "$in/Combine/random/*/*_Key_RNA_gene.xls";
my @node_xls = glob "$in/Combine/random/*/*_node.xls";
my @Kegg_pics = glob "$in/Combine/random/*/*KEGG_pathway_enrich_*";
my @KEGG_list = glob "$in/Combine/random/*/*KEGG_pathway_enrich.list";
foreach my $deg(@diff){
	`mkdir $odir/BMK_3_ceRNA_random/$deg` unless(-d "$odir/BMK_3_ceRNA_random/$deg");
	`mkdir $odir/BMK_3_ceRNA_random/$deg/BMK_1_KEGG` unless(-d "$odir/BMK_3_ceRNA_random/$deg/BMK_1_KEGG");
	my $random_dir = "$odir/BMK_3_ceRNA_random/$deg";
	`cp $dir/$deg/sub* $random_dir`;
	`cp $dir/$deg/*network.* $random_dir`;
	foreach my $a1(@edge_color_xls){
		if($a1=~/$deg/){`cp $a1 $random_dir`;}
	}
	foreach my $a2(@Key_gene_xls){
		if($a2=~/$deg/){`cp $a2 $random_dir`}
	}
	foreach my $a3(@node_xls){
		if($a3=~/$deg/){`cp $a3 $random_dir`;}
	}
	foreach my $a4(@Kegg_pics){
		if($a4=~/$deg/){`cp $a4 $random_dir/BMK_1_KEGG`;}
	}
	foreach my $a5(@KEGG_list){
		if($a5=~/$deg/){
	                my $KEGG_xls = basename $a5;
        	        $KEGG_xls =~ s/\.list$/.xls/;
                	`cp $a5 $random_dir/BMK_1_KEGG/$KEGG_xls`;
		}
        }
}


$dir = "$in/Combine/cytoscape";
`mkdir $odir/BMK_4_cytoscape` unless (-d "$odir/BMK_4_cytoscape");
`cp $Bin/readme/BMK_7.7_cytoscape.pdf $odir/BMK_4_cytoscape/Readme.pdf`;
my $cyto_sRNA = "$dir/sRNA";
my $cyto_lncRNA = "$dir/lncRNA";
my $cyto_circRNA = "$dir/circRNA";
my $cyto_gene = "$dir/gene";
my $lnc_sm = "$dir/lncRNA-miRNA-mRNA";
my $circ_sm = "$dir/circRNA-miRNA-mRNA";
`cp -r $cyto_sRNA $cyto_lncRNA $cyto_circRNA $cyto_gene $lnc_sm $circ_sm $odir/BMK_4_cytoscape/`;
`rm $odir/BMK_4_cytoscape/*/*r`;

#`cp -r $in/Combine/cytoscape/* $odir/BMK_4_cytoscape`;
#`rm $odir/BMK_4_cytoscape/*all.DEG.xls`;
#`rm $odir/BMK_4_cytoscape/*/*.r`;

$dir="$in/Combine/coexpression";
if(-d $dir){
	`mkdir $odir/BMK_5_coexp`	unless(-d "$odir/BMK_5_coexp");
	`cp $Bin/readme/BMK_7.6_coexp.pdf $odir/BMK_5_coexp/Readme.pdf`;
	`cp $dir/circRNA_miRNA/Sig_co_expression.txt $odir/BMK_5_coexp/Sig.circRNA_miRNA.xls`;
	`cp $dir/circRNA_miRNA/Diff_Sig_co_expression.txt $odir/BMK_5_coexp/Diff.Sig.circRNA_miRNA.xls`;
	`cp $dir/lncRNA_mRNA/Sig_co_expression.txt $odir/BMK_5_coexp/Sig.lncRNA_mRNA.xls`;
	`cp $dir/lncRNA_mRNA/Diff_Sig_co_expression.txt $odir/BMK_5_coexp/Diff.Sig.lncRNA_mRNA.xls`;
}
}

if(exists $steps{8}){
	if(exists $config{medical}){
		my $symbol;
		my $genome_dir = dirname $config{Ref_seq};
		my $list1 = "$genome_dir/id_name.list";
		my $list2 = "$genome_dir/lncRNA.id_name.list";
		if(exists $config{Symbol}){
			$symbol = $config{Symbol};
		}else{
		        if(-e $list1 && -e $list2){
                		&cmd("cat $list1 $list2 > $od/../work_sh/all.id_name.list");
				$symbol = "$od/../work_sh/all.id_name.list";
		        }elsif(-e $list1 && !-e $list2){
                		&cmd("cp $list1 $od/../work_sh/all.id_name.list");
				$symbol = "$od/../work_sh/all.id_name.list";
		        }else{
				print "You have no symbol file!\n";die;
			}
		}
		my @file=&getFiles;
		open(SH,">$od/../work_sh/8.1.add_symbol.sh")||die $!;
		###generate full_table.xls
		#wholetrans_mRNA_full_table.xls
		&cmd("perl $Bin/make_full_table.mRNA.pl -degdir $in/DEG_Analysis/gene -count $in/Basic_Analysis/Hisat_Stringtie/prepDE/All_gene_counts.list -fpkm $in/Basic_Analysis/Hisat_Stringtie/prepDE/All_gene_fpkm.list -anno $in/Anno_Integrate/Allgene_Anno/Result/Integrated_Function.annotation.xls -out $od/wholetrans_mRNA_full_table.xls -gnfile $symbol -cfg $detail_cfg");
		my $lnc_cmd = "perl $Bin/make_full_table.lncRNA.pl -degdir $in/DEG_Analysis/lncRNA -count $in/Basic_Analysis/Hisat_Stringtie/prepDE/All_lncRNA_counts.list -fpkm $in/Basic_Analysis/Hisat_Stringtie/prepDE/All_lncRNA_fpkm.list -anno $in/Anno_Integrate/Allgene_Anno/Result/Integrated_Function.annotation.xls -cis $in/Basic_Analysis/LncRNA_Analysis/Lnc_target_predict/Cis_target_gene.xls -out $od/wholetrans_lncRNA_full_table.xls -gnfile $symbol -cfg $detail_cfg";
		$lnc_cmd .= "-trans $Bin/Basic_Analysis/LncRNA_Analysis/Lnc_target_predict/Trans/Trans_target_gene.xls" if(-e "$Bin/Basic_Analysis/LncRNA_Analysis/Lnc_target_predict/Trans/Trans_target_gene.xls");
		&cmd("$lnc_cmd");
		&cmd("perl $Bin/make_full_table.circRNA.pl -degdir $in/DEG_Analysis/circRNA -count $in/Basic_Analysis/circRNA_analysis/expression/All_gene_counts.list -fpkm $in/Basic_Analysis/circRNA_analysis/expression/All_gene_expression.list -anno $in/Anno_Integrate/Allgene_Anno/Result/Integrated_Function.annotation.xls -host $in/Basic_Analysis/circRNA_analysis/expression/All_gene_counts_detail.xls -target $in/miRNA_Target/circRNA/mir2target.list -out $od/wholetrans_circRNA_full_table.xls -gnfile $symbol -cfg $detail_cfg -type $config{normalization}");
		&cmd("perl $Bin/make_full_table.sRNA.pl -degdir $in/DEG_Analysis/sRNA -count $in/Basic_Analysis/sRNA_Analysis/miRDeep2/miRNA_Quantify/All_miRNA.count.list -exp $in/Basic_Analysis/sRNA_Analysis/miRDeep2/miRNA_Quantify/All_miRNA_expression.list -anno $in/Anno_Integrate/Allgene_Anno/Result/Integrated_Function.annotation.xls -target $in/miRNA_Target/gene/mir2target.list -fa $in/Basic_Analysis/sRNA_Analysis/miRDeep2/miRNA_Quantify/All_miRNA.expressed.fa -out $od/wholetrans_sRNA_full_table.xls -gnfile $symbol -cfg $detail_cfg -type TPM");
		foreach my $f(@file){
			print SH "perl $Bin/convert_symbol.pl -id $symbol -i $f \n";
		}	
		close(SH);
		open (SH1,">$od/../work_sh/8.2.add_symbol.sh");
		print SH1 "perl $Bin/convert_symbol_other.pl -idir $od -name2list $symbol \n";
		close(SH1);
		`sh $od/../work_sh/8.1.add_symbol.sh > $od/../work_sh/8.1.add_symbol.sh.o`;
		`sh $od/../work_sh/8.2.add_symbol.sh > $od/../work_sh/8.2.add_symbol.sh.o`;	
		}
}

##########################################################
sub getFiles{
	my @files=();
	push @files,"$od/BMK_2_LncRNA/BMK_3_LncRNA_Expression/lncRNA_counts.xls";
	push @files,"$od/BMK_2_LncRNA/BMK_3_LncRNA_Expression/lncRNA_expression.xls";
	push @files,glob("$od/BMK_2_LncRNA/BMK_5_DEG_Analysis/BMK_1_All_DEG/All.DEG_final.xls");
	push @files,glob("$od/BMK_2_LncRNA/BMK_5_DEG_Analysis/*_vs_*/BMK_1_Statistics_Visualization/*.all.xls");
        push @files,glob("$od/BMK_2_LncRNA/BMK_5_DEG_Analysis/*_vs_*/BMK_1_Statistics_Visualization/*.DEG_final.xls");
	push @files,glob("$od/BMK_2_LncRNA/BMK_5_DEG_Analysis/*_vs_*/*Cis_Anno_enrichment/*enrich/*_vs_*_enrich.list.xls");
	push @files,glob("$od/BMK_2_LncRNA/BMK_5_DEG_Analysis/*_vs_*/*Trans_Anno_enrichment/*enrich/*_vs_*_enrich.list.xls");
	push @files,glob("$od/BMK_2_LncRNA/BMK_6_Lnc_disease/LncRNA_disease.xls");
	push @files,glob("$od/BMK_2_LncRNA/BMK_7_Lnc_precursor/lncRNA_precursor.xls");
	push @files,glob("$od/BMK_3_mRNA/BMK_2_geneExpression/gene_counts.xls");
	push @files,glob("$od/BMK_3_mRNA/BMK_2_geneExpression/gene_expression.xls");
	push @files,glob("$od/BMK_3_mRNA/BMK_2_geneExpression/Kmeans_Cluster/kmeans_cluster.xls");
	push @files,glob("$od/BMK_3_mRNA/BMK_3_DEG_Analysis/BMK_1_All_DEG/All.DEG_final.xls");
	push @files,glob("$od/BMK_3_mRNA/BMK_3_DEG_Analysis/BMK_1_All_DEG/All.DEG_final_anno.xls");
	push @files,glob("$od/BMK_3_mRNA/BMK_3_DEG_Analysis/*_vs_*/BMK_1_Statistics_Visualization/*_vs_*.all.xls");
	push @files,glob("$od/BMK_3_mRNA/BMK_3_DEG_Analysis/*_vs_*/BMK_1_Statistics_Visualization/*_vs_*.DEG_final.xls");
	push @files,glob("$od/BMK_3_mRNA/BMK_3_DEG_Analysis/*_vs_*/*Anno_enrichment/BMK_1_Annotation/*_vs_*.annotation.xls");
	push @files,glob("$od/BMK_3_mRNA/BMK_3_DEG_Analysis/*_vs_*/*Anno_enrichment/BMK_1_Annotation/*_vs_*.GO.list.xls");
	push @files,glob("$od/BMK_3_mRNA/BMK_3_DEG_Analysis/*_vs_*/*Anno_enrichment/*enrich/*_vs_*_enrich.list.xls");
	push @files,glob("$od/BMK_3_mRNA/BMK_3_DEG_Analysis/*_vs_*/*Anno_enrichment/*GSEA/*GSEA.xls");
	push @files,glob("$od/BMK_4_circRNA/BMK_1_circRNA_Prediction/BMK_1_CIRC_info/BMK_*_CIRI/*CIRI_Predict.xls");
	push @files,glob("$od/BMK_4_circRNA/BMK_1_circRNA_Prediction/BMK_1_CIRC_info/*intersect.xls");
	push @files,glob("$od/BMK_4_circRNA/BMK_3_DEG_Analysis/BMK_1_All_DEG/All.DEG_final_target.xls");
	push @files,glob("$od/BMK_4_circRNA/BMK_3_DEG_Analysis/BMK_1_All_DEG/All.DEG_final_anno.xls");
	push @files,glob("$od/BMK_4_circRNA/BMK_3_DEG_Analysis/*_vs_*/BMK_1_Statistics_Visualization/*_vs_*.DEG_final.Target.xls");
	push @files,glob("$od/BMK_4_circRNA/BMK_3_DEG_Analysis/*_vs_*/*Anno_enrichment/BMK_1_Annotation/*_vs_*.annotation.xls");
	push @files,glob("$od/BMK_4_circRNA/BMK_3_DEG_Analysis/*_vs_*/*Anno_enrichment/*enrich/*_vs_*_enrich.list.xls");
	push @files,glob("$od/BMK_4_circRNA/BMK_1_circRNA_Prediction/overlap_alitisplice.xls");
	push @files,glob("$od/BMK_*_mRNA/BMK_3_DEG_Analysis/*/TFs_activity_grn.xls");
        push @files,glob("$od/BMK_*_mRNA/BMK_3_DEG_Analysis/*/each_geneRes/*txt");
        push @files,glob("$od/BMK_*_mRNA/BMK_3_DEG_Analysis/*/newGenes_TFBS_predictRes.xls");
	push @files,glob("$od/BMK_*_mRNA/BMK_3_DEG_Analysis/*_vs_*/BMK_3_diff_AS_analysis/*.JC.xls");
#	push @files,glob("$od/BMK_4_circRNA/BMK_1_circRNA_Prediction/circRNA_newname.xls");
	push @files,glob("$od/BMK_5_miRNA/BMK_3_DEG_Analysis/BMK_1_All_DEG/All.DEG_final_target.xls");
	push @files,glob("$od/BMK_5_miRNA/BMK_3_DEG_Analysis/BMK_1_All_DEG/All.DEG_final_anno.xls");
	push @files,glob("$od/BMK_5_miRNA/BMK_3_DEG_Analysis/*_vs_*/BMK_1_Statistics_Visualization/*_vs_*.DEG_final.Target.xls");
	push @files,glob("$od/BMK_5_miRNA/BMK_3_DEG_Analysis/*_vs_*/*Anno_enrichment/BMK_1_Annotation/*_vs_*.annotation.xls");
	push @files,glob("$od/BMK_5_miRNA/BMK_3_DEG_Analysis/*_vs_*/*Anno_enrichment/*enrich/*_vs_*_enrich.list.xls");
	push @files,glob("$od/BMK_6_Structure/BMK_1_SNP_Analysis/final.*.anno.gatk.all.list.xls");
	push @files,glob("$od/BMK_6_Structure/BMK_2_miRNA_Target/gene.mir2target.list.xls");
	push @files,glob("$od/BMK_6_Structure/BMK_2_miRNA_Target/gene.miRNA_target_info.xls");
	push @files,glob("$od/BMK_6_Structure/BMK_2_miRNA_Target/lncRNA.mir2target.list.xls");
	push @files,glob("$od/BMK_6_Structure/BMK_2_miRNA_Target/lncRNA.miRNA_target_info.xls");
	push @files,glob("$od/BMK_6_Structure/BMK_3_Gene_Structure_Optimize/*.geneStructure.optimize.xls");
	push @files,glob("$od/BMK_6_Structure/BMK_4_Alt_splice/*.AS.list.xls");
	push @files,glob("$od/BMK_7_Combine/BMK_2_conservation/lncRNA.PhastCons.score.xls");
	push @files,glob("$od/BMK_7_Combine/BMK_2_conservation/gene.PhastCons.score.xls");
	push @files,glob("$od/BMK_7_Combine/BMK_3_ceRNA_random/*/BMK_1_KEGG/*enrich.list.xls");
	return(@files);
}
sub getDensity{
	my ($density_dir,$png_dir)=@_;
	`cp $density_dir/*png $png_dir/`;
	`cp $density_dir/*pdf $png_dir/`;
	`cp $density_dir/sample_coefficient.txt $png_dir/sample_coefficient.xls`;
}


sub getDEG{
	my ($dir,$odir,$type)=@_;
	my $all_deg="$dir/All_DEG/All.DEG_final.xls";	
	if(-e $all_deg){
		opendir(DIR,$dir)||die $!;
		my @vs=grep{/_vs_/ && -d "$dir/$_"}readdir(DIR);
		closedir(DIR);
		`mkdir $odir`                unless(-d $odir);
		`rm -r $odir/BMK_1_All_DEG` if(-d "$odir/BMK_1_All_DEG");
		`cp -r $dir/All_DEG $odir/BMK_1_All_DEG`;
		my $s = 1;
		###PPI
		my @ppi=glob("$dir/DEG_PPI/*_vs_*sif");
		if(@ppi>0){
			$s++;
			`mkdir $odir/BMK_$s\_DEG_PPI` unless(-d "$odir/BMK_$s\_DEG_PPI");
			`cp $dir/DEG_PPI/*_vs_*sif $odir/BMK_$s\_DEG_PPI/`;
			my @txt = glob "$dir/DEG_PPI/*_vs_*txt";
			foreach my $TXT(@txt){
				my $a = basename $TXT;
				$a =~s/\.txt//;
				`cp $TXT $odir/BMK_$s\_DEG_PPI/$a.xls`;
			}
		}
		###TFBS analysis
		if(-d "$dir/TFBS_Analysis"){
			$s++;
			`mkdir $odir/BMK_$s\_TFBS_Analysis` unless(-d "$odir/BMK_$s\_TFBS_Analysis");
			`mkdir $odir/BMK_$s\_TFBS_Analysis/each_DEgeneRes` unless(-d "$odir/BMK_$s\_TFBS_Analysis/each_DEgeneRes");
			if(-e "$dir/TFBS_Analysis/allGenes_TFBS_predictRes.txt"){
				`cp $dir/TFBS_Analysis/allGenes_TFBS_predictRes.txt $odir/BMK_$s\_TFBS_Analysis/allGenes_TFBS_predictRes.xls`;
			}elsif(-e "$dir/TFBS_Analysis/newGenes_TFBS_predictRes.txt"){
				`cp $dir/TFBS_Analysis/newGenes_TFBS_predictRes.txt $odir/BMK_$s\_TFBS_Analysis/newGenes_TFBS_predictRes.xls`;
			}
			`mkdir $odir/BMK_$s\_TFBS_Analysis/DEG_seqLogo` unless(-d "$odir/BMK_$s\_TFBS_Analysis/DEG_seqLogo");
			`cp $dir/TFBS_Analysis/DEG_seqLogo/*.p* $odir/BMK_$s\_TFBS_Analysis/DEG_seqLogo`;
			`cp $dir/TFBS_Analysis/each_DEgeneRes/* $odir/BMK_$s\_TFBS_Analysis/each_DEgeneRes`;
			`cp $Bin/readme/BMK_3.7_TFBS_analysis_readme.pdf $odir/BMK_$s\_TFBS_Analysis/Readme.pdf`;
		}
		###TF_activity for diff gene
		if(-d "$dir/TF_activity"){
                        $s++;
                        `mkdir $odir/BMK_$s\_TF_activity` unless (-d "$odir/BMK_$s\_TF_activity");
                        `cp $dir/TF_activity/TFs_network.p* $odir/BMK_$s\_TF_activity/`;
                        `cp $dir/TF_activity/TFs_influences_heatmap.p* $odir/BMK_$s\_TF_activity/`;
                        `cp $dir/TF_activity/TFs_cornet.txt $odir/BMK_$s\_TF_activity/TFs_cornet.xls`;
                        `cp $dir/TF_activity/TFs_activity_grn.txt $odir/BMK_$s\_TF_activity/TFs_activity_grn.xls`;
                        `cp $dir/TF_activity/TFs_influence.txt $odir/BMK_$s\_TF_activity/TFs_influence.xls`;
			`cp $Bin/readme/BMK_3.8_TF_activity_readme.pdf $odir/BMK_$s\_TF_activity/Readme.pdf`;
                }

		$s++;
		open(STAT,">$odir/DEG.stat.xls")||die $!;
		print STAT "#DEG Set\tAll DEG\tup-regulated\tdown-regulated\n";
		for(my $i=0;$i<@vs;$i++){
        		my $id=$i+$s;
		        my $vdir="$odir/BMK_${id}_$vs[$i]";
			`rm -r $vdir` if(-d $vdir);
        		`mkdir $vdir` unless(-d $vdir);
			`mkdir $vdir/BMK_1_Statistics_Visualization` unless (-d "$vdir/BMK_1_Statistics_Visualization");
	        	`cp $dir/$vs[$i]/*.xls $vdir/BMK_1_Statistics_Visualization && cp $dir/$vs[$i]/$vs[$i].all $vdir/BMK_1_Statistics_Visualization/$vs[$i].all.xls`;
		        `cp $dir/$vs[$i]/*.png $vdir/BMK_1_Statistics_Visualization && cp $dir/$vs[$i]/*.pdf $vdir/BMK_1_Statistics_Visualization`;

		        `cp -r $dir/$vs[$i]/*Anno_enrichment $vdir`;
			`rm -rf $dir/$vs[$i]/*Anno_enrichment/pathway/kegg_enrichment`;
			my @enrich_dir = glob "$vdir/*Anno_enrichment";
			my $j=1;
			foreach my $dir_enrich(@enrich_dir){
				$j++;
				my $dirname = dirname $dir_enrich;
				my $newname = (split /\//,$dir_enrich)[-1];
				`rm -r $dirname/BMK_${j}_$newname` if(-d "$dirname/BMK_${j}_$newname");
				`mv $dir_enrich $dirname/BMK_${j}_$newname`;
				$dir_enrich = "$dirname/BMK_${j}_$newname";
				`rm $dir_enrich/anno/*.GO.list.txt`;
				`mv $dir_enrich/anno $dir_enrich/BMK_1_Annotation`;
				`mkdir $dir_enrich/BMK_2_GO_enrich` unless (-d "$dir_enrich/BMK_2_GO_enrich");
				`mkdir $dir_enrich/BMK_3_KEGG_enrich` unless (-d "$dir_enrich/BMK_3_KEGG_enrich");
				if(-d "$dir_enrich/enrich/GSEA"){
					`mkdir $dir_enrich/BMK_4_GSEA` unless (-d "$dir_enrich/BMK_4_GSEA");
					`mv $dir_enrich/enrich/GSEA $dir_enrich/BMK_4_GSEA/BMK_1_Plot`;
					`mv $dir_enrich/enrich/*GSEA.xls $dir_enrich/BMK_4_GSEA/`;
					`mv $dir_enrich/enrich/*Phase* $dir_enrich/BMK_3_KEGG_enrich`;
				}
				`mv $dir_enrich/enrich/*Biological_Process_enrich*.p* $dir_enrich/BMK_2_GO_enrich`;
				`mv $dir_enrich/enrich/*Cellular_Component_enrich*.p* $dir_enrich/BMK_2_GO_enrich`;
				`mv $dir_enrich/enrich/*Molecular_Function_enrich*.p* $dir_enrich/BMK_2_GO_enrich`;
				my @go_list = glob "$dir_enrich/enrich/*Biological_Process_enrich.list";
				push @go_list,glob "$dir_enrich/enrich/*Cellular_Component_enrich.list";
				push @go_list,glob "$dir_enrich/enrich/*Molecular_Function_enrich.list";
				foreach my $golist (@go_list){
					my $file = basename $golist;
					`mv $golist $dir_enrich/BMK_2_GO_enrich/$file.xls`;
				}
				`mv $dir_enrich/enrich/*KEGG_pathway_enrich*.p* $dir_enrich/BMK_3_KEGG_enrich`;
				my @kegg_list = glob "$dir_enrich/enrich/*KEGG_pathway_enrich.list";
				foreach my $kegglist(@kegg_list){
					my $file = basename $kegglist;
					`mv $kegglist $dir_enrich/BMK_3_KEGG_enrich/$file.xls`;
				}
				`mv $dir_enrich/pathway/kegg_map $dir_enrich/BMK_3_KEGG_enrich/BMK_1_Keggmap`;
				`rm -r $dir_enrich/pathway`;
				`rm -r $dir_enrich/enrich`;
			}
			my $all=`wc -l $dir/$vs[$i]/$vs[$i].DEG_final.xls`;
			chomp($all);$all=(split(/\s+/,$all))[0]-1;
			my $up=`awk '{print \$NF}' $dir/$vs[$i]/$vs[$i].DEG_final.xls|grep up|wc -l`;
			chomp($up);$up=(split(/\s+/,$up))[0];
			my $down=$all-$up;
			print STAT "$vs[$i]\t$all\t$up\t$down\n";
			if($type eq "gene"){
				###rMATs diff AS
				my @diff_AS = glob "$in/differential_AS_analysis/*_vs_*";
				if(@diff_AS>0){
					my $vs = $vs[$i];
					foreach my $diffas (@diff_AS){
						if($diffas =~/$vs/){
							`mkdir $vdir/BMK_3_diff_AS_analysis` unless (-d "$vdir/BMK_3_diff_AS_analysis");
							`sed -i 's/\"//g' $diffas/*.JC.txt `;
							`cut -f2,4-23 $diffas/A3SS.MATS.JC.txt > $vdir/BMK_3_diff_AS_analysis/A3SS.MATS.JC.xls`;
							`cut -f2,4-23 $diffas/A5SS.MATS.JC.txt > $vdir/BMK_3_diff_AS_analysis/A5SS.MATS.JC.xls`;
							`cut -f2,4-23 $diffas/MXE.MATS.JC.txt > $vdir/BMK_3_diff_AS_analysis/MXE.MATS.JC.xls`;
							`cut -f2,4-23 $diffas/RI.MATS.JC.txt > $vdir/BMK_3_diff_AS_analysis/RI.MATS.JC.xls`;
							`cut -f2,4-23 $diffas/SE.MATS.JC.txt > $vdir/BMK_3_diff_AS_analysis/SE.MATS.JC.xls`;
							`cp $Bin/readme/BMK_3.6_DEGanalysis_rMATS.pdf $vdir/BMK_3_diff_AS_analysis/Readme.pdf`;
						}
					}
					&getDiff_AS("$in/differential_AS_analysis","$odir/Diff.AS.stat.xls");
				}
			}
		}
		close(STAT);
	}
}



sub DEG_anno_stat{
	my ($dir,$path,$out)=@_;
#L01_vs_L11/Anno_enrichment/anno/L01_vs_L11.annotation.xls
	my @annos=("COG_class_annotation","GO_annotation","KEGG_annotation","KOG_class_annotation","Swiss-Prot_annotation","eggNOG_class_annotation","NR_annotation");
	my @files=glob("$dir/*_vs_*/$path/anno/*.annotation.xls");
	if(@files>0){
		open(OUT,">$out")||die $!;
		print OUT "#DEG Set\tAnnotated\t",join("\t",map{(split(/_annotation/,$_))[0]}@annos),"\n";
		foreach my $file(@files){
			my %info=();my %genes=();
			open(ANNO,$file)||die $!;
			my $head=<ANNO>;chomp($head);my @header=split(/\t/,$head);
			my %index=();
			for(my $i=0;$i<@header;$i++){
				foreach my $a(@annos){
					$index{$a}=$i if($header[$i] eq $a);
				}
			}
			while(<ANNO>){
				chomp;next if($_=~/^#/);
				my @tmp=split(/\t/,$_);
				for(my $i=0;$i<@annos;$i++){
					if($tmp[$index{$annos[$i]}] ne "--"){
						$info{$annos[$i]}{$tmp[0]}++;
						$genes{$tmp[0]}++;
					}	
				}		
			}
			close(ANNO);
			my @values=();
			for(my $i=0;$i<@annos;$i++){
				push @values, scalar(keys %{$info{$annos[$i]}});
			}
			my $base=basename $file;$base=(split(/\./,$base))[0];
			print OUT "$base\t",scalar(keys %genes),"\t",join("\t",@values),"\n";
		}
		close(OUT);
	}
}

sub getDiff_AS{
	my ($dir,$o)=@_;
	opendir(DIR,$dir)||die $!;
	my @vs=grep{/_vs_/ && -d "$dir/$_"}readdir(DIR);
	closedir(DIR);
	open(STAT,">$o")||die $!;
	print STAT "DEG Set\tA3SS\tA5SS\tMXE\tRI\tSE\n";
	for(my $i=0;$i<@vs;$i++){
		my $A3SS=`wc -l $dir/$vs[$i]/A3SS.MATS.JC.txt`;
		$A3SS=(split(/\s+/,$A3SS))[0]-1;
		my $A5SS=`wc -l $dir/$vs[$i]/A5SS.MATS.JC.txt`;
		$A5SS=(split(/\s+/,$A5SS))[0]-1;
		my $MXE=`wc -l $dir/$vs[$i]/MXE.MATS.JC.txt`;
		$MXE=(split(/\s+/,$MXE))[0]-1;
		my $RI=`wc -l $dir/$vs[$i]/RI.MATS.JC.txt`;
		$RI=(split(/\s+/,$RI))[0]-1;
		my $SE=`wc -l $dir/$vs[$i]/SE.MATS.JC.txt`;
		$SE=(split(/\s+/,$SE))[0]-1;
		print STAT "$vs[$i]\t$A3SS\t$A5SS\t$MXE\t$RI\t$SE\n";
	}
	close(STAT);
}

######################################################################
#			Sub function
######################################################################
##############################basic sub function
sub read_Config{
	my $cfg = shift; 
	open(IN,$detail_cfg) or die $!;
	while(<IN>){
		chomp;
		next if(/^#|^$|^\s$/);
		my @tmp = split /\s+/,$_;
		if($tmp[0] eq "Sep"){
			$tmp[1]=~ s/;/_vs_/g;
			$tmp[1]=~ s/,/_/g;
			push @diff,$tmp[1];
		}elsif($tmp[1] eq "Com"){
			$tmp[1]=~ s/,/_vs_/g;
			push @diff,$tmp[1];
		}else{
			$config{$tmp[0]}=$tmp[1];
		}
	}
	close(IN);
#	my $d=Config::General->new(-ConfigFile => "$configFile");
#	my %config=$d->getall;	
#	return %config;
}

sub qsub(){
        my $shfile= shift;
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $shfile";
        &run_or_die($cmd);              
        return ;
}
sub run_or_die(){
        my ($cmd) = @_ ;
        &show_log($cmd);
        my $flag = system($cmd) ;
        if ($flag != 0){
                &show_log("Error: command fail: $cmd");
                exit(1);
        }
        &show_log("done.");
        return ;
}
sub show_log(){
        my ($txt) = @_ ;
        my $time = time();
        my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime($time);
        $wday = $yday = $isdst = 0;
        my $Time=sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
        print "$Time:\t$txt\n" ;
        return ($time) ;
}
sub cmd{
	my $shell = shift;
	print "$shell\n";
	my $flag = system("$shell");
	if($flag ==0 ){
		print "success!\n";
	}else{
		print "fail,please Check,good luck!\n";die;
	}
}


sub USAGE{
	my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
	The config file for TGS data and SGS data in $Bin/Config can be refered!

	Options:
	-in	path	analysis output path
	-od 	path	analysis report path, default in/BMK_Result
	-step 	str	default, 1,2,3,4,5,6,7,8
			1. raw data assess, map
			2. lncRNA(predict,exp,target,DEG)
			3. mRNA(new_anno,exp,DEG)
			4. circRNA(predict,exp,DEG)
			5. miRNA(predict,exp,DEG,Family,edit)
			6. Structure (SNP,Fusion,AS,Optimize,miRNA_Target)
			7. Combine (Global, Global_diff, ceRNA, coexp)
			8. generate *full_table.xls and add symbol
				
	-cfg1	file	data.cfg
	-cfg2	file	detail.cfg
			The two cfg file are only needed by the step 8.
	-h	Help

Example: perl $0 -in Analysis -od Analysis/BMK_Result -cfg1 data.cfg -cfg2 detail.cfg -step 1,2,3,4,5,6,7,8 

USAGE
	print $usage;
	exit;
}


