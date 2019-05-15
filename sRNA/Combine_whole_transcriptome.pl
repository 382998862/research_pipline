use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path getcwd);
use threads;
use newPerlBase;
my ($data_cfg, $detail_cfg, $od, $step,$switch,$exosome);
GetOptions(
        "h|?"           =>\&USAGE,
	"cfg1:s"	=>\$data_cfg,
	"cfg2:s"	=>\$detail_cfg,
	"o:s"     	=>\$od,
	"switch:s"	=>\$switch,
	"step:s"	=>\$step,
	"exosome"	=>\$exosome,
)or &USAGE;
&USAGE unless ($data_cfg and $detail_cfg and $od);

$switch||="on";
#################调用脚本
my %sys=&readConfig("$Bin/script.cfg");
my $python		=$sys{python};
my $Rscript		=$sys{Rscript};
my $lncRNA_data_assess	=$sys{lncRNA_data_assess};

my $gene_anno		="$Bin/lncRNA/annotation/Gene_Anno/Gene_Func_Anno_Pipline.pl";
my $anno_integrated	="$Bin/lncRNA/annotation/annotation_integrate/anno_integrated2.pl";
my $PPI			="$Bin/tools/ppi/PPI.py";
my $miRNA_data_assess   ="$Bin/sRNA/data_assess/data_assess.pl";
my $miRNA_map		="$Bin/sRNA/alignment/alignment.pl";
my $miRNA_identify	="$Bin/sRNA/basic_analysis/miRNA_indentify.pl";
my $miRNA_target	="$Bin/sRNA/target/miRNA_target.pl";
my $rMATs		="$Bin/differential_AS_analysis/diff_AS_analysis.pl";
my $miRNA_family	=$sys{miRNA_family};
my $miRNA_edit		="$Bin/sRNA/edit/miRNA_Edit.pl";
my $HMDD		="/share/nas2/genome/bmksoft/pipeline/sRNA_process/v2.5/bin/HMDD_Annotation/HMDD_annotation.pl";

my $circRNA_bwa		="$Bin/circRNA/bwa/bwa_process.pl";
my $circRNA_map		=$sys{circRNA_map};
my $circRNA_pre		="$Bin/circRNA/circRNA_identify/circRNA_analysis.pl";
my $count2len 		="$Bin/circRNA/circRNA_identify/bin/count_to_length.pl";

my $DEG_Analysis	="$Bin/DEG_Analysis/DEG_Analysis.pl";
my $gene_fusion		="$Bin/gene_fusion/fusionmap3.pl";
my $optimize		="/share/nas2/genome/bmksoft/pipeline/LncRNA_pipeline/v3.1.6/bin/gene_optimize/gene_structure_optimize.pl";
my $altsplice		=$sys{altsplice};
my $SNP			=$sys{SNP};
#################
$step||=join(",",1..8);
my %steps=();
foreach my $s (split(/,/,$step)){
	$steps{$s}++;
}

`mkdir -p $od`	unless(-d $od);
`mkdir -p $od/Config`  unless(-d "$od/Config");

$od		=abs_path($od);
$data_cfg	=abs_path($data_cfg);
$detail_cfg	=abs_path($detail_cfg);


my %config=();
&readcfg($detail_cfg);
my %sample=&relate($data_cfg);
my $index		=$config{Project_key};
my $known_unigene	=$config{Known_unigene};
my $known_anno		=$config{Known_anno};

my %sRNA_samples=();
my $miRNA_cfg=&getMiRNAcfg($data_cfg,$detail_cfg);

my $work_sh="$od/work_sh";
`mkdir $work_sh`		unless(-d $work_sh);
`mkdir $work_sh/lncRNA`		unless(-d "$work_sh/lncRNA");
`mkdir $work_sh/circRNA`	unless(-d "$work_sh/circRNA");
`mkdir $work_sh/miRNA`		unless(-d "$work_sh/miRNA");
my ($odir,$lines,$cmd,$miRNA_line);

my $qsubcmd="sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent --queue $config{Queue_type} ";
my $qsubcmd_once="sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --independent --queue $config{Queue_type} ";
$odir="$od/Basic_Analysis";
my($lnc_fpkm,$lnc_count,$circ_exp,$circ_count,$gene_fpkm,$gene_count,$miRNA_exp,$miRNA_count);
$gene_fpkm      ="$odir/Hisat_Stringtie/prepDE/All_gene_fpkm.list";
$gene_count     ="$odir/Hisat_Stringtie/prepDE/All_gene_counts.list";
$circ_exp       ="$odir/circRNA_analysis/expression/All_gene_expression.list";
$circ_count     ="$odir/circRNA_analysis/expression/All_gene_counts.list";
$lnc_fpkm       ="$odir/Hisat_Stringtie/prepDE/All_lncRNA_fpkm.list";
$lnc_count      ="$odir/Hisat_Stringtie/prepDE/All_lncRNA_counts.list";
$miRNA_exp      ="$odir/sRNA_Analysis/miRDeep2/miRNA_Quantify/All_miRNA_expression.list";
$miRNA_count    ="$odir/sRNA_Analysis/miRDeep2/miRNA_Quantify/All_miRNA.count.list";

my($miRNA_fa,$circRNA_fa,$lncRNA_fa,$mRNA_fa,$novel_unigene);
$miRNA_fa    ="$odir/sRNA_Analysis/miRDeep2/miRNA_Quantify/All_miRNA.expressed.fa";
$circRNA_fa  ="$odir/circRNA_analysis/circRNA.fa";
$lncRNA_fa   ="$odir/LncRNA_Analysis/lnc_predict/All_filter_final.fa";
$mRNA_fa     ="$odir/Hisat_Stringtie/genePredict/All.longest_transcript.fa";

my($new_gene_gff,$gene_gff,$lnc_gff);
$new_gene_gff	="$odir/Hisat_Stringtie/genePredict/final_track/$index.newGene_final.filtered.gff";
$gene_gff	="$odir/Hisat_Stringtie/genePredict/All_gene_final.gff";
$lnc_gff	="$odir/LncRNA_Analysis/lnc_predict/All_filter_final.gff";

my $anno	="$od/Anno_Integrate/Allgene_Anno/Result";

############data_assess
if(exists $steps{1}){
	$odir="$od/Data_Assess";
	&writeSH("perl $lncRNA_data_assess -config $data_cfg -detailfig $detail_cfg -outdir $odir/lncRNA ","$work_sh/lncRNA/s1.data_assess.sh");
	&writeSH("perl $miRNA_data_assess -Q $config{Q} -main_cfg $data_cfg -word_cfg $detail_cfg -od $odir/miRNA -min $config{min_len} -max $config{max_len}","$work_sh/miRNA/s1.data_assess.sh");
	$lines ="$work_sh/lncRNA/s1.data_assess.sh";
	&run($lines,"$work_sh/s1.data_assess.sh");
	&qsubCheck("$work_sh/s1.data_assess.sh");
}

############Basic_Analysis
if(exists $steps{2}){
	$odir="$od/Basic_Analysis";
	`mkdir $odir` unless (-d $odir);
	#LncRNA_mapping and miRNA data_assess if $switch eq "on"
	&writeSH("perl $Bin/lncRNA/hisat_string/hisat2_stringtie.pl -cfg1 $data_cfg -cfg2 $detail_cfg -od $odir/Hisat_Stringtie","$work_sh/lncRNA/s2.1.lnc_map.sh");
	if($switch eq "on"){
		$lines ="$work_sh/miRNA/s1.data_assess.sh";
		$lines .=",$work_sh/lncRNA/s2.1.lnc_map.sh";
	}else{
		$lines = "$work_sh/lncRNA/s2.1.lnc_map.sh";
	}
	&run($lines,"$work_sh/s2.1.lncRNA_map_sRNA_data_assess.sh");
	&qsubCheck("$work_sh/s2.1.lncRNA_map_sRNA_data_assess.sh");
	#LncRNA_map stat
	&writeSH("perl $Bin/lncRNA/hisat_string/map_stat/Statistics_draw_hisat.pl -cfg2 $detail_cfg -in $odir/Hisat_Stringtie -od $odir/Hisat_Stringtie/Mapped ","$work_sh/lncRNA/s2.2.map_stat.sh");

	#circRNA_mapping
	`mkdir $odir/circRNA_Bwa`	unless(-d "$odir/circRNA_Bwa");
	`echo "Ref_seq\t$config{Ref_seq}\nRef_ann\t$gene_gff" >$odir/circRNA_Bwa/detail.cfg`;
	`echo "medical\t$config{medical}" >>$odir/circRNA_Bwa/detail.cfg`	if(exists $config{medical});
	&writeSH("perl $circRNA_bwa -cfg1 $data_cfg -cfg2 $odir/circRNA_Bwa/detail.cfg -od $odir/circRNA_Bwa && perl $circRNA_map -id $odir/circRNA_Bwa/Map_Stat/ -o $odir/circRNA_Bwa/Map_Stat/All.mappedStat.xls && perl $count2len -cfg $data_cfg -o $odir/circRNA_Bwa/Map_Stat/All.readLength.xls","$work_sh/circRNA/s2.map.sh");

	#sRNA_mapping
	`mkdir $odir/sRNA_Alignment`	unless(-d "$odir/sRNA_Alignment");
	&writeSH("perl $miRNA_map  -od $odir/sRNA_Alignment -id $od/Data_Assess/miRNA -main_cfg $detail_cfg -data_cfg $data_cfg -filter 1 -genome $config{Ref_seq} -gff $gene_gff ","$work_sh/miRNA/s2.alignment.sh");

	$lines ="$work_sh/circRNA/s2.map.sh,$work_sh/lncRNA/s2.2.map_stat.sh";
	$lines.=",$work_sh/miRNA/s2.alignment.sh"	if($switch eq "on");
	&run($lines,"$work_sh/s2.2.circRNA_sRNA_map.sh");
	&qsubCheck("$work_sh/s2.2.circRNA_sRNA_map.sh");
}

$lines = "";
if(exists $steps{3}){
	###lncRNA identify
	&writeSH("perl $Bin/lncRNA/lncRNA_Analysis/Lncrna_analysis.v3.pl -cfg $detail_cfg -gtf $odir/Hisat_Stringtie/LncPredict/merged_filter.gtf -gff $gene_gff -od $odir/LncRNA_Analysis -prep $odir/Hisat_Stringtie/prepDE -type $config{type} -db $config{db}","$work_sh/lncRNA/s3.identify_target.sh");

	###circRNA identify
	&writeSH("perl $circRNA_pre -data_cfg $data_cfg -detail_cfg $detail_cfg -od $odir/circRNA_analysis -gff $gene_gff ","$work_sh/circRNA/s3.identify.sh");

	##miRNA identify
	$config{miRBase}=$config{Species_type}  if(!exists $config{miRBase});
	$cmd="perl $miRNA_identify -mature /share/nas2/database/miRBase/v21/mature.fa -miR $config{miRBase} -idir $odir/sRNA_Alignment -od $odir/sRNA_Analysis -maincfg $miRNA_cfg -min $config{min_len} -max $config{max_len} ";
	if(exists $config{medical}){
			$cmd.="&& perl $HMDD -in $odir/sRNA_Analysis/miRDeep2/miRNA_Quantify/All_miRNA_expression.list -od $odir/sRNA_Analysis/HMDD"	if($config{medical} eq "GRCh37" || $config{medical} eq "GRCh38");
	}
	&writeSH($cmd,"$work_sh/miRNA/s3.identify.sh");

	$lines ="$work_sh/lncRNA/s3.identify_target.sh,$work_sh/circRNA/s3.identify.sh";
	$miRNA_line ="$work_sh/miRNA/s3.identify.sh";

	if(!exists $steps{4}){
		$lines.="$work_sh/miRNA/s3.identify.sh"                 if($switch eq "on");
		&run($lines,"$work_sh/s3_4.identify_structure.sh");
		&qsubCheck("$work_sh/s3_4.identify_structure.sh");
	}
}

if(exists $steps{4}){
	############Independent Structure Analysis
	###Function analysis
	#`mkdir -p $od/Anno_Integrate/New_Anno`  unless(-d "$od/Anno_Integrate/New_Anno");
	$novel_unigene="$odir/Hisat_Stringtie/genePredict/final_track/$index.newGene.longest_transcript.fa";
	my $anno_cfg="$od/Config/anno.cfg"; &annoCfg($anno_cfg);

	my $cmd	= "perl $gene_anno  --cfg $anno_cfg --od $od/Anno_Integrate/New_Anno --queue $config{Queue_type} --eggNOG --pfam --kog ";
	$cmd.= "--nr "          if (exists $config{nr});
	$cmd.= "--swissprot "   if (exists $config{Swissprot});
	$cmd.= "--kegg "        if (exists $config{Kegg});
	$cmd.= "--GO "          if (exists $config{nr});
	$cmd.= "--cog "         if (exists $config{Cog});
	$cmd.= "&& perl $anno_integrated -i ${index}_Unigene -gene $config{Known_unigene},$novel_unigene -anno $config{Known_anno},$od/Anno_Integrate/New_Anno -cfg $detail_cfg -od $od/Anno_Integrate/Allgene_Anno";
	&writeSH($cmd,"$work_sh/lncRNA/s4_Anno_Integrate.sh");

	###SNP
	&writeSH("perl $SNP -cfg $detail_cfg -hisat $odir/Hisat_Stringtie/Hisat -gff $new_gene_gff -od $od/SNP_Analysis","$work_sh/lncRNA/s5.SNP.sh");

	###AS
	&writeSH("perl $altsplice -cfg $detail_cfg -cufflinks $odir/Hisat_Stringtie/StringTie -Ref_Genome $odir/Hisat_Stringtie/Ref_Genome -od $od/Alitsplice_Analysis","$work_sh/lncRNA/s6.AS.sh");

	###Gene_Structure_Optimize
	&writeSH("perl $optimize -gtf $odir/Hisat_Stringtie/Compare/gffcmp.annotated.gtf -gff $config{Ref_ann} -od $od/Gene_Structure_Optimize -index $index","$work_sh/lncRNA/s7.Gene_Structure_Optimize.sh");
	###PPI
	&writeSH("$python $PPI -fa $novel_unigene -cfg $detail_cfg -of $od/PPI/PPI.txt -q medical.q -p $config{Project_key}","$work_sh/lncRNA/s8.newGene_PPI.sh");

	###rMATs diff AS
	my $rMATs_cfg = "$od/Config/rMATs.cfg";&rMATsCfg($rMATs_cfg);
	&writeSH("perl $rMATs -cfg $rMATs_cfg -bamdir $odir/Hisat_Stringtie/Hisat -od $od/differential_AS_analysis","$work_sh/lncRNA/s9.different_AS.sh");

	$lines.=",$work_sh/lncRNA/s4_Anno_Integrate.sh,$work_sh/lncRNA/s5.SNP.sh,$work_sh/lncRNA/s6.AS.sh,$work_sh/lncRNA/s7.Gene_Structure_Optimize.sh,$work_sh/lncRNA/s8.newGene_PPI.sh,$work_sh/lncRNA/s9.different_AS.sh";

	###Gene fusion
	if(exists $config{medical}){
		my $fusion_db=&get_FusionType($config{medical});
		&writeSH("perl $gene_fusion --indir $odir/Hisat_Stringtie/Hisat --type $fusion_db --od $od/Gene_Fusion --script_cfg $Bin/gene_fusion/script.cfg --cfg2 $detail_cfg","$work_sh/lncRNA/s10.fusion.sh");
		$lines.=",$work_sh/lncRNA/s10.fusion.sh";
	}
	###miRNA edit family
	my $sam=join(",",keys %sRNA_samples);
	$cmd ="perl $miRNA_family -i $miRNA_fa  -type $config{Species_type} -od $od/miRNA_Family ";
	$cmd.="&& perl $miRNA_edit -idir $odir -od $od/miRNA_edit -sample $sam ";
	&writeSH("perl $miRNA_family -i $miRNA_fa  -type $config{Species_type} -od $od/miRNA_Family && perl $miRNA_edit -idir $odir -od $od/miRNA_edit -sample $sam ","$work_sh/miRNA/s4.miRNA_family_edit.sh");

	$miRNA_line .= " && $work_sh/miRNA/s4.miRNA_family_edit.sh";
	$miRNA_line  =~s/^&&//;

	$lines.=",$miRNA_line"	if($switch eq "on");
	$lines=~s/^,//;
	&run($lines,"$work_sh/s3_4.identify_structure.sh");
	&qsubCheck("$work_sh/s3_4.identify_structure.sh");
}


if(exists $steps{5}){
	############miRNA target prediction
	if(-e $miRNA_fa && -e $mRNA_fa && -e $circRNA_fa && -e $lncRNA_fa){
		$cmd  ="perl $miRNA_target -miRNA $miRNA_fa -target $mRNA_fa -od $od/miRNA_Target/gene -cfg $detail_cfg\n";
		$cmd .="perl $miRNA_target -miRNA $miRNA_fa -target $circRNA_fa -od $od/miRNA_Target/circRNA -cfg $detail_cfg\n";
		$cmd .="perl $miRNA_target -miRNA $miRNA_fa -target $lncRNA_fa -od $od/miRNA_Target/lncRNA -cfg $detail_cfg";
		&writeSH($cmd,"$work_sh/s5.targetPrediction.sh");
		&qsub("$work_sh/s5.targetPrediction.sh");
		&qsubCheck("$work_sh/s5.targetPrediction.sh");
	}else{
		print "Please Check your previous analysis steps and results!!!\n";
		die;
	}
}

if(exists $steps{6}){
	my $human||="no";
	if($config{Project_key} eq "Human" || $config{Project_key} eq "human" || $config{medical}=~/GRCh/){
		$human = "yes";
	}
	###########DEG analysis 
	my $deg="$od/DEG_Analysis";
	`mkdir $deg`	unless(-d $deg);	
	&getDEGcfg($deg);
	if(!-d $anno){print "Please Check Your Annotation Result!!!\n";}
	$cmd ="perl $DEG_Analysis -conf $deg/gene.cfg -count $gene_count -fpkm $gene_fpkm -type FPKM -anno $anno -human $human -od $deg/gene -ppi $od/PPI/PPI.txt -gff $new_gene_gff ";	##gene DEG
	$cmd.="\nperl $DEG_Analysis -conf $deg/lncRNA.cfg -count $lnc_count -fpkm $lnc_fpkm -type FPKM -anno $anno -human $human -od $deg/lncRNA -cis $odir/LncRNA_Analysis/Lnc_target_predict/Cis_target_gene.xls ";
	$cmd.=" -trans $odir/LncRNA_Analysis/Lnc_target_predict/Trans/Trans_target_gene.xls" if(-e "$odir/LncRNA_Analysis/Lnc_target_predict/Trans/Trans_target_gene.xls");			##lncRNA DEG
	$cmd.="\nperl $DEG_Analysis -conf $deg/miRNA.cfg -count $miRNA_count -fpkm $miRNA_exp -type TPM -anno $anno -human $human -od $deg/sRNA -target $od/miRNA_Target/gene/mir2target.list"	if($switch eq "on");	##miRNA DEG
	$cmd.="\nperl $DEG_Analysis -conf $deg/circRNA.cfg -count $circ_count -fpkm $circ_exp -type $config{normalization} -anno $anno -human $human -od $deg/circRNA -target $odir/circRNA_analysis/new_name/Circ_source_gene.xls";	###circRNA DEG
	&writeSH($cmd,"$work_sh/s6.DEG.sh");
	&qsub("$work_sh/s6.DEG.sh");
	&qsubCheck("$work_sh/s6.DEG.sh");
}

############Combine analysis
if(exists $steps{7}){
	`mkdir $od/Combine`	unless(-d "$od/Combine");
	my $input_dir = "$od/Combine/Input"; mkdir $input_dir unless (-d "$input_dir");
	my $ceRNA_xls = "$od/Combine/ceRNA/ceRNA_pair_adjust_p_Sig.txt";
	&cmd_call("perl $Bin/combine/extract_needFiles.pl -idir $od -data_cfg $data_cfg -detail_cfg $detail_cfg -od $input_dir -all ");
	###Global circos
	my $loc="-lnc_gff $lnc_gff -gene_gff $gene_gff -mi_loc $odir/sRNA_Analysis/miRNA_loc/miRNA_pos.list -ideogram $input_dir/ideogram.info";
	$cmd  = "perl $Bin/combine/Global/Global.pl -lnc $lnc_fpkm -gene $gene_fpkm -circ $circ_exp -mi $miRNA_exp -relate $data_cfg -od $od/Combine/Global $loc \n";

	###Diff stat
	$cmd .= "perl $Bin/combine/Diff_stat/Diff_statis.pl $loc -cfg1 $data_cfg -cfg2 $detail_cfg -DEG $od/DEG_Analysis -od $od/Combine/Diff_stat \n";

	###co-expression
	my $sample_num = &get_data_num($data_cfg);
	if($sample_num >=5){
		$ceRNA_xls = "$od/Combine/ceRNA/coexp_ceRNA_pair_adjust_p_Sig.txt";
		$cmd .="perl $Bin/combine/co_express/co_expression.pl -circRNA $circ_exp -miRNA $miRNA_exp -lnc_mRNA $odir/LncRNA_Analysis/Lnc_target_predict/Trans/co_expression.mRNA_lncRNA.Total.xls -mRNA $gene_fpkm -lncRNA $lnc_fpkm -od $od/Combine/coexpression -cor $config{coexp_cor} -pvalue $config{coexp_p} -method $config{coexp_method} -diff $input_dir/All_DEG.xls && perl $Bin/combine/co_express/merge_coexp.pl -in $od/Combine/coexpression/ -out $od/Combine/coexpression/All_coexpression_sig.xls \n";
	}
	###ceRNA
	$cmd.="perl $Bin/combine/ceRNA/list_to_ceRNA.pl -od $od/Combine/ceRNA -num $config{ceRNA_common} -pvalue $config{ceRNA_pvalue} -FDR $config{ceRNA_fdr} -mi2c $input_dir/circRNA.mir2target.list -mi2l $input_dir/lncRNA.mir2target.list -mi2g $input_dir/gene.mir2target.list";
	$cmd.="-coexp $od/Combine/coexpression/All_coexpression_sig.xls "  if($sample_num >=5);

	###random network
        $cmd .="&& perl $Bin/combine/random/random_diff_ceRNA_main.pl -deg $input_dir -cfg $input_dir/combine.cfg -ceRNA $ceRNA_xls -anno $od/Anno_Integrate/Allgene_Anno/Result/ -od $od/Combine/random && perl $Bin/combine/cytoscape/get_ceRNA_cytoscape.pl -in $input_dir -ce $ceRNA_xls -od $od/Combine/cytoscape -cfg $input_dir/combine.cfg \n";

	###conservation
	if(exists $config{medical} || exists $config{phastCons}){
		my $dir="$od/Combine/conservation";	`mkdir $dir`	unless(-d $dir);
		&convert_conservation($lnc_gff,"$input_dir/lncRNA.bed.temp");
		my $con="perl $Bin/conservation/conservation.pl -od $dir";
		if(exists $config{phastCons} ){
			$con.=" -bw $config{phastCons} ";
		}else{	$con.=" -db $config{medical} ";}
		$cmd .="$con -gff $gene_gff -key gene && $con -bed $input_dir/lncRNA.bed.temp -key lncRNA && $con -bed $odir/circRNA_analysis/new_name/Circ.bed -key circRNA && awk '{print \$(NF-1)\"\\tcircRNA\"}' $dir/circRNA.PhastCons.score >$dir/all_type && awk '{print \$(NF-1)\"\\tlncRNA\"}' $dir/lncRNA.PhastCons.score >>$dir/all_type && awk '{print \$(NF-1)\"\\tgene\"}' $dir/gene.PhastCons.score >>$dir/all_type && grep -v '^mean0' $dir/all_type >$dir/all_type.score && sed -i 1\"imean0\\ttype\" $dir/all_type.score && rm $dir/all_type";
		$cmd .=" && $Rscript $Bin/conservation/phastcons_density_box.r $dir/all_type.score $dir All \n";
	}
	
	###lncRNA precursor
	$cmd .="perl $Bin/combine/lncRNA_precursor/Precursor_analysis.pl -fa $lncRNA_fa -key $config{miRBase} -od $od/Combine/lncRNA_precursor \n";
	
	###cytoscape analysis
	$cmd .="perl $Bin/combine/cytoscape/combine_cytoscape.pl -cfg $input_dir/combine.cfg -in $input_dir -od $od/Combine/cytoscape ";
	$cmd .="-list $input_dir/all.id_name.list " if(-e "$input_dir/all.id_name.list");
	#####################
	&writeSH($cmd,"$work_sh/s7.combine.sh");
	&qsub("$work_sh/s7.combine.sh");
	&qsubCheck("$work_sh/s7.combine.sh");
}

############Web_report
if(exists $steps{8}){
	my $qc_report_script = "/share/nas1/niepy/Tools/module_wholetrans_qc_report/make_qc_report.pl";
	###QC_Report
	my $qc_report_cmd = "perl $qc_report_script -cfg1 $data_cfg -cfg2 $detail_cfg -idir $od -odir $od/QC_Report -library lsc ";
	if(defined $exosome){$qc_report_cmd .= "-exosome ";}
	$qc_report_cmd .= "\n";
	&writeSH($qc_report_cmd,"$work_sh/s8.1.qc_report.sh");
	system("sh $work_sh/s8.1.qc_report.sh > $work_sh/s8.1.qc_report.sh.log");
	
	###produce xml/html for local
	$cmd ="perl $Bin/report/extract_result.pl -in $od -od $od/Web_Report -cfg1 $data_cfg -cfg2 $detail_cfg ";
	$cmd.="&& perl $Bin/report/extract_info.pl -in $od/Web_Report -od $od/Web_Report/HTML -cfg1 $data_cfg  -cfg2 $detail_cfg ";
	$cmd.="&& perl $Bin/report/build_xml.pl -i $od/Web_Report -cfg1 $data_cfg -cfg2 $detail_cfg -o $od/Web_Report/configtest_local.xml ";
	$cmd.="&& $python $Bin/report/bin/xml2HtmlConverter.py -i $od/Web_Report/configtest_local.xml -o $od/Web_Report ";
	###produce xml/html for biocloud
	$cmd.="&& perl $Bin/report/extract_info.pl -in $od/Web_Report -od $od/Web_Report/HTML -cfg1 $data_cfg  -cfg2 $detail_cfg -cloud -step 4 ";
	$cmd.="&& perl $Bin/report/build_xml.pl -i $od/Web_Report -cfg1 $data_cfg -cfg2 $detail_cfg -o $od/Web_Report/configtest.xml -cloud ";
	$cmd.="&& perl $Bin/report/biocloud_report.pl -in $od -data_cfg $data_cfg -detail_cfg $detail_cfg ";
	
	&writeSH($cmd,"$work_sh/s8.web_report.sh");
	&qsub("$work_sh/s8.web_report.sh");
	&qsubCheck("$work_sh/s8.web_report.sh");

	###backup
	my $cmd_biocloud = "perl $Bin/report/backup.pl -in $od";
	&writeSH($cmd_biocloud,"$work_sh/s8.2.biocloud_backup.sh");
	system("sh $work_sh/s8.2.biocloud_backup.sh > $work_sh/s8.2.biocloud_backup.sh.log");
}
################################
#	sub function
################################
sub convert_conservation{            
	my ($i,$o)=@_;
	open (IN,$i);
        open (OUT,">$o");
        my %hash;my %h;
        while (<IN>){
        	chomp;next if($_=~/^#/);
        	my($chr,$source,$type,$start,$end,$tmp1,$strand,$tmp2,$info)=split(/\t/,$_);
        	if($type=~/mRNA|^lincRNA$|^lnc_RNA$|^lncRNA$|^transcript$/){
                	$info=~/ID=(.*?)$/;
                        my $id=(split(/;/,$1))[0];
                        my $parent=(split(/;/,$1))[1];
                        $id=~s/transcript://;
                        $id=~s/mRNA://;
                        print OUT "chr$chr\t$start\t$end\t$id\n";
                }
	}
        close(IN);
        close(OUT);
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

sub get_FusionType{
	my $db=shift;
	my %genotype=(
		"GRCh37"        =>"Human.hg19",
		"GRCh38"        =>"Human.B38",
		"GRCm38"        =>"Mouse.B38",
		"Rnor_6.0"      =>"Rat.B6.0",
		"Rnor_5.0"      =>"Rat.B5.0",
	);
	return $genotype{$db};
}
sub run{
        my ($line,$shell)=@_;
        `rm -rf $shell*`        if(-e $shell);
        my @lines=split(/\,/,$line);
        open(SH,">$shell")||die $!;
        foreach my $l(@lines){
		my @multi=split(/&&/,$l);
		print SH join("&&",map{"$qsubcmd $_"}@multi),"\n";
        }
        close(SH);
	my $cmd = "$qsubcmd_once $shell";
        &run_or_die($cmd);
}

sub writeSH{
	my ($cmd,$file)=@_;
	if(-e $file){
		`rm $file`;
	}
	open(OUT,">$file")||die $!;
	print OUT "$cmd\n";
	close(OUT);
}
sub getMiRNAcfg{
	my ($data,$detail)=@_;
	`cp $data $od/Config/sRNA.cfg`;
	open(CFG,$detail)||die $!;
	open(OUT,">>$od/Config/sRNA.cfg")||die $!;
	while(<CFG>){
		chomp;
		next if($_=~/^#|^$|^\s+/);
		my @tmp=split(/\s+/,$_);
		print OUT "GENOME\t$tmp[1]\n" if($tmp[0] eq "Ref_seq");
		print OUT "$tmp[0]\t$tmp[1]\n"	if($tmp[0] eq "SPECIES_NAME"||$tmp[0] eq "SPECIES_TYPE"||$tmp[0] eq "max_map");
	}
	close(OUT);
	close(CFG);
	open(CFG,$data)||die $!;
	while(<CFG>){
		chomp;my @tmp=split(/\s+/,$_);
		$sRNA_samples{$tmp[1]}++        if($tmp[0] eq "NAME");
	}
	close(CFG);
	return "$od/Config/sRNA.cfg";
}
sub annoCfg{
	my $anno_cfg=shift;
        open(CFG,$detail_cfg)||die $!;
        open(OUT,">$anno_cfg")||die $!;
        while(<CFG>){
                chomp;
                next if($_!~/^nr|^Cog|^Kog|^Pfam|^Swissprot|^Kegg|^nt|^TrEMBL|^blast_cpu|^hmmscan_cpu|^blast_e|^blast_cut/);
                print OUT "$_\n";
        }
        print OUT "mRNA\t$novel_unigene";
        close(OUT);
        close(CFG);
}

sub rMATsCfg{
        my $rMAT_cfg = shift;
        open(CFG,$detail_cfg)||die $!;
        open(OUT,">$rMAT_cfg")||die $!;
        while(<CFG>){
                chomp;
                my @tmp = split /\s+/,$_;
                if($tmp[0] eq "Ref_ann" || $tmp[0] eq "Lib_type1" || $tmp[0] eq "ReadLength" || $tmp[0] eq "Readtype"){
                        print OUT "$tmp[0]\t$tmp[1]\n";
                }
                if($tmp[0] eq "Com"){
                        my($g1,$g2)=split /,/,$tmp[1];
                        print OUT "Com\t$sample{$g1}{gene},$sample{$g2}{gene}\n";
                }
                if($tmp[0] eq "Sep"){
                        my($g1,$g2)=split /;/,$tmp[1];
                        my @g1 = split /,/,$g1;my @g2 = split /,/,$g2;
                        print OUT "Sep\t",join(",",map{"$sample{$_}{gene}"}@g1),"\;",join(",",map{"$sample{$_}{gene}"}@g2),"\n";
                }
       }
	close(OUT);
	close(CFG);
	return($rMAT_cfg);
}

sub getDEGcfg{
	my ($od)=@_;
	my @RNAs=("lncRNA","circRNA","miRNA","gene");
	my @para=("Method_RE","Method_NR","fold","FDR","PValue");
	foreach my $rna(@RNAs){
		open(O,">$od/$rna.cfg")||die $!;
		if($rna eq "gene"){
			print O "medical\t$config{medical}\n" if(exists $config{medical});
			print O "Ref_seq\t$config{Ref_seq}\n";
			print O "Ref_ann\t$config{Ref_ann}\n";
			print O "spe_id\t$config{spe_id}\n";
			print O "TFDB\t$config{TFDB}\n";
			print O "Project_key\t$config{Project_key}\n";
			print O "score\t$config{score}\n";
		}

		foreach my $pa(@para){
			if(exists $config{"$rna.$pa"}){
				print O "$pa\t",$config{"$rna.$pa"},"\n";
			}elsif(exists $config{$pa}){
				print O "$pa\t$config{$pa}\n";
			}
		}	
		my @Com=@{$config{Com}};	
		foreach my $c(@Com){
			my ($g1,$g2)=split(/,/,$c);
			print O "Com\t$sample{$g1}{$rna},$sample{$g2}{$rna}\n";
		}
		my @Sep=@{$config{Sep}};
                foreach my $c(@Sep){
                        my ($g1,$g2)=split(/\;/,$c);
			my @g1=split(/,/,$g1);my @g2=split(/,/,$g2);
			print O "Sep\t",join(",",map{"$sample{$_}{$rna}"}@g1),"\;",join(",",map{"$sample{$_}{$rna}"}@g2),"\n";
                }
		print O "filter\t$config{filter}\n";
		close(O);
	}
}

sub relate{     ##获取不同RNA样品的对应关系
        my $file=shift;
        open(REL,$file)||die $!;
	my %sample=();
        while(<REL>){
                chomp;
                next if($_ !~/^SampleID/);
                my @tmp=split(/\s+/,$_);
                $sample{$tmp[-1]}{miRNA}=$tmp[2];
                $sample{$tmp[-1]}{lncRNA}=$tmp[3];
                $sample{$tmp[-1]}{circRNA}=$tmp[3];
                $sample{$tmp[-1]}{gene}=$tmp[3];
        }
        close(REL);
        return %sample;
}

sub converExp{
	my ($i,$o,$type)=@_;
	open(EXP,$i)||die $!;
	open(OUT,">$o")||die $!;
	my $head=<IN>;	chomp($head);	my @header=split(/\t/,$head);
	my @index=();
        foreach my $s(sort{$a cmp $b}keys %sample){
                for(my $i=0;$i<@header;$i++){
                        push @index,$i  if($header[$i] eq $sample{$s}{$type});
                }
        }

	my @samples=sort{$a cmp $b} keys %sample;
	print OUT "#ID\t",join("\t",map{ $sample{$_}{$type} }@samples),"\n";	
	while(<IN>){
		
	}
	close(OUT);
	close(EXP);
}

sub readcfg{
        my $cfg=shift;
        open(CFG,$cfg)||die $!;
        while(<CFG>){
                next if($_=~/^#/);
                chomp;my @tmp=split(/\s+/,$_);
		if($tmp[0] eq "Com"||$tmp[0] eq "Sep"){
			push @{$config{$tmp[0]}},$tmp[1];
		}else{
			$config{$tmp[0]}=$tmp[1];
		}
        }
        close(CFG);
}

sub getmiRNApara{
	my($mi,$cfg)=@_;
	open(MI,$mi)||die $!;
	while(<MI>){
		chomp;next if($_=~/^#|^$/);
		`sed -i "1i$_" $cfg`	if($_=~/^RNAhybrid|^miRanda|^TargetFinder/);
		my ($key,$con)=split(/\s+/,$_,2);
		`sed -i "1itype\t$con" $cfg`	if($key eq "SPECIES_TYPE");
	}
	close(MI);
}

sub unique{
	my $s=shift;
	my @set=@{$s};
	my %hash=();
	foreach my $s(@set){
		$hash{$s}++;
	}
	return keys %hash;
}



sub extractKey{
	my ($node,$type,$ratio)=@_;
	$ratio||=$config{keyratio};
	##node score: gene	score
	##node type:  rna1	mRNA
	open(TYPE,$type)||die $!;
	my %node_types=();
	while(<TYPE>){
		chomp;
		my @tmp=split(/\s+/,$_);
		$node_types{$tmp[0]}=$tmp[1];
	}
	close(TYPE);
	print "$node\n";
	my $dir=dirname $node;
	my $num=&readLine($node);
	my $keyNum=int($num*$ratio);
	$keyNum=$config{keymin} if($keyNum < $config{keymin});
	print "$keyNum\n";
	my $i=0;
	my $j=0;
	open(SCORE,$node)||die $!;
	open(KEY,">$dir/key_RNA.txt")||die $!;
	open(KEYM,">$dir/key_mRNA.txt")||die $!;
	print KEYM "#Key\n";
	print KEY  "#Key\n";
	while(<SCORE>){
		chomp;
		my ($rna,$score)=split(/\t/,$_);
		print KEY "$rna\n"; 
		if($node_types{$rna} eq "mRNA"){
			print KEYM "$rna\n";
			$j++;	
		}
		$i++;
		last if($i>=$keyNum);
	}
	close(SCORE);
	close(KEY);
	close(KEYM);
	return $j;
}
sub readLine{
	my $file=shift;
	open(FILE,$file)||die $!;
	my @tmp=<FILE>;
	close(FILE);
	return scalar(@tmp);
}
#############################################################################
#############################################################################
sub cmd_call {
        print "@_\n";
        system(@_) == 0 or print "Error: system @_ failed: $?\n";
}

sub readConfig{
	my $configFile=shift;
	my $d=Config::General->new(-ConfigFile => "$configFile");
	my %config=$d->getall;	
	return %config;
}

sub qsub(){
        my $shfile= shift;
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $shfile --queue $config{Queue_type}";
        &run_or_die($cmd);              
}

sub qsubCheck{
        my $sh = shift;
        my @Check_file = glob "$sh*.qsub/*.Check";
        my @sh_file    = glob "$sh*.qsub/*.sh";
        if ( $#sh_file != $#Check_file ) {
                print "Their Some Error Happend in $sh qsub, Please Check..\n";
                die;
        }else {
                print "$sh qsub is Done!\n";
        }
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
}
sub show_log(){
        my ($txt) = @_ ;
        my $time = time();
        my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime($time);
        $wday = $yday = $isdst = 0;
        my $Time=sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
        print "$Time:\t$txt\n" ;
}
sub USAGE{
	my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
	Options:
	-cfg1	<file>	input data config file, forced
	-cfg2	<file>	input detail config file, forced
	-o	<path>	output path, forced
	-step	<str>	step number,separated by comma (default all)
			1. lncRNA data_assess
			2. 1) sRNA data_assess, lncRNA map
			   2) circRNA, miRNA mapping
			3. lncRNA, circRNA, miRNA predict
			4. Structure Analysis(AS,SNP,Fusion,miRNA edit/Family,new Gene Function)
			5. miRNA target prediction
			6. DEG analysis and anno
			7. Combine analysis (ceRNA, co-expression, global/diff circos, conservation, key RNA)
			8. extract result and produce web_report
				
	-switch	<str>	default on, if off will not process miRNA data
	-exosome	whether is exosome or not

	-h	Help

Example: perl Combine_whole_transcriptome.pl -cfg1 config/data.cfg -cfg2 config/detail.cfg -o Analysis -step 1,2,3,4,5,6,7,8

USAGE
	print $usage;
	exit;
}


