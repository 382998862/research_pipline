use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path getcwd);
use threads;
use newPerlBase;
my ($data_cfg, $detail_cfg,$od,$step,$onluy_lnc,$only_circ,$only_sRNA,$exosome);
GetOptions(
        "h|?"           =>\&USAGE,
	"cfg1:s"	=>\$data_cfg,
	"cfg2:s"	=>\$detail_cfg,
	"o:s"     	=>\$od,
	"step:s"	=>\$step,
	"lnc:s"		=>\$only_lnc,
	"circ:s"	=>\$only_circ,
	"sRNA:s"	=>\$only_sRNA,
	"exosome"	=>\$exosome,
)or &USAGE;
&USAGE unless ($data_cfg and $detail_cfg and $od and (defined $only_lnc || defined $only_circ || defined $only_sRNA));

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
my $rMATs               ="$Bin/differential_AS_analysis/diff_AS_analysis.pl";
my $miRNA_family	=$sys{miRNA_family};
my $miRNA_edit		="$Bin/sRNA/edit/miRNA_Edit.pl";
my $HMDD		="/share/nas2/genome/bmksoft/pipeline/sRNA_process/v2.5/bin/HMDD_Annotation/HMDD_annotation.pl";

my $circRNA_bwa		="$Bin/circRNA/bwa/bwa_process.pl";
my $circRNA_map		=$sys{circRNA_map};
my $circRNA_pre		="$Bin/circRNA/circRNA_identify/circRNA_analysis.pl";
my $count2len 		="$Bin/circRNA/circRNA_identify/bin/count_to_length.pl";

my $DEG_Analysis	="$Bin/DEG_Analysis/DEG_Analysis.pl";
my $gene_fusion		="/share/nas1/wenyh/develop/lncRNA/medical_lncRNA_v1.3/gene_fusion/fusionmap3.pl";
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
if(defined $only_lnc){
	`mkdir $work_sh/lncRNA`		unless(-d "$work_sh/lncRNA");
}
if(defined $only_circ){
	`mkdir $work_sh/circRNA`	unless(-d "$work_sh/circRNA");
}
if(defined $only_sRNA){
	`mkdir $work_sh/miRNA`		unless(-d "$work_sh/miRNA");
}
my ($odir,$lines,$cmd,$miRNA_line);

my($lnc_fpkm,$lnc_count,$circ_exp,$circ_count,$gene_fpkm,$gene_count,$miRNA_exp,$miRNA_count);
$gene_fpkm      ="$od/Basic_Analysis/Hisat_Stringtie/prepDE/All_gene_fpkm.list";
$gene_count     ="$od/Basic_Analysis/Hisat_Stringtie/prepDE/All_gene_counts.list";
$circ_exp       ="$od/Basic_Analysis/circRNA_analysis/expression/All_gene_expression.list";
$circ_count     ="$od/Basic_Analysis/circRNA_analysis/expression/All_gene_counts.list";
$lnc_fpkm       ="$od/Basic_Analysis/Hisat_Stringtie/prepDE/All_lncRNA_fpkm.list";
$lnc_count      ="$od/Basic_Analysis/Hisat_Stringtie/prepDE/All_lncRNA_counts.list";
$miRNA_exp      ="$od/Basic_Analysis/sRNA_Analysis/miRDeep2/miRNA_Quantify/All_miRNA_expression.list";
$miRNA_count    ="$od/Basic_Analysis/sRNA_Analysis/miRDeep2/miRNA_Quantify/All_miRNA.count.list";
$odir="$od/Basic_Analysis";

my($miRNA_fa,$circRNA_fa,$lncRNA_fa,$mRNA_fa,$novel_unigene);
$circRNA_fa  ="$od/Basic_Analysis/circRNA_analysis/circRNA.fa";
$lncRNA_fa   ="$od/Basic_Analysis/LncRNA_Analysis/lnc_predict/All_filter_final.fa";

my($new_gene_gff,$gene_gff,$lnc_gff,$anno);
$new_gene_gff	="$od/Basic_Analysis/Hisat_Stringtie/genePredict/final_track/$index.newGene_final.filtered.gff";
$lnc_gff        ="$od/Basic_Analysis/LncRNA_Analysis/lnc_predict/All_filter_final.gff";
if(defined $only_sRNA || (defined $only_circ && !defined $only_lnc)){
	$gene_gff = $config{Ref_ann};
	$anno = "$config{Known_anno}/Result";
	$mRNA_fa = $config{Known_unigene};
}else{
	$gene_gff	="$od/Basic_Analysis/Hisat_Stringtie/genePredict/All_gene_final.gff";
	$anno        ="$od/Anno_Integrate/Allgene_Anno/Result";
	$mRNA_fa     ="$od/Basic_Analysis/Hisat_Stringtie/genePredict/All.longest_transcript.fa";
}

my $configdir = dirname $data_cfg;
if(defined $only_sRNA ){
	$miRNA_fa ="$od/Basic_Analysis/sRNA_Analysis/miRDeep2/miRNA_Quantify/All_miRNA.expressed.fa";
}else{
	if(exists $config{miRNA_fa}){
		$miRNA_fa = $config{miRNA_fa};
	}else{
		system("perl $Bin/sRNA/basic_analysis/bin/select_fa.pl -i $config{miRBase} -fa /share/nas2/database/miRBase/v21/mature.fa -o $configdir/$config{miRBase}.mature.fa");
		print "perl $Bin/sRNA/basic_analysis/bin/select_fa.pl -i $config{miRBase} -fa /share/nas2/database/miRBase/v21/mature.fa -o $configdir/$config{miRBase}.mature.fa\n";
		$miRNA_fa = "$configdir/$config{miRBase}.mature.fa";
	}
}

############data_assess

if(exists $steps{1}){
	$odir="$od/Data_Assess";
	&writeSH("perl $lncRNA_data_assess -config $data_cfg -detailfig $detail_cfg -outdir $odir ","$work_sh/lncRNA/s1.data_assess.sh") if (defined $only_lnc);
	&writeSH("perl $lncRNA_data_assess -config $data_cfg -detailfig $detail_cfg -outdir $odir ","$work_sh/circRNA/s1.data_assess.sh") if (defined $only_circ && !defined $only_lnc);
	&writeSH("perl $miRNA_data_assess -Q $config{Q} -main_cfg $data_cfg -word_cfg $detail_cfg -od $odir -min $config{min_len} -max $config{max_len}","$work_sh/miRNA/s1.data_assess.sh") if(defined $only_sRNA);
	$lines ="$work_sh/lncRNA/s1.data_assess.sh" if (defined $only_lnc);
	$lines ="$work_sh/circRNA/s1.data_assess.sh" if(!defined $only_lnc && defined $only_circ);
	$lines ="$work_sh/miRNA/s1.data_assess.sh" if(defined $only_sRNA);
	&run($lines,"$work_sh/s1.data_assess.sh");
	&qsubCheck("$work_sh/s1.data_assess.sh");
}

############Basic_Analysis
if(exists $steps{2}){
	$odir="$od/Basic_Analysis";
	`mkdir $odir` unless (-d $odir);
	if(defined $only_lnc){
		###lncRNA map
		&writeSH("perl $Bin/lncRNA/hisat_string/hisat2_stringtie.pl -cfg1 $data_cfg -cfg2 $detail_cfg -od $odir/Hisat_Stringtie && perl $Bin/lncRNA/hisat_string/map_stat/Statistics_draw_hisat.pl -cfg2 $detail_cfg -in $odir/Hisat_Stringtie -od $odir/Hisat_Stringtie/Mapped ","$work_sh/lncRNA/s2.lncRNA_map.sh");
	}
	
	if(defined $only_circ){
		###circRNA map
		`mkdir $odir/circRNA_Bwa` unless(-d "$odir/circRNA_Bwa");
		`echo "Ref_seq\t$config{Ref_seq}\nRef_ann\t$gene_gff" >$odir/circRNA_Bwa/detail.cfg`;
		`echo "medical\t$config{medical}" >>$odir/circRNA_Bwa/detail.cfg`	if(exists $config{medical});
		&writeSH("perl $circRNA_bwa -cfg1 $data_cfg -cfg2 $odir/circRNA_Bwa/detail.cfg -od $odir/circRNA_Bwa && perl $circRNA_map -id $odir/circRNA_Bwa/Map_Stat/ -o $odir/circRNA_Bwa/Map_Stat/All.mappedStat.xls && perl $count2len -cfg $data_cfg -o $odir/circRNA_Bwa/Map_Stat/All.readLength.xls","$work_sh/circRNA/s2.cicRNA_map.sh");
	}
	#########
	if(defined $only_lnc && defined $only_circ){
		$lines ="$work_sh/lncRNA/s2.lncRNA_map.sh,$work_sh/circRNA/s2.cicRNA_map.sh";
		&run($lines,"$work_sh/s2.lncRNA_library_map.sh");
	}elsif(defined $only_lnc && !defined $only_circ){
		$lines ="$work_sh/lncRNA/s2.lncRNA_map.sh";
		&run($lines,"$work_sh/s2.lncRNA_library_map.sh");
	}elsif(!defined $only_lnc && defined $only_circ){
		$lines ="$work_sh/circRNA/s2.cicRNA_map.sh";
		&run($lines,"$work_sh/s2.lncRNA_library_map.sh");
	}
	#########
	if(defined $only_sRNA){
		###sRNA map
		`mkdir $odir/sRNA_Alignment`	unless(-d "$odir/sRNA_Alignment");
		&writeSH("perl $miRNA_map -od $odir/sRNA_Alignment -id $od/Data_Assess -main_cfg $detail_cfg -data_cfg $data_cfg -filter 1 -genome $config{Ref_seq} -gff $gene_gff ","$work_sh/miRNA/s2.alignment.sh");
		$lines ="$work_sh/miRNA/s2.alignment.sh";
		&run($lines,"$work_sh/s2.sRNA_library_map.sh");
	}
	&qsubCheck("$work_sh/s2.lncRNA_library_map.sh") if(-e "$work_sh/s2.lncRNA_library_map.sh");
	&qsubCheck("$work_sh/s2.sRNA_library_map.sh") if(-e "$work_sh/s2.sRNA_library_map.sh");
}

if(exists $steps{3}){
	if(defined $only_lnc){
		###lncRNA identify
		&writeSH("perl $Bin/lncRNA/lncRNA_Analysis/Lncrna_analysis.v3.pl -cfg $detail_cfg -gtf $odir/Hisat_Stringtie/LncPredict/merged_filter.gtf -gff $gene_gff -od $odir/LncRNA_Analysis -prep $odir/Hisat_Stringtie/prepDE -type $config{type} -db $config{db}","$work_sh/lncRNA/s3.identify_target.sh");
		###Function analysis
		$novel_unigene="$odir/Hisat_Stringtie/genePredict/final_track/$index.newGene.longest_transcript.fa";
		my $anno_cfg="$od/Config/anno.cfg"; &annoCfg($anno_cfg);
		my $cmd = "perl $gene_anno  --cfg $anno_cfg --od $od/Anno_Integrate/New_Anno --queue $config{Queue_type} --eggNOG --pfam --kog ";
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
		#my $rMATs_cfg = "$od/Config/rMATs.cfg";&rMATsCfg($rMATs_cfg);
		&writeSH("perl $rMATs -cfg $detail_cfg -bamdir $odir/Hisat_Stringtie/Hisat -od $od/differential_AS_analysis","$work_sh/lncRNA/s9.different_AS.sh");
		
		$lines ="$work_sh/lncRNA/s3.identify_target.sh,$work_sh/lncRNA/s4_Anno_Integrate.sh,$work_sh/lncRNA/s5.SNP.sh,$work_sh/lncRNA/s6.AS.sh,$work_sh/lncRNA/s7.Gene_Structure_Optimize.sh,$work_sh/lncRNA/s8.newGene_PPI.sh,$work_sh/lncRNA/s9.different_AS.sh";
		
		###Gene fusion
		if(exists $config{medical}){
			my $fusion_db=&get_FusionType($config{medical});
			&writeSH("perl $gene_fusion --indir $odir/Hisat_Stringtie/Hisat --type $fusion_db --od $od/Gene_Fusion --script_cfg  /share/nas1/wenyh/develop/lncRNA/medical_lncRNA_v1.3/gene_fusion/script.cfg --cfg2 $detail_cfg","$work_sh/lncRNA/s10.fusion.sh");
			$lines.=",$work_sh/lncRNA/s10.fusion.sh";
        	}

	}
	if(defined $only_circ){
		##circRNA identify
		&writeSH("perl $circRNA_pre -data_cfg $data_cfg -detail_cfg $detail_cfg -od $odir/circRNA_analysis -gff $gene_gff ","$work_sh/circRNA/s3.identify.sh");
	}

	#########
	if(defined $only_lnc && defined $only_circ){
		$lines .= ",$work_sh/circRNA/s3.identify.sh";
		&run($lines,"$work_sh/s3.identify_structure.sh");
	}elsif(defined $only_lnc && defined !$only_circ){
		&run($lines,"$work_sh/s3.identify_structure.sh");
	}elsif(!defined $only_lnc && defined $only_circ){
		$lines = "$work_sh/circRNA/s3.identify.sh";
		&run($lines,"$work_sh/s3.identify_structure.sh");
	}
	#########
	if(defined $only_sRNA){
		##miRNA identify
		$config{miRBase}=$config{Species_type}  if(!exists $config{miRBase});
		$cmd="perl $miRNA_identify -mature /share/nas2/database/miRBase/v21/mature.fa -miR $config{miRBase} -idir $odir/sRNA_Alignment -od $odir/sRNA_Analysis -maincfg $miRNA_cfg -min $config{min_len} -max $config{max_len} ";
		if(exists $config{medical}){
			$cmd.="&& perl $HMDD -in $odir/sRNA_Analysis/miRDeep2/miRNA_Quantify/All_miRNA_expression.list -od $odir/sRNA_Analysis/HMDD"	if($config{medical} eq "GRCh37" || $config{medical} eq "GRCh38");
		}
		&writeSH($cmd,"$work_sh/miRNA/s3.identify.sh");
		##miRNA edit and family
		my $sam=join(",",keys %sRNA_samples);
		$cmd ="perl $miRNA_family -i $miRNA_fa -type $config{Species_type} -od $od/miRNA_Family ";
		$cmd.="&& perl $miRNA_edit -idir $odir -od $od/miRNA_edit -sample $sam ";
		&writeSH("perl $miRNA_family -i $miRNA_fa  -type $config{Species_type} -od $od/miRNA_Family && perl $miRNA_edit -idir $odir -od $od/miRNA_edit -sample $sam ","$work_sh/miRNA/s4.miRNA_family_edit.sh");
		$lines = "$work_sh/miRNA/s3.identify.sh && $work_sh/miRNA/s4.miRNA_family_edit.sh";
		&run($lines,"$work_sh/s3.identify_structure.sh");
	}
	&qsubCheck("$work_sh/s3.identify_structure.sh");
}

if(exists $steps{4}){
	############miRNA target prediction
	my $cmd=();
	if(defined $only_sRNA){
		$cmd ="perl $miRNA_target -miRNA $miRNA_fa -target $mRNA_fa -od $od/Target_Predict/ -cfg $detail_cfg\n";
	}else{
		$cmd ="perl $miRNA_target -miRNA $miRNA_fa -target $mRNA_fa -od $od/Target_Predict/miRNA-mRNA -cfg $detail_cfg\n";
	}
	my $cmd1 ="perl $miRNA_target -miRNA $miRNA_fa -target $circRNA_fa -od $od/Target_Predict/miRNA-circRNA -cfg $detail_cfg\n";
	my $cmd2 ="perl $miRNA_target -miRNA $miRNA_fa -target $lncRNA_fa -od $od/Target_Predict/miRNA-lncRNA -cfg $detail_cfg";
	
	if(defined $only_lnc && defined $only_circ){
		$cmd .= "$cmd1"."$cmd2";
	}elsif(defined $only_lnc && !defined $only_circ){
		$cmd .= $cmd2;
	}elsif(!defined $only_lnc && defined $only_circ){
		$cmd .= $cmd1;
	}

	&writeSH($cmd,"$work_sh/s4.targetPrediction.sh");
	&qsub("$work_sh/s4.targetPrediction.sh");
	&qsubCheck("$work_sh/s4.targetPrediction.sh");
}

if(exists $steps{5}){
	###########DEG analysis 
	my $deg="$od/DEG_Analysis";
	`mkdir $deg`	unless(-d $deg);
	my @TYPE;
	if(defined $only_sRNA){push @TYPE,"miRNA";}
	if(defined $only_lnc){push @TYPE,"lncRNA";push @TYPE,"gene";}
	if(defined $only_circ){push @TYPE,"circRNA";}
	my $TYpe=join(",",@TYPE);
	&getDEGcfg($deg,$TYpe);
	my $cmd =();
	my $human||="no";
        if($config{Project_key} eq "Human" || $config{Project_key} eq "human" || $config{medical}=~/GRCh/){
                $human = "yes";
        }
	if(defined $only_sRNA){
	###miRNA DEG
		$cmd = "perl $DEG_Analysis -conf $deg/miRNA.cfg -count $miRNA_count -fpkm $miRNA_exp -type TPM -anno $anno -human $human -od $deg -target $od/Target_Predict/mir2target.list";
	}
	###gene DEG
	my $cmd1 ="perl $DEG_Analysis -conf $deg/gene.cfg -count $gene_count -fpkm $gene_fpkm -type FPKM -anno $anno -human $human -gff $new_gene_gff -od $deg/gene -ppi $od/PPI/PPI.txt \n";
	###lncRNA DEG
	my $cmd2 ="perl $DEG_Analysis -conf $deg/lncRNA.cfg -count $lnc_count -fpkm $lnc_fpkm -type FPKM -anno $anno -human $human -od $deg/lncRNA -cis $odir/LncRNA_Analysis/Lnc_target_predict/Cis_target_gene.xls ";
	$cmd2.="-trans $odir/LncRNA_Analysis/Lnc_target_predict/Trans/Trans_target_gene.xls " if(-e "$odir/LncRNA_Analysis/Lnc_target_predict/Trans/Trans_target_gene.xls");
	###circRNA DEG
	my $cmd3;
	if(defined $only_lnc && defined $only_circ){
		$cmd3="perl $DEG_Analysis -conf $deg/circRNA.cfg -count $circ_count -fpkm $circ_exp -type $config{normalization} -anno $anno -human $human -od $deg/circRNA -target $odir/circRNA_analysis/new_name/Circ_source_gene.xls ";
		$cmd = "$cmd1"."$cmd2"."\n"."$cmd3";
	}elsif(defined $only_lnc && !defined $only_circ){
		$cmd = "$cmd1"."$cmd2";
	}elsif(defined !$only_lnc && defined $only_circ){
		$cmd3 = "perl $DEG_Analysis -conf $deg/circRNA.cfg -count $circ_count -fpkm $circ_exp -type $config{normalization} -anno $anno -human $human -od $deg -target $odir/circRNA_analysis/new_name/Circ_source_gene.xls ";
		$cmd = $cmd3;
	}
	&writeSH($cmd,"$work_sh/s5.DEG.sh");
	&qsub("$work_sh/s5.DEG.sh");
	&qsubCheck("$work_sh/s5.DEG.sh");
}

############Combine analysis
if(exists $steps{6} && (defined $only_lnc || defined $only_circ)){
	`mkdir $od/Personalization`	unless(-d "$od/Personalization");
	my $cmd=();
	if(defined $only_lnc){
		`mkdir $od/Personalization/lncRNA_precursor` unless(-d "$od/Personalization/lncRNA_precursor");
		$cmd ="perl $Bin/combine/lncRNA_precursor/Precursor_analysis.pl -fa $lncRNA_fa -key $config{miRBase} -od $od/Personalization/lncRNA_precursor \n";
	}
	
	if(exists $config{medical} || exists $config{phastCons}){
                my $dir="$od/Personalization/conservation";
		`mkdir $dir` unless(-d $dir);
                my $con="perl $Bin/conservation/conservation.pl -od $dir";
		
                open (IN,"$lnc_gff");
		print "$lnc_gff\n";
                open (OUT,">$dir/lncRNA.bed.temp");
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

                if(exists $config{phastCons} ){
                        $con.=" -bw $config{phastCons} ";
                }else{
			$con.=" -db $config{medical} ";
		}
		my $cmd1 = "$con -gff $gene_gff -key gene && $con -bed $dir/lncRNA.bed.temp -key lncRNA";
		my $cmd2 = "$con -bed $odir/circRNA_analysis/new_name/Circ.bed -key circRNA";
		if(defined $only_lnc && defined $only_circ){
			$cmd .= "$cmd1 && $cmd2 && awk '{print \$(NF-1)\"\\tcircRNA\"}' $dir/circRNA.PhastCons.score >$dir/all_type && awk '{print \$(NF-1)\"\\tlncRNA\"}' $dir/lncRNA.PhastCons.score >>$dir/all_type && awk '{print \$(NF-1)\"\\tgene\"}' $dir/gene.PhastCons.score >>$dir/all_type && grep -v '^mean0' $dir/all_type >$dir/all_type.score && sed -i 1\"imean0\\ttype\" $dir/all_type.score && rm $dir/all_type && $Rscript $Bin/conservation/phastcons_density_box.r $dir/all_type.score $dir All ";
		}elsif(defined $only_lnc && !defined $only_circ){
			$cmd .= "$cmd1 && awk '{print \$(NF-1)\"\\tlncRNA\"}' $dir/lncRNA.PhastCons.score > $dir/all_type && awk '{print \$(NF-1)\"\\tgene\"}' $dir/gene.PhastCons.score >>$dir/all_type && grep -v '^mean0' $dir/all_type >$dir/all_type.score && sed -i 1\"imean0\\ttype\" $dir/all_type.score && rm $dir/all_type && $Rscript $Bin/conservation/phastcons_density_box.r $dir/all_type.score $dir All ";
		}elsif(!defined $only_lnc && defined $only_circ){
			$cmd .= "$cmd2 && awk '{print \$(NF-1)\"\\tcircRNA\"}' $dir/circRNA.PhastCons.score >$dir/all_type && grep -v '^mean0' $dir/all_type >$dir/all_type.score && sed -i 1\"imean0\\ttype\" $dir/all_type.score && rm $dir/all_type && $Rscript $Bin/conservation/phastcons_density_box.r $dir/all_type.score $dir All "
		}
	}
		&writeSH($cmd,"$work_sh/s6.personalize.sh");
        	&qsub("$work_sh/s6.personalize.sh");
		&qsubCheck("$work_sh/s6.personalize.sh");
}
############Web_report
if(exists $steps{7}){
	###QC_Report
	my $qc_report_script = "/share/nas1/niepy/Tools/module_wholetrans_qc_report/make_qc_report.pl";
	my $qc_report_cmd = "perl $qc_report_script -cfg1 $data_cfg -cfg2 $detail_cfg -idir $od -odir $od/QC_Report ";
	if(defined $only_sRNA){	$qc_report_cmd .= "-library s ";}
	if(defined $only_lnc && defined $only_circ){
		$qc_report_cmd .= "-library lc ";
	}elsif(defined $only_lnc && !defined $only_circ){
		$qc_report_cmd .= "-library l ";
	}elsif(defined $only_circ && !defined $only_lnc){
		$qc_report_cmd .= "-library c ";
	}
	if(defined $exosome){ $qc_report_cmd = "-exosome ";}
	$qc_report_cmd .= "\n";
	&writeSH($qc_report_cmd,"$work_sh/s7.1.qc_report.sh");
	system("sh $work_sh/s7.1.qc_report.sh > $work_sh/s7.1.qc_report.sh.log");
	###produce xml/html in Web_Report
	my $cmd ="perl $Bin/report/extract_result.1.pl -in $od -od $od/Web_Report -cfg1 $data_cfg -cfg2 $detail_cfg ";
	$cmd.="&& perl $Bin/report/extract_info.1.pl -in $od/Web_Report -od $od/Web_Report/HTML -cfg1 $data_cfg  -cfg2 $detail_cfg ";
	if(defined $only_sRNA){
		$cmd.="&& perl $Bin/report/build_xml.sRNA.pl -i $od/Web_Report -cfg1 $data_cfg -cfg2 $detail_cfg -o $od/Web_Report/configtest_local.xml ";
	}
	if(defined $only_circ || defined $only_lnc){
		$cmd.="&& perl $Bin/report/build_xml.lncRNA.pl -i $od/Web_Report -cfg1 $data_cfg -cfg2 $detail_cfg -o $od/Web_Report/configtest_local.xml ";
	}
	$cmd.="&& $python $Bin/report/bin/xml2HtmlConverter.py -i $od/Web_Report/configtest_local.xml -o $od/Web_Report ";
	$cmd.="&& perl $Bin/report/extract_info.1.pl -in $od/Web_Report -od $od/Web_Report/HTML -cfg1 $data_cfg  -cfg2 $detail_cfg -cloud -step 4 ";
	if(defined $only_sRNA){
		$cmd.="&& perl $Bin/report/build_xml.sRNA.pl -i $od/Web_Report -cfg1 $data_cfg -cfg2 $detail_cfg -o $od/Web_Report/configtest.xml -cloud";
	}
	if(defined $only_circ || defined $only_lnc){
		$cmd.="&& perl $Bin/report/build_xml.lncRNA.pl -i $od/Web_Report -cfg1 $data_cfg -cfg2 $detail_cfg -o $od/Web_Report/configtest.xml -cloud ";
	}
	$cmd.= "&& perl $Bin/report/biocloud_report.pl -in $od -data_cfg $data_cfg -detail_cfg $detail_cfg ";

	&writeSH($cmd,"$work_sh/s7.web_report.sh");
	&qsub("$work_sh/s7.web_report.sh");

	my $cmd_biocloud; 
	if(defined $only_lnc || defined $only_circ){
		$cmd_biocloud .= "perl $Bin/report/backup.lncRNA.pl -in $od ";
	}elsif(defined $only_sRNA){
		$cmd_biocloud .= "perl $Bin/report/backup.sRNA.pl -in $od ";
	}
	
        &writeSH($cmd_biocloud,"$work_sh/s7.2.biocloud_backup.sh");
	system("sh $work_sh/s7.2.biocloud_backup.sh > $work_sh/s7.2.biocloud_backup.sh.log");
}
################################
#	sub function
################################

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
	my $qsubcmd="sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent --queue $config{Queue_type} ";
	my $qsubcmd_once="sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --independent --queue $config{Queue_type} ";
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
        my $rMATs_cfg = shift;
        open(CFG,$detail_cfg)||die $!;
        open(OUT,">$rMATs_cfg")||die $!;
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
	return($rMATs_cfg);
}

sub getDEGcfg{
	my $od=shift;
	my $Type=shift;
#	my @RNAs=("lncRNA","circRNA","miRNA","gene");
	my @RNAs=split /,/,$Type;
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
			print O "score\t$config{score}";
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
#			print O "Com\t$sample{$g1}{$rna},$sample{$g2}{$rna}\n";
			print O "Com\t$g1,$g2\n";
		}
		my @Sep=@{$config{Sep}};
                foreach my $c(@Sep){
                        my ($g1,$g2)=split(/\;/,$c);
			my @g1=split(/,/,$g1);my @g2=split(/,/,$g2);
#			print O "Sep\t",join(",",map{"$sample{$_}{$rna}"}@g1),"\;",join(",",map{"$sample{$_}{$rna}"}@g2),"\n";
			print O "Sep\t$g1;$g2\n";
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
			1. lncRNA/circ/sRNA data_assess
			2. lncRNA/circ/sRNA map
			3. lncRNA/circRNA/miRNA predict and Structure Analysis(AS,SNP,Fusion,miRNA edit/Family,new Gene Function)
			4. target prediction
			5. DEG analysis and anno
			6. lncRNA/circRNA conservation;lncRNA precursor
			7. extract result and produce web_report	
	-lnc	defined, only analysis lncRNA and creat lncRNA report,can exist with "-circ" and creat lncRNA and circRNA report separately
	-circ	defined, only analysis circRNA and creat circRNA report,can exist with "-lnc" and creat lncRNA and circRNA report separately
	-sRNA	defined, only analysis sRNA and creat sRNA report
	-exosmoe	whether is exosome or not
	
	-h	Help
Example:
	only lncRNA: perl $0 -cfg1 config/data.cfg -cfg2 config/detail.cfg -o Analysis -step 1,2,3,4,5,6,7 -lnc
	only circRNA: perl $0 -cfg1 config/data.cfg -cfg2 config/detail.cfg -o Analysis -step 1,2,3,4,5,6,7 -circ
	only sRNA: perl $0 -cfg1 config/data.cfg -cfg2 config/detail.cfg -o Analysis -step 1,2,3,4,5,7 -sRNA
	lncRNA and circRNA: perl $0 -cfg1 config/data.cfg -cfg2 config/detail.cfg -o Analysis -step 1,2,3,4,5,6,7 -lnc -circ

USAGE
	print $usage;
	exit;
}
