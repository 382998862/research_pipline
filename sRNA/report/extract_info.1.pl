use Getopt::Long;
use Getopt::Std;
use Config::General;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);

my ($in,$od,$data_cfg,$detail_cfg,$cloud,$step);
GetOptions(
	"h|?"   	        =>\&USAGE,
	"od:s"     		=>\$od,
	"in:s"			=>\$in,
	"cfg1:s"		=>\$data_cfg,
	"cfg2:s"		=>\$detail_cfg,
	"cloud"			=>\$cloud,
	"step:s"		=>\$step,
)or &USAGE;
&USAGE unless ($in and $data_cfg and $detail_cfg);
$od||="$in/HTML";
`mkdir $od` unless(-d $od);
$in =abs_path($in);
$od =abs_path($od);
$data_cfg =abs_path($data_cfg);
$detail_cfg =abs_path($detail_cfg);

`mkdir $od/template` unless (-d "$od/template");
my $temlate_dir = "$od/template";

my %config=&readConfig($detail_cfg);
my $deg_flag =1;
if(!exists $config{Sep} && !exists $config{Com}){
	$deg_flag=0;
}

my %step;
$step||="1,2,3";
my @step = split /,/,$step;
foreach my $s(@step){
	$step{$s}=1;
}

if(!exists $step{4} && defined $cloud){
        print "Step 4 is for generate xml for cloud ,if you want to do this ,please set you options if '-step 4 -cloud'\n";
        die;
}

my $lnc_dir = (glob "$in/BMK_*_LncRNA")[0];
my $mRNA_dir = (glob "$in/BMK_*_mRNA")[0];
my $circ_dir = (glob "$in/BMK_*_circRNA")[0];
my $sRNA_dir = (glob "$in/BMK_*_miRNA")[0];
my $struc_dir = (glob "$in/BMK_*_Structure")[0];

if(defined $lnc_dir && defined $circ_dir){
	`cp $Bin/readme/lncRNA_circRNA_materials_methods.pdf $in/BMK_1_rawData/materials_and_methods.pdf`;
}elsif(defined $lnc_dir && !defined $circ_dir){
	`cp $Bin/readme/lncRNA_materials_methods.pdf $in/BMK_1_rawData/materials_and_methods.pdf`;
}elsif(!defined $lnc_dir && defined $circ_dir){
	`cp $Bin/readme/circRNA_materials_methods.pdf $in/BMK_1_rawData/materials_and_methods.pdf`;
}

if(defined $sRNA_dir){
	`cp $Bin/readme/sRNA_materials_methods.pdf $in/BMK_1_rawData/materials_and_methods.pdf`;
}

my $disease_dir = (glob "$in/BMK_*LncRNA/BMK_*_Lnc_disease")[0];
my $precursor_dir = (glob "$in/BMK_*_LncRNA/BMK_*_Lnc_precursor")[0];

my $family_dir = (glob "$in/BMK_*_miRNA/BMK_*_miRNA_Family")[0];
my $hmdd_dir = (glob "$in/BMK_*_miRNA/BMK_*_HMDD")[0];
my $edit_dir = (glob "$in/BMK_*_miRNA/BMK_*_miRNA_Edit")[0];



if(exists $step{1}){
	print "\n------extract informations Start-----\n";
`cp -r $Bin/images/* $temlate_dir`;
`cp -r $Bin/demo $temlate_dir`;
`cp -r $Bin/demo/*png $temlate_dir`;

my %sample=&relate($data_cfg);

if(-e "$in/BMK_1_rawData/BMK_1_Data_Assess/nonRib.AllSample_GC_Q.stat.xls"){
	open(IN,"$in/BMK_1_rawData/BMK_1_Data_Assess/nonRib.AllSample_GC_Q.stat.xls")||die $!;
	open(OUT,">$in/BMK_1_rawData/BMK_1_Data_Assess/nonRib.AllSample_GC_Q.stat.xls.temp")||die $!;
	my %lncRNA_data=();
	my $lncRNA_total;
	my $lncRNA_Q30=100;
	while(<IN>){
		chomp;
		if($_=~/^#B/){
			s/^#//;
			print OUT "#Sample_ID\t$_\n";next;
		}
		my @tmp=split;
		if(exists $sample{$tmp[0]}){
			print OUT "$sample{$tmp[0]}\t$_\n";
		}else{
			print "Check Whether Your CustomerID is correct or not !!!\n";
		}
		$lncRNA_total=$lncRNA_total+$tmp[2];
		$lncRNA_Q30=$tmp[-1]	if($lncRNA_Q30>$tmp[-1]);
		$lncRNA_data{$tmp[0]}=$tmp[2];
	}
	close(IN);
	close(OUT);
	`mv $in/BMK_1_rawData/BMK_1_Data_Assess/nonRib.AllSample_GC_Q.stat.xls.temp $in/BMK_1_rawData/BMK_1_Data_Assess/nonRib.AllSample_GC_Q.stat.xls`;
	my @lncRNA_data=values %lncRNA_data;
	my $lncRNA_min=&min(\@lncRNA_data);
	$lncRNA_total=sprintf("%.2f",$lncRNA_total/1000000000);
	$lncRNA_min=sprintf("%.2f",$lncRNA_min/1000000000);

	open(OUT,">$temlate_dir/stat.info")||die $!;
	print OUT "#ID\tnum\n";
	if(defined $lnc_dir){
		print "$lnc_dir\n$mRNA_dir\n$struc_dir\n";
		`head -n 4 $lnc_dir/BMK_3_LncRNA_Expression/lncRNA_expression.xls >$temlate_dir/lncRNA_exp.xls`;
		`head -n 4 $mRNA_dir/BMK_2_geneExpression/gene_expression.xls >$temlate_dir/gene_exp.xls`;
		`head -n 4 $mRNA_dir/BMK_1_NewGene/$config{Project_key}.newGene_final.filtered.gff >$temlate_dir/new_gene.xls`;
		`head -n 4 $lnc_dir/BMK_4_LncRNA_Target/Cis_target_gene.xls >$temlate_dir/cis.xls`;
		`head -n 4 $lnc_dir/BMK_4_LncRNA_Target/Trans_target_gene.xls >$temlate_dir/trans.xls`   if(-e "$lnc_dir/BMK_4_LncRNA_Target/Trans_target_gene.xls");
		`head -n 4 $disease_dir/LncRNA_disease.xls >$temlate_dir/lncRNA_disease.xls` if(-e "$disease_dir/LncRNA_disease.xls");
		`head -n 4 $precursor_dir/lncRNA_precursor.xls >$temlate_dir/lncRNA_precursor.xls`;
		
		`head -n 4 $struc_dir/BMK_2_Target_Predict/miRNA-lncRNA/miRNA_target_info.xls >$temlate_dir/miRNA_lncRNA.xls`;
		`head -n 4 $struc_dir/BMK_2_Target_Predict/miRNA-mRNA/miRNA_target_info.xls >$temlate_dir/miRNA_gene.xls`;
		`awk '{print \$1"\\t"\$2}' $mRNA_dir/BMK_1_NewGene/BMK_1_NewGene_Anno/Function_Annotation.stat.xls >$temlate_dir/New_gene_anno.stat`;
		my $new=&cmd("grep '^>' $mRNA_dir/BMK_1_NewGene/$config{Project_key}.newGene.longest_transcript.fa |wc -l");

		`head -n 4 $struc_dir/BMK_3_Gene_Structure_Optimize/$config{Project_key}.geneStructure.optimize.xls >$temlate_dir/geneStructure.optimize.xls`;
		my $optimization=&cmd("wc -l $struc_dir/BMK_3_Gene_Structure_Optimize/$config{Project_key}.geneStructure.optimize.xls") -1;
		###rMATs diff AS
		my @diff_AS = glob "$mRNA_dir/BMK_3_DEG_Analysis/*_vs_*/BMK_3_diff_AS_analysis";
		if(@diff_AS>0){
			`head -n 6 $diff_AS[0]/SE.MATS.JC.xls > $temlate_dir/diff_SE_AS.xls`;
		}
		###TFBS analysis
		my @tfbs_xls = glob "$mRNA_dir/BMK_3_DEG_Analysis/BMK_3_TFBS_Analysis/*Genes_TFBS_predictRes.xls";
		if(@tfbs_xls>0){
			`head -n 6 $tfbs_xls[0] > $temlate_dir/tfbs.xls`;
		}
		###TF_activity analysis
		my @tf_activity = glob "$mRNA_dir/BMK_3_DEG_Analysis/BMK_4_TF_activity/TFs_*.xls";
		if(@tf_activity>0){
			for (my $i=0;$i<@tf_activity;$i++){
		                my $tmp=$tf_activity[$i];
                		if ($tmp=~m/influence/) {
                	        	`head -n 6 $tmp > $temlate_dir/tfs_influence.xls`;
		                }elsif($tmp=~m/cornet/){
                        		`head -n 6 $tmp > $temlate_dir/tfs_cornet.xls`;
                		}else{
                        		`head -n 6 $tmp > $temlate_dir/tfs_activity_grn.xls`;
                		}
        		}
		`cp $mRNA_dir/BMK_3_DEG_Analysis/BMK_4_TF_activity/TFs_influences_heatmap.png $temlate_dir/TFs_influences_heatmap.png`;
		`cp $mRNA_dir/BMK_3_DEG_Analysis/BMK_4_TF_activity/TFs_network.png $temlate_dir/TFs_network.png`;
		}

		my $lncRNA_num=&cmd("wc -l $lnc_dir/BMK_3_LncRNA_Expression/lncRNA_expression.xls") -1;
		my $lncRNA_new_num=&cmd("grep -v '#' $lnc_dir/BMK_2_LncRNA_Prediction/BMK_1_Software_Result/Software_veen.xls |grep -v 'NA'|wc -l");
		my $gene_num=&cmd("wc -l $mRNA_dir/BMK_2_geneExpression/gene_expression.xls") -1;
		my $anno_newgene_num=&cmd("grep -v '^#' $mRNA_dir/BMK_1_NewGene/BMK_1_NewGene_Anno/Integrated_Function.annotation.xls|wc -l");
	
		my ($lnc_map,%map);
		if(defined $circ_dir){
			$lnc_map = "$in/BMK_1_rawData/BMK_2_Mapped_Statistics/lncRNA/All.mappedStat.xls";
			`head -n 4 $struc_dir/BMK_2_Target_Predict/miRNA-circRNA/miRNA_target_info.xls >$temlate_dir/miRNA_circRNA.xls`;
		}else{
			$lnc_map = "$in/BMK_1_rawData/BMK_2_Mapped_Statistics/All.mappedStat.xls";
		}
		
		open(IN,"$lnc_map") or die $!;
		while(<IN>){
			chomp;
			if(/BMK_ID/){next;}
			my @tmp = split;
			if($tmp[2]=~/.*?\((.*?)\)/){
				$map{$tmp[0]}=$1;
			}
		}
		close(IN);
		my @map = values %map;
		my $map_min = &min(\@map);
		my $map_max = &max(\@map);
		
		my $lncRNA_known=$lncRNA_num-$lncRNA_new_num;
		my $gene_known=$gene_num-$new;
		print OUT "gene_new\t$new\n";
		print OUT "anno_newgene_num\t$anno_newgene_num\n";
		print OUT "optimization\t$optimization\n";
		print OUT "lncRNA_Q30\t$lncRNA_Q30\n";
		print OUT "lncRNA_total\t$lncRNA_total\n";
		print OUT "lncRNA_min\t$lncRNA_min\n";
		print OUT "lncRNA_num\t$lncRNA_num\n";
		print OUT "lncRNA_new_num\t$lncRNA_new_num\n";
		print OUT "lncRNA_known_num\t$lncRNA_known\n";
		print OUT "gene_num\t$gene_num\n";
		print OUT "gene_known_num\t$gene_known\n";
		print OUT "lnc_map_min\t$map_min\n";
		print OUT "lnc_map_max\t$map_max\n";
	}
	if(defined $circ_dir){
		`head -n 4 $circ_dir/BMK_2_circRNA_Expression/circRNA_expression.xls >$temlate_dir/circRNA_exp.xls`;
		`head -n 4 $circ_dir/BMK_1_circRNA_Prediction/circRNA_newname.xls >$temlate_dir/circRNA_new.xls`;
		`head -n 4 $circ_dir/BMK_1_circRNA_Prediction/overlap_alitisplice.xls > $temlate_dir/alitisplice.xls`;
		if(exists $config{medical} && $config{medical}=~/GRCh37|GRCh38/){
			`head -n 4 $circ_dir/BMK_1_circRNA_Prediction/BMK_3_Known/circRNA_Circ2disease_annotation.xls > $temlate_dir/circ2disease.xls`;
			`head -n 4 $circ_dir/BMK_1_circRNA_Prediction/BMK_3_Known/circRNA_circBank_annotation.xls > $temlate_dir/circbank.xls`;
		}
		my $circRNA_num=&cmd("wc -l $circ_dir/BMK_2_circRNA_Expression/circRNA_counts.xls") -1;
		my $circRNA_new_num=&cmd("wc -l $circ_dir/BMK_1_circRNA_Prediction/BMK_3_Known/new.list.xls");
		my $circRNA_known_num=$circRNA_num-$circRNA_new_num;

		my ($circ_map,%cmap);
		if(defined $lnc_dir){
			$circ_map = "$in/BMK_1_rawData/BMK_2_Mapped_Statistics/circRNA/All.mappedStat.xls";
		}else{
			$circ_map = "$in/BMK_1_rawData/BMK_2_Mapped_Statistics/All.mappedStat.xls";
			`head -n 4 $in/BMK_3_Target_Predict/miRNA-circRNA/miRNA_target_info.xls >$temlate_dir/miRNA_circRNA.xls`;
			print OUT "lncRNA_total\t$lncRNA_total\n";
			print OUT "lncRNA_Q30\t$lncRNA_Q30\n";
			print OUT "lncRNA_min\t$lncRNA_min\n";
		}
		
		open(IN,"$circ_map") or die $!;
		while(<IN>){
			chomp;
			if(/BMK_ID/){next;}
			my @tmp = split /\t/,$_;
			if($tmp[2]=~/.*?\((.*?)\)/){
				$cmap{$tmp[0]}=$1;
			}
		}
		my @cmap = values %cmap;
		my $cmap_min = &min(\@cmap);
		my $cmap_max = &max(\@cmap);
		
		print OUT "circRNA_num\t$circRNA_num\n";
		print OUT "circRNA_new_num\t$circRNA_new_num\n";
		print OUT "circRNA_known_num\t$circRNA_known_num\n";
		print OUT "circ_map_min\t$cmap_min\n";
		print OUT "circ_map_max\t$cmap_max\n";
	}
	close(OUT);
}	
if(-e "$in/BMK_1_rawData/BMK_1_Data_Assess/All_sample_filter.stat.xls"){
	my %sRNA_data=();
	my $sRNA_total;
	my $sRNA_Q30=100;
	open(IN,"$in/BMK_1_rawData/BMK_1_Data_Assess/All_sample_filter.stat.xls")||die $!;
	open(OUT,">$in/BMK_1_rawData/BMK_1_Data_Assess/All_sample_filter.stat.xls.temp")||die $!;
	while(<IN>){
		chomp;
		if(/^#/){
			s/^#//;
			print OUT "#Sample_ID\t$_\n";next;
		}
		my @tmp=split /\t/,$_;
		if(exists $sample{$tmp[0]}){
			print OUT "$sample{$tmp[0]}\t$_\n";
		}else{
			print "Check Whether Your CustomerID is correct or not !!!\n";
		}
		$sRNA_Q30=$tmp[-1]    if($sRNA_Q30>$tmp[-1]);
		$sRNA_total=$sRNA_total+$tmp[-2];
		$sRNA_data{$tmp[0]}=$tmp[-2];
	}
	close(IN);
	close(OUT);
	`mv $in/BMK_1_rawData/BMK_1_Data_Assess/All_sample_filter.stat.xls.temp $in/BMK_1_rawData/BMK_1_Data_Assess/All_sample_filter.stat.xls`;
	my @sRNA_data=values %sRNA_data;
	my $sRNA_min=&min(\@sRNA_data);
	$sRNA_total=sprintf("%.2f",$sRNA_total/1000000);
	$sRNA_min=sprintf("%.2f",$sRNA_min/1000000);
	$sRNA_dir = (glob "$in/BMK_*_miRNA")[0];
	my $target_dir = (glob "$in/BMK_*_Target_Predict")[0];
	`head -n 4 $sRNA_dir/BMK_2_miRNA_Expression/miRNA_expression.xls >$temlate_dir/miRNA_exp.xls`;
	`head -n 4 $hmdd_dir/HMDD_related_miRNA_with_target_gene.xls > $temlate_dir/hmdd.xls` if (-e "$hmdd_dir/HMDD_related_miRNA_with_target_gene.xls");
	my $mf ="$family_dir/Family_In_miR.xls";
	if(defined $mf){
        	open (IN,"$mf");
	        open (OUT,">$mf.temp");
        	print OUT "#miRNA\tFamily\n";
	        while(<IN>){
        	        chomp;
	                next if(/^#/);
                	print OUT "$_\n";
        	}
	}
	`mv $mf.temp $mf`;
	`head -n 4 $family_dir/Family_In_miR.xls > $temlate_dir/miRNA_family.xls`;
	`head -n 4 $edit_dir/miRNAEdit_cutoff_20.xls > $temlate_dir/miRNAEdit.xls`;
	`head -n 4 $target_dir/miRNA_target_info.xls >$temlate_dir/miRNA_gene.xls`;
	my $miRNA_num=&cmd("grep '^>' $sRNA_dir/BMK_1_miRNA_Prediction/BMK_1_miR_Seq/All_miRNA.expressed.fa |wc -l");
	my $miRNA_known_num=&cmd("grep '^>' $sRNA_dir/BMK_1_miRNA_Prediction/BMK_1_miR_Seq/Known_mature_expressed.fa |wc -l");
	my $miRNA_new_num=&cmd("grep '^>' $sRNA_dir/BMK_1_miRNA_Prediction/BMK_1_miR_Seq/Novel_mature_expressed.fa |wc -l");
	my $target_gene_num=&cmd("grep '^target_gene_num' $target_dir/mir2target.stat.xls|awk '{print \$2}'|cat");
	open(OUT,">$temlate_dir/stat.info")||die $!;
	print OUT "#ID\tnum\n";
	print OUT "sRNA_Q30\t$sRNA_Q30\n";
	print OUT "total_clean_reads\t$sRNA_total\n";
	print OUT "min_clean_reads\t$sRNA_min\n";
	print OUT "miRNA_num\t$miRNA_num\n";
	print OUT "miRNA_known_num\t$miRNA_known_num\n";
	print OUT "miRNA_new_num\t$miRNA_new_num\n";
	print OUT "target_gene_num\t$target_gene_num\n";
	close(OUT);
	print "------extract informations End-----\n\n";
}
}
if(exists $step{2}){
	print "------dataassess html Start-----\n";
	my $cmd = "perl $Bin/build_qc_xml.pl -in $in -o $od/dataassess.xml";
	&run_or_die($cmd);
	print "------dataassess html End-----\n\n";
}

if(exists $step{3}){
	print "------loacl xml Start-----\n";
	if(defined $lnc_dir && defined $mRNA_dir){
		my $cmd="perl $Bin/build_deg_xml.pl -density $lnc_dir/BMK_3_LncRNA_Expression -cfg $detail_cfg -type lncRNA -o $od/lncRNA_local.xml ";
		$cmd.= "-deg $lnc_dir/BMK_5_DEG_Analysis " if($deg_flag==1);
		$cmd.="&& /share/nas2/genome/biosoft/Python/2.7.8/bin/python $Bin/bin/xml2HtmlConverter.py -i $od/lncRNA_local.xml -o $od/lncRNA -n lncRNA ";
		&run_or_die($cmd);

		$cmd="perl $Bin/build_deg_xml.pl -density $mRNA_dir/BMK_2_geneExpression -cfg $detail_cfg -type gene -o $od/gene_local.xml ";
		$cmd.="-deg $mRNA_dir/BMK_3_DEG_Analysis " if($deg_flag==1);
		$cmd.="&& /share/nas2/genome/biosoft/Python/2.7.8/bin/python $Bin/bin/xml2HtmlConverter.py -i $od/gene_local.xml -o $od/gene -n gene ";
		&run_or_die($cmd);
	}
	if(defined $circ_dir){
		$cmd="perl $Bin/build_deg_xml.pl -density $circ_dir/BMK_2_circRNA_Expression -cfg $detail_cfg -type circRNA -o $od/circRNA_local.xml ";
		$cmd.="-deg $circ_dir/BMK_3_DEG_Analysis " if($deg_flag==1);
		$cmd.="&& /share/nas2/genome/biosoft/Python/2.7.8/bin/python $Bin/bin/xml2HtmlConverter.py -i $od/circRNA_local.xml -o $od/circRNA -n circRNA ";
		&run_or_die($cmd);
	}
	if(defined $sRNA_dir){
		$cmd="perl $Bin/build_deg_xml.pl -density $sRNA_dir/BMK_2_miRNA_Expression -cfg $detail_cfg -type miRNA -o $od/miRNA_local.xml ";
		$cmd.="-deg $sRNA_dir/BMK_3_DEG_Analysis " if($deg_flag==1);
		$cmd.="&& /share/nas2/genome/biosoft/Python/2.7.8/bin/python $Bin/bin/xml2HtmlConverter.py -i $od/miRNA_local.xml -o $od/miRNA -n miRNA";
		&run_or_die($cmd);
	}
	print "------loacl xml End-----\n\n";
}

if(exists $step{4}){
	print "------biocloud xml Start-----\n";
	if(defined $lnc_dir && defined $mRNA_dir){
		my $cmd="perl $Bin/build_deg_xml.pl -density $lnc_dir/BMK_3_LncRNA_Expression -cfg $detail_cfg -type lncRNA -o $od/lncRNA.xml -cloud ";
		$cmd.="-deg $lnc_dir/BMK_5_DEG_Analysis " if($deg_flag==1);
		&run_or_die($cmd);
		$cmd="perl $Bin/build_deg_xml.pl -density $mRNA_dir/BMK_2_geneExpression -cfg $detail_cfg -type gene -o $od/gene.xml -cloud ";
		$cmd.="-deg $mRNA_dir/BMK_3_DEG_Analysis " if($deg_flag==1);
		&run_or_die($cmd);
	}
	if(defined $circ_dir){
		$cmd="perl $Bin/build_deg_xml.pl -density $circ_dir/BMK_2_circRNA_Expression -cfg $detail_cfg -type circRNA -o $od/circRNA.xml -cloud ";
		$cmd.="-deg $circ_dir/BMK_3_DEG_Analysis " if($deg_flag==1);
		&run_or_die($cmd);
	}
	if(defined $sRNA_dir){
		$cmd="perl $Bin/build_deg_xml.pl -density $sRNA_dir/BMK_2_miRNA_Expression -cfg $detail_cfg -type miRNA -o $od/miRNA.xml -cloud ";
		$cmd.="-deg $sRNA_dir/BMK_3_DEG_Analysis " if($deg_flag==1);
		&run_or_die($cmd);
	}
	print "------biocloud xml Start-----\n\n";
}
##############################basic sub function###################################
sub cmd{
        my $cmd=shift;
        my $line=`$cmd`;
        $line=(split(/\s+/,$line))[0];
        return $line;
}

sub min{
	my $s=shift;
	my @set=sort{$a <=> $b} @{$s};
	return $set[0];	
}

sub max{
        my $s=shift;
        my @set=sort{$a <=> $b} @{$s};
        return $set[-1];
}

sub relate{
        my $file=shift;
        open(REL,$file)||die $!;
        my %sample=();
        while(<REL>){
                chomp;
                next if($_ !~/^SampleID/);
                my @tmp=split(/\s+/,$_);
        	$sample{$tmp[1]}=$tmp[2];
	}
        close(REL);
        return %sample;
}

sub readConfig{
	my $configFile=shift;
	my $d=Config::General->new(-ConfigFile => "$configFile");
	my %config=$d->getall;	
	return %config;
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


sub USAGE{
	my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: niepy\@biomarker.com.cn
Usage:	
	Options:
	-in	path	analysis input dir, eg: BMK_Result
	-od 	path	analysis output dir, default in/Web_Report/HTML
	-cfg1	file	data_cfg
	-cfg2	file	detail_cfg
	-cloud	null	generate xml for biocloud
	-step	num	default is 1,2,3 ; step4 is for generating xml for biocloud ,if step 4 ,you must "-cloud -step 4 "
	-h	Help

Example: perl $0 -in Web_Report -od Web_Report/HTML -cfg1 data.cfg -cfg2 detail.cfg
	 perl $0 -in Web_Report -od Web_Report/HTML -cfg1 data.cfg -cfg2 detail.cfg -cloud -step 4

USAGE
	print $usage;
	exit;
}


