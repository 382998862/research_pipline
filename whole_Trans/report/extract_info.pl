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
$od	=abs_path($od);
$in	=abs_path($in);
$data_cfg	=abs_path($data_cfg);
$detail_cfg	=abs_path($detail_cfg);

`mkdir $od/template` unless (-d "$od/template");
my $temlate_dir = "$od/template";
`cp $Bin/readme/wholetrans_materials_methods.pdf $in/BMK_1_rawData/materials_and_methods.pdf`;

my %config=&readConfig($detail_cfg);

my %step;
$step ||= "1,2,3";
my @step = split /,/,$step;
foreach my $s (@step){
	$step{$s}=1;
}

if(!exists $step{4} && defined $cloud){
	print "Step 4 is for generate xml for cloud ,if you want to do this ,please set you options if '-step 4 -cloud'\n";
	die;
}

if(exists $step{1}){
	print "\n------extract informations Start-----\n";
`cp -r $Bin/images/* $temlate_dir`;
`cp -r $Bin/demo/*png $temlate_dir`;
`cp -r $Bin/demo $temlate_dir`;

my %sample=&relate($data_cfg);
open(RELATE,">$temlate_dir/relate.info")||die $!;
print RELATE "#客户ID\tBMK_ID\tlncRNA_ID\tmRNA_ID\tcircRNA_ID\tmiRNA_ID\n";
foreach my $s(sort{$a cmp $b}keys %sample){
	print RELATE "$sample{$s}{user}\t$s\t$sample{$s}{lncRNA}\t$sample{$s}{mRNA}\t$sample{$s}{circRNA}\t$sample{$s}{sRNA}\n";
}
close(RELATE);

##BMK_1_rawData/BMK_1_Data_Assess/nonRib.AllSample_GC_Q.stat
open(IN,"$in/BMK_1_rawData/BMK_1_Data_Assess/nonRib.AllSample_GC_Q.stat.xls")||die $!;
my %lncRNA_data=();
my $lncRNA_Q30=100;
while(<IN>){
	chomp;	next if($_=~/^#/);
	my @tmp=split;
	$lncRNA_total=$lncRNA_total+$tmp[2];
	$lncRNA_Q30=$tmp[-1]	if($lncRNA_Q30>$tmp[-1]);
	$lncRNA_data{$tmp[0]}=$tmp[2];
}
close(IN);
my $sRNA_Q30=100;
open(IN,"$in/BMK_1_rawData/BMK_1_Data_Assess/miRNA.AllSample_GC_Q.stat.xls")||die $!;
while(<IN>){
	chomp;	next if($_=~/^#/);
	my @tmp=split;
	next if($tmp[-1]=~/Q30/);
	$sRNA_Q30=$tmp[-1]    if($sRNA_Q30>$tmp[-1]);
}
close(IN);

my @lncRNA_data=values %lncRNA_data;
my $lncRNA_min=&min(\@lncRNA_data);
$lncRNA_total=sprintf("%.2f",$lncRNA_total/1000000000);
$lncRNA_min=sprintf("%.2f",$lncRNA_min/1000000000);

`head -n 4 $in/BMK_2_LncRNA/BMK_3_LncRNA_Expression/lncRNA_expression.xls >$temlate_dir/lncRNA_exp.xls`;
`head -n 4 $in/BMK_3_mRNA/BMK_2_geneExpression/gene_expression.xls >$temlate_dir/gene_exp.xls`;

if(exists $config{medical} && $config{medical}=~/GRCh37|GRCh38/){
	`head -n 4 $in/BMK_4_circRNA/BMK_1_circRNA_Prediction/BMK_3_Known/circRNA_Circ2disease_annotation.xls > $temlate_dir/circ2disease.xls`;
	`head -n 4 $in/BMK_4_circRNA/BMK_1_circRNA_Prediction/BMK_3_Known/circRNA_circBank_annotation.xls > $temlate_dir/circbank.xls`;
}
`head -n 4 $in/BMK_4_circRNA/BMK_2_circRNA_Expression/circRNA_expression.xls >$temlate_dir/circRNA_exp.xls`;
`head -n 4 $in/BMK_5_miRNA/BMK_2_miRNA_Expression/miRNA_expression.xls >$temlate_dir/miRNA_exp.xls`;
`head -n 4 $in/BMK_5_miRNA/BMK_6_HMDD/HMDD_related_miRNA_with_target_gene.xls > $temlate_dir/hmdd.xls` if (-e "$in/BMK_5_miRNA/BMK_6_HMDD/HMDD_related_miRNA_with_target_gene.xls");
my $mf ="$in/BMK_5_miRNA/BMK_4_miRNA_Family/Family_In_miR.xls";
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
`head -n 4 $in/BMK_5_miRNA/BMK_4_miRNA_Family/Family_In_miR.xls	> $temlate_dir/miRNA_family.xls`;
`head -n 4 $in/BMK_5_miRNA/BMK_5_miRNA_Edit/miRNAEdit_cutoff_20.xls > $temlate_dir/miRNAEdit.xls`;
`head -n 4 $in/BMK_4_circRNA/BMK_1_circRNA_Prediction/circRNA_newname.xls >$temlate_dir/circRNA_new.xls`;
`head -n 4 $in/BMK_6_Structure/BMK_2_miRNA_Target/circRNA.miRNA_target_info.xls >$temlate_dir/miRNA_circRNA.xls`;
`head -n 4 $in/BMK_6_Structure/BMK_2_miRNA_Target/lncRNA.miRNA_target_info.xls >$temlate_dir/miRNA_lncRNA.xls`;
`head -n 4 $in/BMK_6_Structure/BMK_2_miRNA_Target/gene.miRNA_target_info.xls >$temlate_dir/miRNA_gene.xls`;

`head -n 4 $in/BMK_3_mRNA/BMK_1_NewGene/$config{Project_key}.newGene_final.filtered.gff >$temlate_dir/new_gene.xls`;
`head -n 4 $in/BMK_2_LncRNA/BMK_4_LncRNA_Target/Cis_target_gene.xls >$temlate_dir/cis.xls`;
`head -n 4 $in/BMK_2_LncRNA/BMK_4_LncRNA_Target/Trans_target_gene.xls >$temlate_dir/trans.xls`	if(-e "$in/BMK_2_LncRNA/BMK_4_LncRNA_Target/Trans_target_gene.xls");
`head -n 4 $in/BMK_2_LncRNA/BMK_6_Lnc_disease/LncRNA_disease.xls >$temlate_dir/lncRNA_disease.xls` if(-e "$in/BMK_2_LncRNA/BMK_6_Lnc_disease/LncRNA_disease.xls");
`head -n 4 $in/BMK_2_LncRNA/BMK_7_Lnc_precursor/lncRNA_precursor.xls >$temlate_dir/lncRNA_precursor.xls`;

`awk '{print \$1"\\t"\$2}' $in/BMK_3_mRNA/BMK_1_NewGene/BMK_1_NewGene_Anno/Function_Annotation.stat.xls >$temlate_dir/New_gene_anno.stat`;
my $new=&cmd("grep '^>' $in/BMK_3_mRNA/BMK_1_NewGene/$config{Project_key}.newGene.longest_transcript.fa |wc -l");
`head -n 4 $in/BMK_6_Structure/BMK_1_SNP_Analysis/final.snp.anno.gatk.all.list.xls > $temlate_dir/snp.xls`;
my @fusion = glob "$in/BMK_6_Structure/BMK_5_Gene_Fusion/*_FusionReport.xls";
`head -n 4 $fusion[0] > $temlate_dir/fusion.xls`;

my $opt=(glob("$in/BMK_6_Structure/BMK_3_Gene_Structure_Optimize/*.geneStructure.optimize.xls"))[0];
`head -n 4 $opt >$temlate_dir/geneStructure.optimize.xls`;
my $optimization=&cmd("wc -l $opt");
$optimization=$optimization-1;

###rMATs diff AS
my @diff_AS = glob "$in/BMK_3_mRNA/BMK_3_DEG_Analysis/*_vs_*/BMK_3_diff_AS_analysis";
if(@diff_AS>0){
	`head -n 4 $diff_AS[0]/SE.MATS.JC.xls > $temlate_dir/diff_SE_AS.xls`;
}
###TFBS analysis
my @tfbs_xls = glob "$in/BMK_3_mRNA/BMK_3_DEG_Analysis/BMK_3_TFBS_Analysis/*Genes_TFBS_predictRes.xls";
if(@tfbs_xls>0){
	`head -n 6 $tfbs_xls[0] > $temlate_dir/tfbs.xls`;
}
###TF_activity analysis
my @tf_activity = glob "$in/BMK_3_mRNA/BMK_3_DEG_Analysis/BMK_4_TF_activity/TFs_*.xls";
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

my $lncRNA_num=&cmd("wc -l $in/BMK_2_LncRNA/BMK_3_LncRNA_Expression/lncRNA_expression.xls");
$lncRNA_num=$lncRNA_num-1;

my $lncRNA_new_num=&cmd("grep -v '#' $in/BMK_2_LncRNA/BMK_2_LncRNA_Prediction/BMK_1_Software_Result/Software_veen.xls |grep -v 'NA'|wc -l");
my $miRNA_num=&cmd("grep '^>' $in/BMK_5_miRNA/BMK_1_miRNA_Prediction/BMK_1_miR_Seq/All_miRNA.expressed.fa |wc -l");
my $miRNA_known_num=&cmd("grep '^>' $in/BMK_5_miRNA/BMK_1_miRNA_Prediction/BMK_1_miR_Seq/Known_mature_expressed.fa |wc -l");
my $miRNA_new_num=&cmd("grep '^>' $in/BMK_5_miRNA/BMK_1_miRNA_Prediction/BMK_1_miR_Seq/Novel_mature_expressed.fa |wc -l");
my $gene_num=&cmd("wc -l $in/BMK_3_mRNA/BMK_2_geneExpression/gene_expression.xls");
$gene_num=$gene_num-1;
my $circRNA_num=&cmd("wc -l $in/BMK_4_circRNA/BMK_2_circRNA_Expression/circRNA_counts.xls");
$circRNA_num=$circRNA_num-1;
my $circRNA_new_num=&cmd("wc -l $in/BMK_4_circRNA/BMK_1_circRNA_Prediction/BMK_3_Known/new.list.xls");
my $circRNA_known_num=$circRNA_num-$circRNA_new_num;

my $ceRNA_file=(glob "$in/BMK_7_Combine/BMK_3_ceRNA_random/*ceRNA_pair_adjust_p_Sig.xls")[0];
print "$ceRNA_file\n";
`head -n 4 $ceRNA_file>$temlate_dir/ceRNA.xls`;

my @Sub = glob "$in/BMK_7_Combine/BMK_3_ceRNA_random/*/sub_*_network.xls";
&select_sub_network(\@Sub);

`perl $Bin/combine.table.pl -data_cfg $data_cfg -detail_cfg $detail_cfg -in $in/BMK_7_Combine/BMK_4_cytoscape/ -od $temlate_dir/`;

my $ceRNA_relate=&cmd("wc -l $ceRNA_file");
$ceRNA_relate=$ceRNA_relate-1;
my ($ceRNA_node,$ceRNA_lncRNA,$ceRNA_circRNA,$ceRNA_gene);
if(exists $config{medical}){
	$ceRNA_node=`awk '{if(\$1!~\/RNA1|^#\/) {print \$1\"\\n\"\$3}}' $ceRNA_file|sort|uniq|wc -l`;
	$ceRNA_lncRNA=`awk '{print \$1\"\\n\"\$3}' $ceRNA_file|sort|uniq|awk -F \":\" '\$1~\/lncRNA\/{a+=1}END{print a}'|cat`;
	$ceRNA_circRNA=`awk '{print \$1\"\\n\"\$3}' $ceRNA_file|sort|uniq|awk -F \":\" '\$1~\/circRNA\/{b+=1}END{print b}'|cat`;
	$ceRNA_gene=`awk '{print \$1\"\\n\"\$3}' $ceRNA_file|sort|uniq|awk -F \":\" '\$1~\/gene\/{d+=1}END{print d}'|cat`;	
}else{
	$ceRNA_node=`awk '{if(\$1!~\/RNA1|^#\/) {print \$1\"\\n\"\$2}}' $ceRNA_file|sort|uniq|wc -l`;
	$ceRNA_lncRNA=`awk '{print \$1\"\\n\"\$2}' $ceRNA_file|sort|uniq|awk -F \":\" '\$1~\/lncRNA\/{a+=1}END{print a}'|cat`;
	$ceRNA_circRNA=`awk '{print \$1\"\\n\"\$2}' $ceRNA_file|sort|uniq|awk -F \":\" '\$1~\/circRNA\/{b+=1}END{print b}'|cat`;
	$ceRNA_gene=`awk '{print \$1\"\\n\"\$2}' $ceRNA_file|sort|uniq|awk -F \":\" '\$1~\/gene\/{d+=1}END{print d}'|cat`;
}
chomp($celine);chomp($ceRNA_node);chomp($ceRNA_lncRNA);chomp($ceRNA_circRNA);chomp($ceRNA_gene);

my $gene_known=$gene_num-$new;
my $lncRNA_known=$lncRNA_num-$lncRNA_new_num;

open(OUT,">$temlate_dir/stat.info")||die $!;
print OUT "#ID\tnum\n";
print OUT "ceRNA_relate\t$ceRNA_relate\n";
print OUT "ceRNA_node\t$ceRNA_node\n";
print OUT "ceRNA_gene\t$ceRNA_gene\n";
print OUT "ceRNA_circRNA\t$ceRNA_circRNA\n";
print OUT "ceRNA_lncRNA\t$ceRNA_lncRNA\n";
print OUT "optimization\t$optimization\n";
print OUT "lncRNA_Q30\t$lncRNA_Q30\n";
print OUT "sRNA_Q30\t$sRNA_Q30\n";

print OUT "lncRNA_total\t$lncRNA_total\n";
print OUT "lncRNA_min\t$lncRNA_min\n";
print OUT "lncRNA_num\t$lncRNA_num\n";
print OUT "lncRNA_new_num\t$lncRNA_new_num\n";
print OUT "lncRNA_known_num\t$lncRNA_known\n";

print OUT "miRNA_num\t$miRNA_num\n";
print OUT "miRNA_known_num\t$miRNA_known_num\n";
print OUT "miRNA_new_num\t$miRNA_new_num\n";

print OUT "gene_new\t$new\n";
print OUT "gene_num\t$gene_num\n";
print OUT "gene_known_num\t$gene_known\n";

print OUT "circRNA_num\t$circRNA_num\n";
print OUT "circRNA_new_num\t$circRNA_new_num\n";
print OUT "circRNA_known_num\t$circRNA_known_num\n";

close(OUT);
	print "------extract informations End-----\n\n";
}

if(exists $step{2}){
	print "------dataassess html Start-----\n";
	my $cmd = "perl $Bin/build_qc_xml.pl -in $in -o $od/dataassess.xml";
	&run_or_die($cmd);
	print "------dataassess html End-----\n\n";
}

if(exists $step{3}){
	print "------loacl xml Start-----\n";
	my $cmd="perl $Bin/build_deg_xml.pl -deg $in/BMK_2_LncRNA/BMK_5_DEG_Analysis -density $in/BMK_2_LncRNA/BMK_3_LncRNA_Expression -cfg $detail_cfg -type lncRNA -o $od/lncRNA_local.xml && /share/nas2/genome/biosoft/Python/2.7.8/bin/python $Bin/bin/xml2HtmlConverter.py -i $od/lncRNA_local.xml -o $od/lncRNA -n lncRNA";
	&run_or_die($cmd);

	$cmd="perl $Bin/build_deg_xml.pl -deg $in/BMK_3_mRNA/BMK_3_DEG_Analysis -density $in/BMK_3_mRNA/BMK_2_geneExpression -cfg $detail_cfg -type gene -o $od/gene_local.xml && /share/nas2/genome/biosoft/Python/2.7.8/bin/python $Bin/bin/xml2HtmlConverter.py -i $od/gene_local.xml -o $od/gene -n gene";
	&run_or_die($cmd);

	$cmd="perl $Bin/build_deg_xml.pl -deg $in/BMK_4_circRNA/BMK_3_DEG_Analysis -density $in/BMK_4_circRNA/BMK_2_circRNA_Expression -cfg $detail_cfg -type circRNA -o $od/circRNA_local.xml && /share/nas2/genome/biosoft/Python/2.7.8/bin/python $Bin/bin/xml2HtmlConverter.py -i $od/circRNA_local.xml -o $od/circRNA -n circRNA";
	&run_or_die($cmd);

	$cmd="perl $Bin/build_deg_xml.pl -deg $in/BMK_5_miRNA/BMK_3_DEG_Analysis -density $in/BMK_5_miRNA/BMK_2_miRNA_Expression -cfg $detail_cfg -type miRNA -o $od/miRNA_local.xml && /share/nas2/genome/biosoft/Python/2.7.8/bin/python $Bin/bin/xml2HtmlConverter.py -i $od/miRNA_local.xml -o $od/miRNA -n miRNA";
	&run_or_die($cmd);

	$cmd="perl $Bin/build_combine_xml.pl -combine $in/BMK_7_Combine -detail_cfg $detail_cfg -data_cfg $data_cfg -o $od/combine_local.xml -type combine && /share/nas2/genome/biosoft/Python/2.7.8/bin/python $Bin/bin/xml2HtmlConverter.py -i $od/combine_local.xml -o $od/combine -n combine";
	&run_or_die($cmd);
	print "------loacl xml End-----\n\n";
}

if(exists $step{4}){
	print "------biocloud xml Start-----\n";
	if(defined $cloud){
		my $cmd="perl $Bin/build_deg_xml.pl -deg $in/BMK_2_LncRNA/BMK_5_DEG_Analysis -density $in/BMK_2_LncRNA/BMK_3_LncRNA_Expression -cfg $detail_cfg -type lncRNA -o $od/lncRNA.xml -cloud ";
		&run_or_die($cmd);
		$cmd="perl $Bin/build_deg_xml.pl -deg $in/BMK_3_mRNA/BMK_3_DEG_Analysis -density $in/BMK_3_mRNA/BMK_2_geneExpression -cfg $detail_cfg -type gene -o $od/gene.xml -cloud ";
		&run_or_die($cmd);
		$cmd="perl $Bin/build_deg_xml.pl -deg $in/BMK_4_circRNA/BMK_3_DEG_Analysis -density $in/BMK_4_circRNA/BMK_2_circRNA_Expression -cfg $detail_cfg -type circRNA -o $od/circRNA.xml -cloud ";
		&run_or_die($cmd);
		$cmd="perl $Bin/build_deg_xml.pl -deg $in/BMK_5_miRNA/BMK_3_DEG_Analysis -density $in/BMK_5_miRNA/BMK_2_miRNA_Expression -cfg $detail_cfg -type miRNA -o $od/miRNA.xml -cloud ";
		&run_or_die($cmd);
		$cmd="perl $Bin/build_combine_xml.pl -combine $in/BMK_7_Combine -detail_cfg $detail_cfg -data_cfg $data_cfg -o $od/combine.xml -type combine -cloud ";
		&run_or_die($cmd);
	}
	print "------biocloud xml End-----\n\n";
}
##############################basic sub function
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
sub relate{     ##获取不同RNA样品的对应关系
        my $file=shift;
        open(REL,$file)||die $!;
        my %sample=();
        while(<REL>){
                chomp;
                next if($_ !~/^SampleID/);
                my @tmp=split(/\s+/,$_);
                $sample{$tmp[-1]}{sRNA}=$tmp[2];
                $sample{$tmp[-1]}{lncRNA}=$tmp[3];
                $sample{$tmp[-1]}{circRNA}=$tmp[3];
                $sample{$tmp[-1]}{mRNA}=$tmp[3];
		$sample{$tmp[-1]}{user}=$tmp[1];
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
sub select_sub_network{
        my $array = shift;
        my (@num1,@num2,@num3,$flag);
        foreach my $sub (@$array){
                my $line = `wc -l $sub`;
                chomp($line);
                if($line>3){push @num3,$sub;}
                if($line>2){push @num2,$sub;}
                if($line>1){push @num1,$sub;}
                if($line=1){$flag=1;}
        }

        if(scalar(@num3)>=1){
                `head -4 $num3[0] > $temlate_dir/sub_network.xls`;
        }else{
                if(scalar(@num2)>=1){
                        `head -3 $num2[0] > $temlate_dir/sub_network.xls`;
                }else{
                        if(scalar(@num1)>=1){
                                `head -2 $num1[0] > $temlate_dir/sub_network.xls`;
                        }else{
                                if($flag=1){
                                        `cp $Bin/demo/demo.sub_network.xls $temlate_dir/sub_network.xls`;
                                }
                        }
                }
        }
}

sub qsub()
{
        my $shfile= shift;
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $shfile";
        &run_or_die($cmd);              
        return ;
}
sub run_or_die()
{
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
sub show_log()
{
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
	-od 	path	analysis output dir, default in/BMK_Result/HTML
	-cfg1	file	data_cfg
	-cfg2	file	detail_cfg
	-cloud	null	xml for cloud
	-step	num	default is 1,2,3 ; 4 if for generating xml for cloud
	-h	Help

Example: perl $0 -in Web_Report -od Web_Report/HTML -cfg1 data.cfg -cfg2 detail.cfg
	 perl $0 -in Web_Report -od Web_Report/HTML -cfg1 data.cfg -cfg2 detail.cfg -cloud -step 4

USAGE
	print $usage;
	exit;
}


