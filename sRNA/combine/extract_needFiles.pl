#!/usr/bin/perl -w
use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path getcwd);
my ($idir,$data_cfg,$detail_cfg,$od,$all,$only_lncRNA,$only_circRNA,$only_sRNA);
GetOptions(
        "h|?"           =>\&USAGE,
	"idir:s"	=>\$idir,
	"data_cfg:s"	=>\$data_cfg,
	"detail_cfg:s"	=>\$detail_cfg,
	"all"		=>\$all,
	"lncRNA"	=>\$only_lncRNA,
	"circRNA"	=>\$only_circRNA,
	"sRNA"		=>\$only_sRNA,
	"od:s"		=>\$od,
)or &USAGE;
&USAGE unless ($idir || $data_cfg and $detail_cfg and $od and ($all || $only_lncRNA || $only_circRNA || $only_sRNA));

print "Prepare files for personal analysis!\n";

$idir = abs_path($idir);
`mkdir -p $od` unless(-d "$od");
$od = abs_path($od);

my (%config,%type);
my $sample_num=0;
&Creat_cfg("$data_cfg","$detail_cfg","$od/combine.cfg");
###for circos
if(exists $config{medical}){
	&cmd_call("cp $Bin/tools/circos/$config{medical}.txt $od/ideogram.info");
}else{
	my $len;
	if(defined $all || defined $only_lncRNA){
		$len = (glob("$idir/Basic_Analysis/Hisat_Stringtie/Ref_Genome/*.fai"))[0];
	}elsif(defined $only_circRNA){
		$len = (glob("$idir/Basic_Analysis/circRNA_Bwa/Ref_Genome/*.fai"))[0];	
	}elsif(defined $only_sRNA){
		$len = (glob("$idir/Basic_Analysis/sRNA_Alignment/Ref_Database/genome_len.txt"))[0];
	}
	&cmd_call("perl $Bin/tools/get_ideogram.pl $len $od/ideogram.info");
}

if(!-e "$od/ideogram.info"){
	&cmd_call("perl $Bin/tools/ref_GC_len.pl -ref $config{Ref_seq} -od $od");
	my $length = (glob "$od/*.len")[0];
	&cmd_call("perl $Bin/tools/get_ideogram.pl $length $od/ideogram.info");
}else{
	print "Your karyotype file is ok!\n";
}

###for symbol
my $genome_dir = dirname $config{Ref_seq};
my $list1 = "$genome_dir/id_name.list";
my $list2 = "$genome_dir/lncRNA.id_name.list";
if(exists $config{Symbol}){
	&cmd_call("cp $config{Symbol} $od/all.id_name.list");
}else{
	if(-e $list1 && -e $list2){
		&cmd_call("cat $list1 $list2 > $od/all.id_name.list");
	}elsif(-e $list1 && !-e $list2){
		&cmd_call("cp $list1 $od/all.id_name.list");
	}else{
		print "You have no symbol file!\n";
	}
}


my $deg_dir= "$idir/DEG_Analysis";
my $target_dir;

open(All,">$od/All_DEG.xls") or die $!;
print All "#ID\tregulated\n";
foreach my $rna (keys %type){
	my @deg = glob "$deg_dir/$rna/*_vs_*/*.DEG_final.xls";
	foreach my $file(@deg){
		my $name = basename $file;
		`cp $file $od/$rna.$name`;
		open(IN,"$file") or die $!;
		while(<IN>){
			chomp; next if(/^#/);
			my @tmp = split /\t/,$_;
			print All "$tmp[0]\t$tmp[-1]\n";
		}
		close(IN);
	}
}
close(All);
if(defined $all){
	$target_dir = "$idir/miRNA_Target";
	`cp $idir/Basic_Analysis/LncRNA_Analysis/Lnc_target_predict/Cis_target_gene.xls $od/Cis_target_gene.xls`;
	`cp $idir/Basic_Analysis/LncRNA_Analysis/Lnc_target_predict/Trans/Trans_target_gene.xls $od/Trans_target_gene.xls` if($sample_num>=5);
	`cp $idir/Basic_Analysis/circRNA_analysis/new_name/Circ_source_gene.xls $od/circRNA.source.xls`;
	`cp $target_dir/lncRNA/mir2target.list $od/lncRNA.mir2target.list`;
	`cp $target_dir/circRNA/mir2target.list $od/circRNA.mir2target.list`;
	`cp $target_dir/gene/mir2target.list $od/gene.mir2target.list`;
}else{
	$target_dir = "$idir/Target_Predict";
	if(defined $only_lncRNA && defined $only_circRNA){
		`cp $idir/Basic_Analysis/LncRNA_Analysis/Lnc_target_predict/Cis_target_gene.xls $od/Cis_target_gene.xls`;
		`cp $idir/Basic_Analysis/LncRNA_Analysis/Lnc_target_predict/Trans/Trans_target_gene.xls $od/Trans_target_gene.xls` if($sample_num >=5);
		`cp $idir/Basic_Analysis/circRNA_analysis/new_name/Circ_source_gene.xls $od/circRNA.source.xls`;
		`cp $target_dir/miRNA-circRNA/mir2target.list $od/circRNA.mir2target.list`;
		`cp $target_dir/miRNA-lncRNA/mir2target.list $od/lncRNA.mir2target.list`;
		`cp $target_dir/miRNA-mRNA/mir2target.list $od/gene.mir2target.list`;
	}elsif(defined $only_lncRNA && !defined $only_circRNA){
		`cp $idir/Basic_Analysis/LncRNA_Analysis/Lnc_target_predict/Cis_target_gene.xls $od/Cis_target_gene.xls`;
		`cp $idir/Basic_Analysis/LncRNA_Analysis/Lnc_target_predict/Trans/Trans_target_gene.xls $od/Trans_target_gene.xls` if($sample_num >=5);
		`cp $target_dir/miRNA-lncRNA/mir2target.list $od/lncRNA.mir2target.list`;
		`cp $target_dir/miRNA-mRNA/mir2target.list $od/gene.mir2target.list`;
	}elsif(!defined $only_lncRNA && defined $only_circRNA){
		`cp $idir/Basic_Analysis/circRNA_analysis/new_name/Circ_source_gene.xls $od/circRNA.source.xls`;
		`cp $target_dir/miRNA-circRNA/mir2target.list $od/circRNA.mir2target.list`;
		`cp $target_dir/miRNA-mRNA/mir2target.list $od/gene.mir2target.list`;
	}elsif(defined $only_sRNA){
		`cp $target_dir/mir2target.list $od/gene.mir2target.list`;
	}
}


###########################	sub function	##############################
#
sub Creat_cfg{
	my $cfg1 = shift;
	my $cfg2 = shift;
	my $outcfg = shift;
	open(OUT,">$outcfg") or die $!;
	if(defined $all){
		print OUT "Sample\tID\tlncRNA\tcircRNA\tsRNA\n";
		$type{"lncRNA"}=1;$type{"circRNA"}=1;$type{"sRNA"}=1;$type{"gene"}=1;
	}elsif(defined $only_lncRNA && defined $only_circRNA && !defined $only_sRNA){
		print OUT "Sample\tID\tlncRNA\tcircRNA\n";
		$type{"lncRNA"}=1;$type{"circRNA"}=1;$type{"gene"}=1;
	}elsif(defined $only_lncRNA && !defined $only_circRNA && !defined $only_sRNA){
		print OUT "Sample\tID\tlncRNA\n";
		$type{"lncRNA"}=1;$type{"gene"}=1;
	}elsif(!defined $only_lncRNA && defined $only_circRNA && !defined $only_sRNA){
		print OUT "Sample\tID\tcircRNA\n";
		$type{"circRNA"}=1;
	}elsif(!defined $only_lncRNA && !defined $only_circRNA && defined $only_sRNA){
		print OUT "Sample\tID\tsRNA\n";
		$type{"sRNA"}=1;
	}else{
		print "which type of RNA?";die;
	}
	open(CFG1,"$cfg1")||die $!;
	while(<CFG1>){
        	chomp;
		next if(/^#|^$|^\s$/);
		if(/^SampleID/){
                	my @tmp = split /\s+/,$_;
			$sample_num++;
			if(defined $all){
				print OUT "Sample\t$tmp[4]\t$tmp[3]\t$tmp[3]\t$tmp[2]\n";
			}elsif(defined $only_lncRNA && defined $only_circRNA && !defined $only_sRNA){
				print OUT "Sample\t$tmp[1]\t$tmp[1]\t$tmp[1]\n";
			}else{
				print OUT "Sample\t$tmp[1]$tmp[1]\n";
			}
		}
	}
	close(CFG1);
	print OUT "##########parameters###########\n";
	
	open(CFG2,"$cfg2") or die $!;
	while(<CFG2>){
		next if(/^#|^$|^\s$/);
		my @tmp1 = split /\s+/,$_;
		if($tmp1[0] eq "medical" || $tmp1[0] eq "Ref_seq" || $tmp1[0] eq "Symbol" || $tmp1[0] eq "label"){
			$config{$tmp1[0]}=$tmp1[1];
			print OUT "$tmp1[0]\t$tmp1[1]\n";
		}
		if($tmp1[0] eq "coexp_cor" || $tmp1[0] eq "coexp_p" || $tmp1[0] eq "coexp_method"){
			print OUT "$tmp1[0]\t$tmp1[1]\n";
		}
		if($tmp1[0] eq "ceRNA_common" || $tmp1[0] eq "ceRNA_fdr" || $tmp1[0] eq "ceRNA_pvalue"){
			print OUT "$tmp1[0]\t$tmp1[1]\n";
		}
		if($tmp1[0] eq "keyratio" || $tmp1[0] eq "keymin" || $tmp1[0]=~/Queue_type/){
			print OUT "$tmp1[0]\t$tmp1[1]\n";
		}
		if($tmp1[0] eq "Sep"){
			$tmp1[1]=~s/;/_vs_/g;
			$tmp1[1]=~s/,/_/g;
			print OUT "Diff\t$tmp1[1]\n";
		}
		if($tmp1[0] eq "Com"){
			$tmp1[1]=~s/,/_vs_/g;
			print OUT "Diff\t$tmp1[1]\n";
		}
	}
	close(CFG2);
	close(OUT);
}

sub cmd_call {
	print "@_\n";
	system(@_) == 0 or print "Error: system @_ failed: $?\n";
}




sub USAGE{
	my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
	Options:
	-idir		<path>	Analysis/
	-data_cfg	<file>	
	-detail_cfg	<file>
	-lncRNA         
	-mRNA           
	-sRNA           
	-circRNA
	-od		<path>	output file,forced

	-h	Help

Example: perl $0 -idir Analysis -data_cfg data.cfg -detail_cfg detail.cfg -od Input -all
	 perl $0 -idir Analysis -data_cfg data.cfg -detail_cfg detail.cfg -od Input -lncRNA
	 perl $0 -idir Analysis -data_cfg data.cfg -detail_cfg detail.cfg -od Input -lncRNA -circRNA
	 perl $0 -idir Analysis -data_cfg data.cfg -detail_cfg detail.cfg -od Input -circRNA
	 perl $0 -idir Analysis -data_cfg data.cfg -detail_cfg detail.cfg -od Input -sRNA

USAGE
	print $usage;
	exit;
}

