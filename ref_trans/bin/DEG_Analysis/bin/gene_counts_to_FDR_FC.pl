use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use newPerlBase;
use Cwd qw(abs_path getcwd);
my ($od,$conf,$count,$filter,$fpkm,$type);
GetOptions(
        "h|?"           =>\&USAGE,
	"od:s"     	=>\$od,
	"count:s"	=>\$count,
	"conf:s"	=>\$conf,
	"filter:s"	=>\$filter,
	"fpkm:s"	=>\$fpkm,
	"type:s"	=>\$type
	
)or &USAGE;
&USAGE unless ($od and $conf and $count);
##################sample corelation



my %script=();
$script{"EBseq"}="$Bin/diff/ebseq_analysis.r";
$script{"DEseq"}="$Bin/diff/deseq_analysis.r";
$script{"edgeR"}="$Bin/diff/edger_analysis.r";
my $Rscript="/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript";


$count=abs_path($count);
$conf=abs_path($conf);
my %config=&readConfig($conf);
my ($Com,$Sep,$PValue,$FDR,$FC);

$PValue=$config{PValue}		if(exists $config{PValue});
$PValue=$config{Pvalue}		if(exists $config{Pvalue});
$FDR=$config{FDR}		if(exists $config{FDR});
$FC=$config{fold}		if(exists $config{fold});
$FC||=2;


my $select="--FC $FC ";##############Diff selection standard
if(!defined $FDR && !defined $PValue){
	$select .="--FDR 0.01";
}elsif(defined $FDR){
	$select .="--FDR $FDR";
}elsif(defined $PValue){
	$select .="--PValue $PValue";
}
### -FC 2 -FDR 0.01
$Com=$config{Method_NR}		if(exists $config{Method_NR});
$Sep=$config{Method_RE}		if(exists $config{Method_RE});

$Com||="EBseq";
$Sep||="DEseq";
$filter||="cpm";

if(defined $PValue && !defined $FDR && $Com eq "EBseq"){
	print "EBseq could not produce PValue! So method for deg group without Biological replicates would be replaced edgeR\n";
	$Com="edgeR";
}


my $flag="";#####contain FPKM info
if(defined $fpkm){
	$fpkm=abs_path($fpkm);
	$type||="FPKM";
	$flag="$fpkm $type";
}

#Method_RE       DEseq           ###DEseq/edgeR
#Method_NR       EBseq           ###EBseq/edgeR

`mkdir $od`	unless(-d $od);
$od=abs_path($od);
mkdir "$od/work_sh" unless (-d "$od/work_sh");

my @coms=();my @seps=();
open(CFG,$conf)||die $!;
while(<CFG>){
	chomp;
	my @tmp=split(/\s+/,$_);
	if($tmp[0] eq "Com"){
		push @coms,$tmp[1];
#		push @degs,$tmp[1];	
	}elsif($tmp[0] eq "Sep"){
		push @seps,$tmp[1];
#		push @degs,$tmp[1];
	}
}
close(CFG);

###############DEG parameter
#"1) read_count.txt: the read count file for RNA-SEQ"
#"2) multi_factor_matrix.txt: the multiple factor matrix"
#"3) out.de: the output of DE"
#"4) filter: the method used for filter low gene expression (cpm/count) "
#"5) expression.list: not must "
#"6) expression.type:(FPKM or TPM) :not must if givn expression.list default FPKM"

open(SH1,">$od/work_sh/s1.select_deg.sh")||die $!;
open(SH2,">$od/work_sh/s2.draw_heatmap_MA_volcano.sh")||die $!;
if(@coms>0){
	foreach my $com(@coms){
		my $vs=&getdecfg($com);				
		my $diff_cmd="$Rscript $script{$Com} $count $od/$vs/$vs.de.config $od/$vs/$vs $filter $flag";
		$diff_cmd .=" 0.16 "	if(exists $config{medical} && $config{medical}=~/GRCh/ && $Com eq "edgeR");
		print SH1 "$diff_cmd && perl $Bin/diff/filter_final_DEG.pl -input $od/$vs/$vs.all -output $od/$vs/$vs.DEG_final.xls $select \n";
		print SH2 "$Rscript $Bin/draw/pheatmap.r --infile $od/$vs/$vs.DEG_final.xls --outfile $od/$vs/$vs.heatmap --scale none --color.type 1 --legend  --is.log  --show_colnames --cluster_rows && $Rscript $Bin/draw/plot_MA_volcano.r --all $od/$vs/$vs.all $select --od $od/$vs/$vs --type $type\n";
	}
}
if(@seps>0){
	foreach my $sep(@seps){
		my $vs=&getdecfg($sep);
		print SH1 "$Rscript $script{$Sep} $count $od/$vs/$vs.de.config $od/$vs/$vs $filter $flag && perl $Bin/diff/filter_final_DEG.pl -input $od/$vs/$vs.all -output $od/$vs/$vs.DEG_final.xls $select \n";
		print SH2 "$Rscript $Bin/draw/pheatmap.r --infile $od/$vs/$vs.DEG_final.xls --outfile $od/$vs/$vs.heatmap --scale none --color.type 1 --legend  --is.log  --show_colnames --cluster_rows && $Rscript $Bin/draw/plot_MA_volcano.r --all $od/$vs/$vs.all $select --od $od/$vs/$vs --type $type\n";
	}
}
close(SH1);
close(SH2);


&qsub("$od/work_sh/s1.select_deg.sh");
&qsub("$od/work_sh/s2.draw_heatmap_MA_volcano.sh");



####################
#/share/nas1/linhj/wenyh/project/BMK171221-I239-01/Analysis_lnc/Lnc_Diff_Analysis/L22_L23_L24_vs_L01_L02_L03/L22_L23_L24_vs_L01_L02_L03.de.config
sub getdecfg{
	my $vs=shift;
	my @g1=();my @g2=();
	if($vs=~/\;/){
		my @tmp=split(/\;/,$vs,2);
		@g1=split(/,/,$tmp[0]);
		@g2=split(/,/,$tmp[1]);
	}else{
		my @tmp=split(/,/,$vs,2);
		push @g1,$tmp[0];
		push @g2,$tmp[1];
	}
	my $base=join("_",@g1)."_vs_".join("_",@g2);

	`mkdir $od/$base`	unless(-d "$od/$base");
	open(DCFG,">$od/$base/$base.de.config")||die $!;
	print DCFG "sample\tprocess\n";	
	foreach my $g(@g1){print DCFG "$g\tone\n";}
	foreach my $g(@g2){print DCFG "$g\ttwo\n";}
	close(DCFG);
	return $base;
}

#################
sub readConfig{
	my $configFile=shift;
	my $d=Config::General->new(-ConfigFile => "$configFile");
	my %config=$d->getall;	
	return %config;
}
sub qsub{
        my $shfile= shift;
	my $queue="medical.q";
	$queue=$config{Queue_type}	if(exists $config{Queue_type});
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub $shfile --independent --queue $queue ";
        &run_or_die($cmd);              
        return ;
}
sub run_or_die{
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
sub show_log{
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
Contact: wenyh\@biomarker.com.cn
Usage:
	
	Options:
	-h	Help
	-count	file	forced		All_gene_counts.list
	-conf	file	forced		DEG parameter(detail.cfg)
        -od	path	forced		output path
	-filter	str	not must	default cpm, (cpm/count)
	-fpkm	file	not must	All_gene_expression.list (can All_gene_fpkm.list or All_gene_tmp.list)
	-type	str	not must	All_gene_expression.list type can FPKM/TPM
					Or any type defined by yourself, if defined -fpkm, default FPKM
Example:

USAGE
	print $usage;
	exit;
}


