use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
use Cwd qw(abs_path getcwd);
my ($fpkm,$order,$num,$anno,$od);
GetOptions(
        "h|?"           =>\&USAGE,
	"fpkm:s"	=>\$fpkm,
	"order:s"	=>\$order,
	"num:s"		=>\$num,
	"anno:s"	=>\$anno,
	"od:s"     	=>\$od,
)or &USAGE;
&USAGE unless ($fpkm and $od);

$num ||= 3;

my $head = `grep '^#' $fpkm|cut -f2-`;
chomp($head);
my @sample = split /\t/,$head;
my $sample = join(",",@sample);
$order ||= $sample;

my $gnum;
if($order=~/;/){
	my @group = split /;/,$order;
	$gnum = @group;
}else{
	my @group = split /,/,$order;
	$gnum =@group;
}

my $Rscript="/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript";
if($gnum >= $num){
	`mkdir -p $od`  unless(-d $od);
	`mkdir -p $od/work_sh`  unless(-d "$od/work_sh");
	my $cmd = "perl $Bin/fpkm_prepare.pl -i $fpkm -o $od/kmean_group_fpkm.list -order $order && $Rscript $Bin/kmeans_scale.R --infile $od/kmean_group_fpkm.list --outdir $od";
	open(SH,">$od/work_sh/kmeans.sh") or die $!;
	print SH "$cmd\n";
	$cmd = "perl $Bin/cluster_gene.pl -in $od/kmeans_cluster.txt -out $od/cluster_gene.xls";
	print SH "$cmd\n";
	if(defined $anno){
		$cmd = "perl $Bin/inter_anno.pl -anno $anno -clu $od/kmeans_cluster.txt -o $od/kmeans_cluster_anno.xls";
		print SH "$cmd\n";
	}
	&run_or_die("sh $od/work_sh/kmeans.sh");
}else{
	print "the outdir is $od\n";
	print "the number of group or sample is less than $num, kmeans will not be analysed!\n";
}
#################
sub qsub(){
        my $shfile= shift;
	my $queue="medical.q";
	$queue=$config{Queue_type}	if(defined $config{Queue_type});
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub $shfile --independent --queue $queue ";
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
	-h	Help
	-fpkm	<file>	All_gene_fpkm.list
	-order	<string>	T01,T02,T03;T04,T05,T06;T07,T08,T09 for biological repeat , T01,T02,T03,T04,T05,T06 for no repeat
	-num	<number>	the minimum number of sample or group can be used to analysis,default:3
	-anno	<file>	Integrated_Function.annotation.xls
	-od	<dir>	output dir

Example:
	perl $0 -fpkm All_gene_fpkm.list -od ./kmeans
 
USAGE
	print $usage;
	exit;
}


