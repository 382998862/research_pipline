use Getopt::Long;
use Getopt::Std;
use POSIX;
use Config::General;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path getcwd);
my ($lncRNA_vs_circRNA,$lncRNA_s_mRNA,$mRNA_vs_circRNA,$coexp_dir,$out);

GetOptions(
        "h|?"           =>\&USAGE,
	"lncRNA_vs_circRNA:s"	=>\$lncRNA_vs_circRNA,
	"lncRNA_vs_mRNA:s"	=>\$lncRNA_s_mRNA,
	"mRNA_vs_circRNA:s"	=>\$mRNA_vs_circRNA,
	"in:s"		=>\$coexp_dir,
	"out:s"     	=>\$out,
)or &USAGE;
&USAGE unless ($out);

$out = abs_path($out);
$coexp_dir = abs_path($coexp_dir);

my (@file,$lc_coexp,$lm_coexp,$mc_coexp);
if(defined $coexp_dir){
	$lc_coexp = "$coexp_dir/lncRNA_circRNA/Sig_co_expression.txt";
	$lm_coexp = "$coexp_dir/lncRNA_mRNA/Sig_co_expression.txt";
	$mc_coexp = "$coexp_dir/mRNA_circRNA/Sig_co_expression.txt";
	push @file,$lc_coexp if(-e $lc_coexp);
	push @file,$lm_coexp if(-e $lm_coexp);
	push @file,$mc_coexp if(-e $mc_coexp);
}elsif(defined $lncRNA_vs_circRNA ||defined $lncRNA_s_mRNA ||defined $mRNA_vs_circRNA){
	$lc_coexp = $lncRNA_vs_circRNA if(defined $lncRNA_vs_circRNA);
	$lm_coexp = $lncRNA_s_mRNA if(defined $lncRNA_s_mRNA);
	$mc_coexp = $mRNA_vs_circRNA if(defined $mRNA_vs_circRNA);

	push @file,$lc_coexp if(defined $lncRNA_vs_circRNA);
	push @file,$lm_coexp if(defined $lncRNA_s_mRNA);
	push @file,$mc_coexp if(defined $mRNA_vs_circRNA);
}else{
	&show_log("Check Your Input!");
	die;
}

my $cmd = "cat ";
foreach my $f(@file){
	$cmd .= "$f ";
}
$cmd .= "> $out";
print "$cmd\n";
my $flag = system($cmd);
my $od = dirname $out;
if($flag ==0){
	&run_or_die("touch $od/coexpression_finish");
}else{
	&show_log("merge file, please check.");
	die;
}


###########################sub function#############
sub qsub()
{
        my $shfile= shift;
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $shfile --queue $queue ";
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
}

##################################################
sub USAGE{
	my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
	Options:
        -lncRNA_vs_circRNA	<fiile>	lncRNA_circRNA coexpression.sig.xls
        -lncRNA_vs_mRNA		<file>	lncRNA_mRNA coexpression.sig.xls
        -mRNA_vs_circRNA	<file>	mRNA_circRNA coexpression.sig.xls
	-in	<path>		coexp ananlysis dir coexpression/
        -out		<path>		output path
	-h	Help

Example: perl $0 -lncRNA_vs_circRNA lncRNA_circRNA/Sig_co_expression.txt -lncRNA_vs_mRNA lncRNA_mRNA/Sig_co_expression.txt -mRNA_vs_circRNA mRNA_circRNA/Sig_co_expression.txt -out All_coexpression_sig.xls
         perl $0 -in coexpression/ -out coexpression/All_coexpression_sig.xls

USAGE
	print $usage;
	exit;
}


