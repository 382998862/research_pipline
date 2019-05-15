use Getopt::Long;
use Getopt::Std;
use POSIX;
use Config::General;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);
my ($lncRNA,$miRNA,$mRNA,$circRNA,$diff,$od,$cor,$pvalue,$method,$queue);
GetOptions(
        "h|?"           =>\&USAGE,
	"lncRNA:s"	=>\$lncRNA,
	"mRNA:s"	=>\$mRNA,
	"circRNA:s"	=>\$circRNA,
	"miRNA:s"	=>\$miRNA,
	"lnc_mRNA:s"	=>\$lnc_mRNA,
	"diff:s"	=>\$diff,
	
	"cor:s"		=>\$cor,
	"pvalue:s"	=>\$pvalue,
	"method:s"	=>\$method,
	"queue:s"	=>\$queue,

	"od:s"     	=>\$od,
)or &USAGE;
&USAGE unless ($od);

`mkdir -p $od`	unless(-d $od);
$od=abs_path($od);

$pvalue||=0.01;
$cor||=0.9;
$method||="pearson";
$queue||="medical.q";
my %rnalist=();
$rnalist{mRNA}=abs_path($mRNA)		if(defined $mRNA);
$rnalist{miRNA}=abs_path($miRNA)	if(defined $miRNA);
$rnalist{lncRNA}=abs_path($lncRNA)	if(defined $lncRNA);
$rnalist{circRNA}=abs_path($circRNA)	if(defined $circRNA);

my $Rscript="/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript";
foreach my $type(keys %rnalist){
	`cp $rnalist{$type} $od/$type.fpkm`;
	$rnalist{$type}="$od/$type.fpkm";
}


$diff=abs_path($diff)		if(defined $diff);

my %list=();
if(defined $diff){##如果定义则只对列表下的RNA进行coexpression分析，如果没有定义则对所有的RNA
	$diff=abs_path($diff);
	open(DIFF,$diff)||die $!;
	while(<DIFF>){
		chomp;next if($_=~/^#/);
		my @tmp=split(/\t/,$_);
		$list{$tmp[0]}++;
	}
	close(DIFF);
}

&coexp($rnalist{circRNA},$rnalist{miRNA},"$od/circRNA_miRNA")	if(exists $rnalist{circRNA} && exists $rnalist{miRNA});
&coexp($rnalist{lncRNA},$rnalist{circRNA},"$od/lncRNA_circRNA")	if(exists $rnalist{circRNA} && exists $rnalist{lncRNA});
&coexp($rnalist{mRNA},$rnalist{circRNA},"$od/mRNA_circRNA") if(exists $rnalist{circRNA} && exists $rnalist{mRNA});

if(defined $lnc_mRNA){
	`mkdir $od/lncRNA_mRNA && cp $lnc_mRNA $od/lncRNA_mRNA/Total_co_expression.txt`;
	&run_or_die("perl $Bin/get_sig.pl $od/lncRNA_mRNA/Total_co_expression.txt $od/lncRNA_mRNA/Sig_co_expression.txt $cor $pvalue");
	&diff_subNet("$od/lncRNA_mRNA/Sig_co_expression.txt","$od/lncRNA_mRNA/Diff_Sig_co_expression.txt")      if(defined $diff);
}

###########sub function
sub coexp{
	my($file1,$file2,$od)=@_;
	`mkdir $od`	unless(-d $od);
	`mkdir $od/tmp`	unless(-d "$od/tmp");
	my $head1=`head -n 1 $file1`;chomp($head1);
	`sed -i '1d' $file1`;
	my $head2=`head -n 1 $file2`;chomp($head2);
        `sed -i '1d' $file2`;
	`split -l 500 $file1 -a 4 -d $od/tmp/file1_ && sed -i '1i$head1'  $od/tmp/file1_*`;
	`split -l 500 $file2 -a 4 -d $od/tmp/file2_ && sed -i '1i$head2'  $od/tmp/file2_*`;
	my @num1=glob("$od/tmp/file1_*");
	my @num2=glob("$od/tmp/file2_*");
	open(SH,">$od/tmp/coexp.sh")||die $!;
	for(my $i=0;$i<@num1;$i++){
		for(my $j=0;$j<@num2;$j++){
			print SH "$Rscript $Bin/co_expression.r --rna1 $num1[$i] --rna2 $num2[$j] --outpath $od/tmp --key $i.$j --method $method\n";	
		}
	}	
	close(SH);
	&qsub("$od/tmp/coexp.sh");	
	&run_or_die("sed -i '1d' $od/tmp/co_expression.*.Sig.xls && cat $od/tmp/co_expression.*.Sig.xls  > $od/Total_co_expression.txt && sed -i \"1i\#RNA1\\tRNA2\\tcoefficient\\tpvalue\" $od/Total_co_expression.txt");

	&run_or_die("perl $Bin/get_sig.pl $od/Total_co_expression.txt $od/Sig_co_expression.txt $cor $pvalue");
	&diff_subNet("$od/Sig_co_expression.txt","$od/Diff_Sig_co_expression.txt")	if(defined $diff);
}
	
sub diff_subNet{
	my ($file,$out,$type1,$type2)=@_;
	open(FILE,$file)||die $!;
	open(OUT,">$out")||die $!;
	print OUT "#RNA1\tRNA2\tcoefficient\tpvalue\n";
	while(<FILE>){
		chomp;my ($rna1,$rna2,$co,$p)=split(/\t/,$_);
		if(exists $list{$rna1} || exists $list{$rna2}){
			print OUT "$rna1\t$rna2\t$co\t$p\n";
		}
	}
	close(OUT);
	close(FILE);
}



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
        -lncRNA		<file>		lncRNA fpkm file
        -mRNA		<file>		mRNA fpkm file
        -circRNA	<file>		circRNA RPM file
        -miRNA		<file>		miRNA	RPM file
	-lnc_mRNA	<file>		if lncRNA mRNA coexpression exists(LncRNA Trans predict contained this file),give Total_coexp.file
        -diff		<file>		diffRNA.list, not must, the first column is RNA id
					
        -od		<path>		output path

	-cor		<float>		default 0.9
	-pvalue		<float>		default 0.01
	-method		<str>		default pearson

	-h	Help

Example:

USAGE
	print $usage;
	exit;
}


