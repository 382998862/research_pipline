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

foreach my $type(keys %rnalist){
	&filter_fpkm($rnalist{$type},"$od/$type.fpkm");
	$rnalist{$type}="$od/$type.fpkm";
}


$diff=abs_path($diff)		if(defined $diff);

my %list=();
if(defined $diff){##如果定义则只对列表下的RNA进行coexpression分析，如果没有定义则对所有的RNA
	$diff=abs_path($diff);
	open(DIFF,$diff)||die $!;
	while(<DIFF>){
		chomp;
		my ($type,$rna,$regulate)=split(/\t/,$_);
		$list{$type}{$rna}++;
	}
	close(DIFF);
}

if(exists $rnalist{circRNA} && exists $rnalist{miRNA}){
	mkdir "$od/circRNA_miRNA"	unless(-d "$od/circRNA_miRNA");
	###后期可优化　　若数目太多则进行拆分共表达分析，后合并结果
	print "$rnalist{circRNA}\n";
	my $circ_head=`head -n 1 $rnalist{circRNA}`;chomp($circ_head);
	print "circ: $circ_head\n";
	`sed -i '1d' $rnalist{circRNA}`;
	my $mi_head=`head -n 1 $rnalist{miRNA}`;chomp($mi_head);
        `sed -i '1d' $rnalist{miRNA}`;
	`split -l 500 $rnalist{circRNA} -a 4 -d $od/circRNA_miRNA/circRNA_ && sed -i '1i$circ_head'  $od/circRNA_miRNA/circRNA_*`;
	`split -l 500 $rnalist{miRNA} -a 4 -d $od/circRNA_miRNA/miRNA_ && sed -i '1i$mi_head' $od/circRNA_miRNA/miRNA_*`;
	my @num1=glob("$od/circRNA_miRNA/circRNA_*");
	my @num2=glob("$od/circRNA_miRNA/miRNA_*");
	open(SH,">$od/circRNA_miRNA/coexp.sh")||die $!;
	for(my $i=0;$i<@num1;$i++){
		for(my $j=0;$j<@num2;$j++){
			print SH "Rscript $Bin/co_expression.r --rna1 $num1[$i] --rna2 $num2[$j] --outpath $od/circRNA_miRNA --key circRNA_miRNA_$i.$j --method $method\n";	
		}
	}	
	close(SH);
	&qsub("$od/circRNA_miRNA/coexp.sh");	
	&run_or_die("sed -i '1d' $od/circRNA_miRNA/co_expression.circRNA_miRNA*.Sig.xls && cat $od/circRNA_miRNA/co_expression.circRNA_miRNA*.Sig.xls  > $od/co_expression.miRNA_circRNA.Total.xls && sed -i \"1i\#RNA1\\tRNA2\\tcoefficient\\tpvalue\" $od/co_expression.miRNA_circRNA.Total.xls")	unless(-e "$od/co_expression.miRNA_circRNA.Total.xls");

}
	
if(exists $rnalist{lncRNA} && exists $rnalist{mRNA}){
	mkdir "$od/lncRNA_mRNA"       unless(-d "$od/lncRNA_mRNA");
	`cd $od/lncRNA_mRNA`;
	my $lnc_head=`head -n 1 $rnalist{lncRNA}`;chomp($lnc_head);
	`sed -i '1d' $rnalist{lncRNA}`;
	my $m_head=`head -n 1 $rnalist{mRNA}`;chomp($m_head);
	`sed -i '1d' $rnalist{mRNA}`;
	`split -l 500 $rnalist{lncRNA} -a 4 -d $od/lncRNA_mRNA/lncRNA_ && sed -i '1i$lnc_head'  $od/lncRNA_mRNA/lncRNA_*`;
	`split -l 500 $rnalist{mRNA} -a 4 -d $od/lncRNA_mRNA/mRNA_ && sed -i '1i$m_head' $od/lncRNA_mRNA/mRNA_*`;
	my @num1=glob("$od/lncRNA_mRNA/lncRNA_*");
	my @num2=glob("$od/lncRNA_mRNA/mRNA_*");
	open(SH,">$od/lncRNA_mRNA/coexp.sh")||die $!;
	for(my $i=0;$i<@num1;$i++){
		for(my $j=0;$j<@num2;$j++){
			print SH "Rscript $Bin/co_expression.r --rna1 $num1[$i] --rna2 $num2[$j] --outpath $od/lncRNA_mRNA --key lncRNA_mRNA_$i.$j --method $method\n";
		}
	}
	close(SH);
	&qsub("$od/lncRNA_mRNA/coexp.sh");
	&run_or_die("sed -i '1d' $od/lncRNA_mRNA/co_expression.lncRNA_mRNA*.Sig.xls && cat $od/lncRNA_mRNA/co_expression.lncRNA_mRNA*.Sig.xls >$od/co_expression.mRNA_lncRNA.Total.xls && sed -i \"1i\#RNA1\\tRNA2\\tcoefficient\\tpvalue\" $od/co_expression.mRNA_lncRNA.Total.xls")	unless(-e "$od/co_expression.mRNA_lncRNA.Total.xls");
	&getSig("$od/co_expression.mRNA_lncRNA.Total.xls","$od/Trans_target_gene.xls",$cor ,$pvalue);
}


sub filter_fpkm{
	my ($fpkm,$out)=@_;
	open(FPKM,$fpkm)||die $!;
	open(OUT,">$out")||die $!;
	while(<FPKM>){
		chomp;
		next if($_=~/^$/);
		if($_=~/^#|^ID/){
			print OUT "$_\n";
		}else{
			my @tmp=split(/\t/,$_);
			my $id=shift @tmp;
			my $sum=0;
			foreach my $s(@tmp){$sum+=$s;}	
			if($sum!=0){
				print OUT "$id\t",join("\t",@tmp),"\n";
			}
		}
	}
	close(OUT);
	close(FPKM);
}
sub getSig{
	my ($input,$output,$cor,$pvalue)=@_;
	open(IN,$input)||die $!;
	my %lncs=();
	while(<IN>){
        	chomp;
	        next if($_=~/^#|RNA1/);
        	my ($lnc,$gene,$c,$p)=split(/\t/,$_);
	        next if($p>$pvalue);
        	next if(abs($c)<$cor);
	        $lncs{$lnc}++;
        	push @{$lnc},$gene;
	}
	close(IN);
	open(OUT,">$output")||die $!;
	print OUT "#ID\tTrans_target_gene\n";
	foreach my $lnc(keys %lncs){
        	print OUT "$lnc\t",join("\;",@{$lnc}),"\n";
	}
	close(OUT);
}



sub qsub{
        my $shfile= shift;
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent --queue $queue $shfile";
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
        -diff		<file>		diffRNA.list, not must, if not exists,will do co-expression for all RNA
					2 column: BRCA1	mRNA
        -od		<path>		output path

	-cor		<float>		default 0.9
	-pvalue		<float>		default 0.01
	-method		<str>		default pearson
	-queue		<>		default medical.q
	-h	Help

Example:

USAGE
	print $usage;
	exit;
}


