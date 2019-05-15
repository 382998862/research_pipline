use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path getcwd);
use newPerlBase;
my ($od,$cfg,$deg,$anno,$trans,$cis,$source,$symbol,$ceRNA,$queue);
GetOptions(
        "h|?"           =>\&USAGE,
	"deg:s"		=>\$deg,
	"od:s"		=>\$od,	
	"cfg:s"		=>\$cfg,
	"ceRNA:s"	=>\$ceRNA,
        "anno:s"        =>\$anno,
        "trans:s"       =>\$trans,
        "cis:s"         =>\$cis,
        "source:s"      =>\$source,
        "symbol:s"      =>\$symbol,
	"queue:s"	=>\$queue,
)or &USAGE;
&USAGE unless ($cfg and $od and $deg);

$deg=abs_path($deg);
`mkdir -p $od` unless (-d "$od");
$od=abs_path($od);
$cfg=abs_path($cfg);
$ceRNA=abs_path($ceRNA);

$cis ||= "$deg/Cis_target_gene.xls";
$trans ||= "$deg/Trans_target_gene.xls";
$source ||= "$deg/circRNA.source.xls";
$symbol ||= "$deg/all.id_name.list";

my @files=glob "$deg/*.DEG_final.xls";
#print "@files\n";
my (@RNAs,%samples,%diffs);
my $sample_num=0;
my %config=&readcfg($cfg);

my $ratio=exists $config{keyratio}?$config{keyratio}:0.05;
my $min=exists $config{keymin}?$config{keymin}:100;
$queue = exists $config{Queue_type}?$config{Queue_type}:"medical.q";

`mkdir $od/work_sh` unless(-d "$od/work_sh"); 
open (SH,">$od/work_sh/random_diff.sh") or die $!;

foreach my $diff(@{$config{Diff}}){
	my @vs=split(/_vs_/,$diff);
	my @g1=split(/_/,$vs[0]);
	my @g2=split(/_/,$vs[1]);
	open(DIFF,">$od/All.$diff.xls")||die $!;
	print DIFF "#ID\tregulated\n";
	close(DIFF);
	foreach my $rna(@RNAs){
		my $ng1=join("_",map{$samples{$_}{$rna}} @g1);
		my $ng2=join("_",map{$samples{$_}{$rna}} @g2);
		$diffs{$diff}{$rna}="${ng1}_vs_$ng2";
		my $f=$rna eq "gene"?&search(\@files,$diffs{$diff}{$rna},"gene"):&search(\@files,$diffs{$diff}{$rna},$rna);
		next if(!-e $f);
		`awk '{print \$1\"\\t\"\$NF}' $f|grep -v '#' >>$od/All.$diff.xls`;
	}
	`mkdir $od/$diff` unless (-d "$od/$diff");
	print SH "perl $Bin/random_diff_ceRNA.pl -ceRNA $ceRNA -diff $od/All.$diff.xls -ratio $ratio -num 5 -anno $anno -od $od/$diff -cis $cis -source $source -cfg $cfg ";
	print SH "-trans $trans " if($sample_num >=5);
	print SH "-symbol $symbol " if(exists $config{medical});
	print SH "\n";
}

close(OUT);
&qsub("$od/work_sh/random_diff.sh");
&qsubCheck("$od/work_sh/random_diff.sh");

#################################sub function####################3

sub search{
	my ($f,$key,$type)=@_;
	my @files=@{$f};
	foreach my $f(@files){
		my $base=basename $f;
		return $f	if($base eq "$type.$key.DEG_final.xls");
	}
	return 0;
}

sub readcfg{
	my $cfg=shift;
	open(CFG,$cfg)||die $!;
	my @rnas=();
	while(<CFG>){
		chomp;next if($_=~/^#|^\s+|^$/);
		my @tmp=split(/\s+/,$_);
		if($tmp[0] eq "Diff"){
			push @{$config{Diff}},$tmp[1];
		}elsif($tmp[0] eq "Sample"){
			shift @tmp;
			if($tmp[0] eq "ID" && scalar(@rnas)==0){
				print "$_\n";
				shift @tmp;@rnas=@tmp;
			}else{
				$sample_num++;
				my $id=shift @tmp;
				for(my $i=0;$i<@tmp;$i++){
					$samples{$id}{$rnas[$i]}=$tmp[$i];
					print "$id\t$rnas[$i]\t$samples{$id}{$rnas[$i]}\t$tmp[$i]\n";
				}
				if(exists $samples{$id}{lncRNA}){
					$samples{$id}{gene} = $samples{$id}{lncRNA};
					push @rnas,"gene";
				}
			}
			
		}else{
			$config{$tmp[0]}=$tmp[1];
		}
	}
	close(CFG);
	my %count;
	@RNAs = grep { ++$count{$_} < 2;} @rnas;
	return %config;
}

sub qsub(){
        my $shfile= shift;
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $shfile --queue $queue";
        &run_or_die($cmd);
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
	-deg	<files>		all_deg files
	-od	<dir>		output dir
	-cfg	<file>		config file
	-ceRNA	<file>
	-anno	<file>
	-trans
	-cis
	-source
	-symbol
	-h	Help

Example:

USAGE
	print $usage;
	exit;
}


