use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use newPerlBase;
use Cwd qw(abs_path getcwd);
my ($input, $od, $ratio,$anno,$type,$num,$symbol);
GetOptions(
        "h|?"           =>\&USAGE,
	"i:s"		=>\$input,
	"type:s"	=>\$type,
	"symbol:s"	=>\$symbol,
	"ratio:s"	=>\$ratio,
	"num:s"		=>\$num,
	"anno:s"	=>\$anno,
	"od:s"     	=>\$od,
)or &USAGE;
&USAGE unless ($od and $input and $anno and $type);
$ratio||="0.1";
$num||=5;

my $cmd="Rscript $Bin/bin/random_walk.r --network $input --outpath $od --from 1 --to 2 ";
system("$cmd");

my $line=`wc -l $od/Net_score.txt`;chomp($line);
$line=(split(/\s+/,$line))[0];

my $top=int($line*$ratio);
my $top_file="$od/Key_RNA_top$ratio.txt";
`head -n $top $od/Net_score.txt >$top_file`;

##############
if(defined $type){
	my %genes=();
	open(TYPE,$type)||die $!;
	while(<TYPE>){
		chomp;my @tmp=split;
		$genes{$tmp[0]}++	if($tmp[1] eq "gene");				
	}
	close(TYPE);
	open(TOP,$top_file)||die $!;
	open(GENE,">$od/Key_gene.txt")||die $!;
	while(<TOP>){
		chomp;my @tmp=split;
		if(exists $genes{$tmp[0]}){
			print GENE join("\t",@tmp),"\n";
		}
	}
	close(TOP);
	close(GENE);
	$top_file="$od/Key_gene.txt";
}
##############
if(defined $anno){

my $kegg=(glob("$anno/*.Kegg.pathway"))[0];
my $ko=(glob("$anno/*Kegg.ko"))[0];
`mkdir $od/KEGG`	unless(-d "$od/KEGG");
&run_or_die("perl $Bin/../../DEG_Analysis/enrich/bin/process_Kegg.pathway.pl $kegg $od/KEGG/KEGG.info && /share/nas2/genome/biosoft/R/3.3.2/bin/Rscript $Bin/../../DEG_Analysis/enrich/enrich.r --deg $top_file --kegg $od/KEGG/KEGG.info --prefix Key  --od $od/KEGG");

my $paths=&getKO("$od/KEGG/Key_KEGG_pathway_enrich.list",$num);
$cmd ="perl $Bin/KEGG/KEGG_network_construct.pl -gene $ko -od $od/KEGG -ko $paths -key $top_file";
$cmd.=" -symbol $symbol"	if(defined $symbol);
&run_or_die($cmd);
}

##############Function sub
sub getKO{
	my ($file,$num)=@_;
	my $count=0;
	my @ko=();
	open(FILE,$file)||die $!;
	while(<FILE>){
		chomp;next if($_=~/^ID/);
		next if($count >= $num);
		my @tmp=split(/\t/,$_);
		push @ko,$tmp[0];
		$count++;
	}
	close(FILE);
	my $kos=join(",",@ko);
	return $kos;
}


##############
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
}


sub USAGE{
	my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
	Options:
	-i	<file>	input network file, forced
	-type	<file>	input node type file, #ID type, not must
	-symbol	<file>	symbol file, eg: id_name.list, not must
	-ratio	<float>	random ratio, default 0.1
	-num	<int>	how much pathways would be intergrated
	-anno	<dir>	anno dir	
	-od	<str>	output dir, forced

	-h	Help

Example:

USAGE
	print $usage;
	exit;
}


