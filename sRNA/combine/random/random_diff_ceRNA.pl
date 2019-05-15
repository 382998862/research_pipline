use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path getcwd);
my ($ceRNA, $diff, $cfg, $od, $trans, $cis, $source, $anno, $symbol,$ratio,$num);
GetOptions(
        "h|?"           =>\&USAGE,
	"ceRNA:s"	=>\$ceRNA,
	"diff:s"     	=>\$diff,
	"ratio:s"	=>\$ratio,
	"num:s"		=>\$num,
	"cfg:s"		=>\$cfg,
	"anno:s"	=>\$anno,
	"trans:s"	=>\$trans,
	"cis:s"		=>\$cis,
	"source:s"	=>\$source,
	"symbol:s"	=>\$symbol,
	"od:s"		=>\$od,
)or &USAGE;
&USAGE unless($ceRNA and $diff);

`mkdir -p $od`	unless(-d $od);
$od=abs_path($od);

my %config=&readConfig($cfg);
my $ratio=exists $config{keyratio}?$config{keyratio}:0.05;
my $min=exists $config{keymin}?$config{keymin}:100;

my %target=();
&readTarget($trans)     if(defined $trans);
&readTarget($cis)       if(defined $cis);
&readTarget($source)    if(defined $source);

my $vs = (split /\./,basename $diff)[1];
open(DIFF,$diff)||die $!;
my %diffs=();
while(<DIFF>){
	chomp;next if($_=~/^#|^$|^\s+/);
	my @tmp=split;
	my ($id,$regulate)=($tmp[0],$tmp[-1]);
	$diffs{$id}=$regulate;
}
close(DIFF);
my %pairs=();
open(CE, $ceRNA)||die $!;
open(NET,">$od/sub_${vs}_network.xls")||die $!;
open(NODE,">$od/sub_${vs}.Attribution.list")||die $!;
open(INT,">$od/sub_${vs}.Interaction.list")||die $!;
#print NET "ceRNA1\tceRNA1_regulated\tceRNA2\tceRNA2_regulated\n";
while(<CE>){
	chomp;next if($_=~/^$|^\s+/);
	if(/^#/){
		my($tit1,$tit2,$tit3)=split/\s+/,$_,3;
		print NET "ceRNA1\tceRNA1_regulated\tceRNA2\tceRNA2_regulated\t$tit3\n";
		next;
	}
	my ($rna1,$rna2,$info)=split /\t/,$_,3;
	my ($t1,$id1)=split /\;|:/,$rna1,2;
	my ($t2,$id2)=split /\;|:/,$rna2,2;
	if(exists $diffs{$id1} && exists $diffs{$id2}){
		print NET "$t1:$id1\t$diffs{$id1}\t$t2:$id2\t$diffs{$id2}\t$info\n";
		print INT "$id1\t$id2\t$t1-$t2\n";
		print NODE "$id1\t$t1\t$diffs{$id1}\n$id2\t$t2\t$diffs{$id2}\n";
	}elsif(exists $diffs{$id1} && !exists $diffs{$id2}){
		print NET "$t1:$id1\t$diffs{$id1}\t$t2:$id2\t--\t$info\n";
		print INT "$id1\t$id2\t$t1-$t2\n";
		print NODE "$id1\t$t1\t$diffs{$id1}\n$id2\t$t2\tnormal\n";
	}elsif(!exists $diffs{$id1} && exists $diff{$id2}){
		print NET "$t1:$id1\t--\t$t2:$id2\t$diffs{$id2}\t$info\n";
		print INT "$id1\t$id2\t$t1-$t2\n";
		print NODE "$id1\t$t1\tnormal\n$id2\t$t2\t$diffs{$id2}\n";
	}else{
		next;
	}
}
	
close(CE);
close(NET);
close(NODE);
close(INT);
my $cmd="Rscript $Bin/bin/random_walk.r --network $od/sub_${vs}_network.xls --outpath $od --from 1 --to 3 --dir FALSE --key $vs";
&run_or_die($cmd);
my $node_num=`wc -l $od/${vs}_score.txt`;
chomp($node_num);
$node_num=(split(/\s+/,$node_num))[0];

my $cut=int($ratio*$node_num)>$min?int($ratio*$node_num):$min;
`head -n $cut $od/${vs}_score.txt >$od/${vs}_Key_RNA.txt`;
open(KEY,"$od/${vs}_Key_RNA.txt")||die $!;
open(OUT,">$od/${vs}_Key_RNA_gene.xls")||die $!;
print OUT "#Gene\ttype\tKey_id\n";
while(<KEY>){
	chomp;
	my @tmp=split;
	my ($type,$id)=split(/\;|:/,$tmp[0],2);
	if($type eq "gene"){
		print OUT "$id\tgene\t$id\n";
	}else{
		if(exists $target{$id}){
			foreach my $g(keys %{$target{$id}}){
				print OUT "$g\t$type\t$id\n";
			}
		}else {
			print OUT "--\t$type\t$id\n";
		}
	}
}
close(OUT);
close(KEY);
#/share/nas2/genome/biosoft/R/3.3.2/bin/Rscript
if(defined $anno){
	my $kegg=(-e (glob("$anno/*.Kegg.pathway"))[0])?((glob("$anno/*.Kegg.pathway"))[0]):(glob("$anno/Result/*.Kegg.pathway"))[0];
	my $ko=(-e (glob("$anno/*Kegg.ko"))[0])?((glob("$anno/*Kegg.ko"))[0]):((glob("$anno/Result/*Kegg.ko"))[0]);
	&KEGGinfo($kegg,"$od/KEGG.info");
	&run_or_die("export PATH=/share/nas2/genome/biosoft/R/3.3.2/bin/:\$PATH && export LD_LIBRARY_PATH=/share/nas2/genome/biosoft/gcc/5.4.0/lib64:/share/nas2/genome/biosoft/gcc/5.4.0/lib:/share/nas2/genome/biosoft/pcre/lib:/share/nas2/genome/biosoft/Anaconda2/current/pkgs/xz-5.2.2-0/lib/:$LD_LIBRARY_PATH && Rscript $Bin/../../DEG_Analysis/enrich/enrich.r --deg $od/${vs}_Key_RNA_gene.xls --kegg $od/KEGG.info --prefix $vs  --od $od");
	my $num=exists $config{random_num}?$config{random_num}:5;
	if(-e "$od/${vs}_KEGG_pathway_enrich.list"){
		my $paths=&getKO("$od/${vs}_KEGG_pathway_enrich.list",$num);
		$cmd ="perl $Bin/KEGG/KEGG_network_construct.pl -gene $ko -od $od -ko $paths -key $od/${vs}_Key_RNA_gene.xls -prefix $vs";
		$cmd.=" -symbol $symbol"	if(defined $symbol);
		&run_or_die($cmd);
	}
}


################SUB FUNCTION
sub KEGGinfo{
	my($i,$o)=@_;
	open(IN,$i)||die $!;
	open(OUT,">$o")||die $!;
	while(<IN>){
		chomp;
		next if($_=~/^#/);
		my($path,$ko,$num,$gene,$K)=split(/\t/,$_); 
		my @genes=split(/\;/,$gene);
		foreach my $g(@genes){
			next if($g=~/^$/);
			print OUT "$ko\t$path\t$g\n";
		}
	}
	close(OUT);
	close(IN);
}
sub getKO{
	my ($file,$num)=@_;
	my $count=0;
	my @ko=();
	print "$file\n";
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

sub readTarget{
	my $file=shift;
	open(FILE,$file)||die $!;
	while(<FILE>){
		chomp;next if($_=~/^#/);
		my @tmp=split(/\s+/,$_);
		my @tars=split(/\;|,/,$tmp[1]);
		foreach (@tars){
			$target{$tmp[0]}{$_}++;
		}
	}
	close(FILE);
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
}
sub readConfig{
	my $cfg=shift;
	my %cfg ;
	open(IN,"$cfg") or die $!;
	while(<IN>){
		chomp;
		next if(/^#|^\s$|^$/);
		next if(/^Sample/);
		my @tmp = split /\s+/,$_;
		if($tmp[0] eq "Diff"){
#		if($tmp[0] eq "Com"){
#			my @sample1 = split /,/,$tmp[1];
#			my $g = join("_vs_",@sample1);
#			push @group,$g;
#		}elsif($tmp[0] eq "Sep"){
#			my ($s1,$s2) = split /;/,$tmp[1];
#			my @sample2 = split /,/,$s1;
#			my @sample3 = split /,/,$s2;
#			my $g1 = join("_",@sample2);
#			my $g2 = join("_",@sample3);
#			my $g3 = join("_vs_",$g1,$g2);
			push @group,$g3;
		}else{
			$cfg{$tmp[0]}=$tmp[1] if($tmp[0] eq "random_ratio" || $tmp[0] eq "random_min");
		}
	}
	close(IN);
	return %cfg;
}

sub USAGE{
	my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
	Options:

	-ceRNA	<file>	ceRNA_final file
	-diff	<file>	merged_diff file contained all diff RNA(lncRNA,circRNA,GENE)
	-cfg	<file>	cfg file	
	-anno	<dir>	anno dir
	-trans	<file>	lncRNA trans target gene file
	-cis	<file>	lncRNA cis target gene file
	-source	<file>	circRNA source gene file
	-symbol	<file>	id_name.list file 
	-od	<dir>	output dir
	
	-h	Help

Example:

USAGE
	print $usage;
	exit;
}


