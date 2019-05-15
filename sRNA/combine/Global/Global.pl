use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use newPerlBase;
use Cwd qw(abs_path getcwd);
my ($lnc,$lnc_gff,$circ,$mi,$mi_loc,$gene,$gene_gff,$relate,$od,$ideogram);
GetOptions(
        "h|?"           =>\&USAGE,
	"lnc:s"		=>\$lnc,
	"circ:s"	=>\$circ,
	"mi:s"		=>\$mi,
	"gene:s"	=>\$gene,
	"ideogram:s"	=>\$ideogram,
	"lnc_gff:s"	=>\$lnc_gff,
	"gene_gff:s"	=>\$gene_gff,
	"mi_loc:s"	=>\$mi_loc,
	"od:s"     	=>\$od,
	"relate:s"	=>\$relate,
)or &USAGE;
&USAGE unless ($od and $relate and $ideogram);

`mkdir -p $od`	unless(-d $od);
$od=abs_path($od);

my %sample=();  #$sample{T01}{mRNA}
my %RNAs=();
my @samples=&relate($relate);

my %locs=();
my %total=();##total diff in all diff group
&readGFF($lnc_gff)	if(defined $lnc_gff);		
&readGFF($gene_gff)	if(defined $gene_gff);

my $Rscript="/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript";
#############################pie
`mkdir $od/pie`	unless(-d "$od/pie");
open(PIE,">$od/pie/Total.stat.xls")||die $!;
if(defined $gene){
	my $line=`wc -l $gene`;
	chomp($line);print PIE "mRNA\t",(split(/\s+/,$line))[0]-1,"\n";
	&getExp($gene,"mRNA");
}
if(defined $lnc){
        my $line=`wc -l $lnc`;
        chomp($line);print PIE "lncRNA\t",(split(/\s+/,$line))[0]-1,"\n";
	&getExp($lnc,"lncRNA");
}
if(defined $mi){
        my $line=`wc -l $mi`;
        chomp($line);print PIE "miRNA\t",(split(/\s+/,$line))[0]-1,"\n";
	&getExp($mi,"miRNA");
	&miLoc($mi_loc);
}

if(defined $circ){
        my $line=`wc -l $circ`;
        chomp($line);print PIE "circRNA\t",(split(/\s+/,$line))[0]-1,"\n";
	&getExp($circ,"circRNA");
	&circLoc($circ);
}
close(PIE);
&run_or_die("$Rscript $Bin/pie_chart_test.r -i $od/pie/Total.stat.xls -s Total -o $od/pie/Total.pie.png");
############################
`mkdir $od/barplot`	unless(-d "$od/barplot");
`mkdir $od/circos`	unless(-d "$od/circos");
my @rnas=("mRNA","lncRNA","circRNA","miRNA");

my @chrs=();
open(IDEOGRAM,$ideogram)||die $!;
while(<IDEOGRAM>){
	chomp;next if($_=~/^#/);next if(/^band/);
	my @tmp=split(/\s+/,$_);
	push @chrs,$tmp[3];
}
close(IDEOGRAM);

foreach my $s(@samples){####plot each sample
	`mkdir $od/circos/$s`	unless(-d "$od/circos/$s");
	`mkdir $od/barplot/$s`   unless(-d "$od/barplot/$s");
	my $cmd="perl $Bin/../tools/circos_v1.pl -chr $ideogram -od $od/circos/$s ";
	open(BAR,">$od/barplot/$s/barplot.txt")||die $!;
	print BAR "chr\t",join("\t",@chrs),"\n";
	foreach my $rna(@rnas){#${$rna}{sample}{gene}
		next if(!exists $sample{$s}{$rna});
		my $id=$sample{$s}{$rna};
		print "$id\n";
		my @genes=keys %{${$rna}{$id}};
		my %chr_exp=();
		open(OUT,">$od/circos/$s/$rna.txt")||die $!;
		foreach my $g(@genes){
			my $e=${$rna}{$id}{$g};
			$chr_exp{$locs{$g}{chr}}=$chr_exp{$locs{$g}{chr}}+$e;
			$e=&log2($e+1);
			$e=int($e);
			next if($e==0);
			print OUT "chr$locs{$g}{chr}\t$locs{$g}{tss}\t$locs{$g}{tts}\t$e\n";
		}
		close(OUT);
		$cmd .=" --circle $od/circos/$s/$rna.txt  --type histogram ";
		my @exps=values %chr_exp;
		my $sum=&sum(\@exps);
		print BAR "$rna";
		foreach my $c(@chrs){
			if(exists $chr_exp{$c}){
				my $ratio=$chr_exp{$c}/$sum;	print BAR "\t$ratio";
			}else{print BAR "\t0";}
		}		
		print BAR "\n";
	}
	close(BAR);
	`awk -F \$"\\t" '{for(i=1;i<=NF;i++){a[FNR,i]=\$i}}END{for(i=1;i<=NF;i++){for(j=1;j<=FNR;j++){printf a[j,i]"\\t"}print ""}}' $od/barplot/$s/barplot.txt >$od/barplot/$s/barplot_t.txt`;
	`sed -i 's/\\t\$//g' $od/barplot/$s/barplot_t.txt`;
	&run_or_die("cd $od/barplot/$s && $Rscript $Bin/multi_RNA_exp_barplot.R $od/barplot/$s/barplot_t.txt $od/barplot/$s/$s.exp_barplot");
	&run_or_die($cmd);
}

###############################
#	Self defined function
###############################
sub relate{     ##获取不同RNA样品的对应关系
        my $file=shift;
        open(REL,$file)||die $!;
	#Sample mRNA    lncRNA  circRNA miRNA
        my @samples=();
        while(<REL>){
                chomp;
                next if($_ !~/^SampleID/);
                my @tmp=split(/\s+/,$_);
                $sample{$tmp[-1]}{miRNA}=$tmp[2]; 
		$sample{$tmp[-1]}{lncRNA}=$tmp[3];
		$sample{$tmp[-1]}{circRNA}=$tmp[3];
		$sample{$tmp[-1]}{mRNA}=$tmp[3];
                push @samples,$tmp[-1];
        }
        close(REL);
        return @samples;
}

sub unique{
	my $s=shift;
	my @set=@{$s};
	my @new=();
	my %hash=();
	for(my $i=0;$i<@set;$i++){
		$hash{$set[$i]}++;
		push @new,$set[$i]	if($hash{$set[$i]}==1);

	}
	return @new;
}
sub readGFF{####get lncRNA/mRNA locs based on GFF files
	my $gff=shift;
	open(GFF,$gff)||die $!;
	while(<GFF>){
		chomp;next if($_=~/^#/);
		my($chr,$source,$type,$tss,$tts,$tmp1,$strand,$tmp2,$info)=split(/\t/,$_);
		$info=~/ID=(.*?)$/;
		my $id=(split(/;/,$1))[0];
		$id=~s/transcript:|gene://;
		$locs{$id}{chr}=$chr;
		$locs{$id}{tss}=$tss;
		$locs{$id}{tts}=$tts;
	}
	close(GFF);		
}
sub miLoc{
	my $loc=shift;
	open(LOC,$loc)||die $!;
	while(<LOC>){
		chomp;next if($_=~/^#/);
		my @tmp=split(/\t/,$_);##precursor loc as miRNA Location
		$tmp[1]=~s/^chr//;
		$locs{$tmp[0]}{chr}=$tmp[1];
		$locs{$tmp[0]}{tss}=$tmp[-3];
		$locs{$tmp[0]}{tts}=$tmp[-2];
	}
	close(LOC);
}


sub circLoc{
	my $circ=shift;
	open(CIRC,$circ)||die $!;
	while(<CIRC>){
		next if($_=~/^#/);
		chomp;my @tmp=split(/\t/,$_);
		my $circID=$tmp[0];
                $circID=~s/\|/:/;
                my ($id,$start,$end)=split(/:/,$circID);
		$locs{$tmp[0]}{chr}=$id;
		$locs{$tmp[0]}{tss}=$start;
		$locs{$tmp[0]}{tts}=$end;
	}
	close(CIRC);
}
sub log2 {
	my $n = shift;
	return log($n)/log(2);
}
sub sum{
	my $s=shift;
	my @set=@{$s};
	my $sum=0;
	foreach my $s(@set){$sum=$sum+$s;}
	return $sum;
}
sub getExp{
	my ($file,$rna)=@_;
	open(FILE,$file)||die $!;
	my $header=<FILE>;
	chomp($header);my @heads=split(/\t/,$header);
	while(<FILE>){
		chomp;
		my @tmp=split(/\t/,$_);
		for(my $i=1;$i<@tmp;$i++){
			${$rna}{$heads[$i]}{$tmp[0]}=$tmp[$i];	#${mRNA}{L01}{GENE1}=35.9;
		}
	}
	close(FILE);
}
###############################################################
sub readConfig{
	my $configFile=shift;
	my $d=Config::General->new(-ConfigFile => "$configFile");
	my %config=$d->getall;	
	return %config;
}
sub qsub()
{
        my $shfile= shift;
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $shf
ile";
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
	-lnc		<file>	lncRNA fpkm.list
	-lnc_gff	<file>	lncRNA gff
	-circ		<file>	circRNA exp.list
	-mi		<file>	miRNA exp.list
	-mi_loc		<file>	miRNA_pos.list
	-gene		<file>	gene fpkm.list
	-gene_gff	<file>	gene gff
	-od		<path>	output path, forced
	-relate		<file>	sample relate info
	-ideogram	<file>	ideogram file
	-h	Help

Example:

USAGE
	print $usage;
	exit;
}


