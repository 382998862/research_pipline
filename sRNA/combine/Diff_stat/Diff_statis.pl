use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
use Cwd qw(abs_path getcwd);
my ($DEG,$od,$detail_cfg,$data_cfg,$lnc_gff,$gene_gff,$mi_loc,$ideogram);
GetOptions(
        "h|?"           =>\&USAGE,
	"lnc_gff:s"	=>\$lnc_gff,
	"gene_gff:s"	=>\$gene_gff,
	"mi_loc:s"	=>\$mi_loc,
	"DEG:s"		=>\$DEG,
	"cfg1:s"	=>\$data_cfg,
	"cfg2:s"	=>\$detail_cfg,
	"od:s"     	=>\$od,
	"ideogram:s"	=>\$ideogram,
)or &USAGE;
&USAGE unless ($DEG and $od and $detail_cfg and $data_cfg);

`mkdir -p $od`	unless(-d $od);
$od=abs_path($od);
if(!defined $ideogram){
	if(-e "$od/../ideogram.info"){
		$ideogram="$od/../ideogram.info"
	}else{
		print "You should give ideogram parameter!\n";
		die;
	}
}
########################################
my %sample=&relate($data_cfg);
my @diff=();
open(CFG,$detail_cfg)||die $!;
while(<CFG>){
	chomp;my @tmp=split(/\s+/,$_);
	if($tmp[0] eq "Com"){
		$tmp[1]=~s/,/_vs_/;
		push @diff,$tmp[1];
	}elsif($tmp[0] eq "Sep"){
		$tmp[1]=~s/,/_/g;
		$tmp[1]=~s/\;/_vs_/;
		push @diff,$tmp[1];
	}
}
close(CFG);

my @RNAs=("sRNA","mRNA","lncRNA","circRNA");
my %locs=();
my %total=();	##total diff in all diff group
&readGFF($lnc_gff)	if(defined $lnc_gff);		
&readGFF($gene_gff)	if(defined $gene_gff);
&miLoc($mi_loc)		if(defined $mi_loc);

my $Rscript = "/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript";
open(TP,">$od/diff_type_num.xls")||die $!;
print TP "group\tRNA_type\tNumber\n";
my %vs=();      #$vs{T01_vs_T02}{mRNA}
foreach my $v(@diff){####plot each diff group
	&trans($v);
	`mkdir $od/$v`	unless(-d "$od/$v");
	my $circos_cmd="perl $Bin/../tools/circos_v1.pl --chr $ideogram --od $od/$v ";
	open(PIE,">$od/$v/$v.stat.xls")||die $!;
	if(exists $vs{$v}{mRNA}){
		my $deg="$DEG/gene/$vs{$v}{mRNA}/$vs{$v}{mRNA}.DEG_final.xls";
		my $tmp2=`wc -l $deg`;
		my $l2=(split(/\t/,$tmp2))[0]-1;
		print PIE "gene\t$l2\n";
		&diff_loc($deg,"$od/$v/gene.txt","gene");
		$circos_cmd .="--circle $od/$v/gene.txt --type histogram ";
		print TP "$v\tgene\t$l2\n";
	}
        if(exists $vs{$v}{lncRNA}){
		my $deg="$DEG/lncRNA/$vs{$v}{lncRNA}/$vs{$v}{lncRNA}.DEG_final.xls";
                my $tmp1=`wc -l $deg`;
                my $l1=(split(/\t/,$tmp1))[0]-1;
                print PIE "lncRNA\t$l1\n";
                &diff_loc($deg,"$od/$v/lncRNA.txt","lncRNA");
                $circos_cmd .="--circle $od/$v/lncRNA.txt --type histogram ";
		print TP "$v\tlncRNA\t$l1\n";
        }
	if(exists $vs{$v}{circRNA}){
		my $deg="$DEG/circRNA/$vs{$v}{circRNA}/$vs{$v}{circRNA}.DEG_final.xls";
		my $tmp=`wc -l $deg`;
		my $l=(split(/\t/,$tmp))[0]-1;
		print PIE "circRNA\t$l\n";
		&diff_loc($deg,"$od/$v/circRNA.txt","circRNA");
		$circos_cmd .="--circle $od/$v/circRNA.txt --type histogram ";
		print TP "$v\tcircRNA\t$l\n";
	}
	if(exists $vs{$v}{sRNA}){
		my $deg="$DEG/sRNA/$vs{$v}{sRNA}/$vs{$v}{sRNA}.DEG_final.xls";
                my $tmp=`wc -l $deg`;
                my $l=(split(/\t/,$tmp))[0]-1;
                print PIE "sRNA\t$l\n";
                print TP "$v\tsRNA\t$l\n";
		&diff_loc($deg,"$od/$v/miRNA.txt","sRNA");
		$circos_cmd .="--circle $od/$v/miRNA.txt --type histogram ";
        }
	close(PIE);
	&run_or_die("$Rscript $Bin/plot/pie_chart_test.r -i $od/$v/$v.stat.xls -s $v -o $od/$v/$v.pie.png");
	&run_or_die("$circos_cmd");
}
close(TP);
=head
################### Total difference group plot
`mkdir $od/Total`       unless(-d "$od/Total");
`mv $od/diff_type_num.xls $od/Total`;
open(DIFF,">$od/Total/Diff_RNA.list")||die $!;
foreach my $RNA(@RNAs){
	open(OUT,">$od/Total/$RNA.txt")||die $!;
	my @degs=keys %{$total{$RNA}};
	foreach my $deg(@degs){
		my @types=&unique(\@{$deg});
		my $col;
		my $final=$types[0];
		if(scalar(@types)==2){
			$col="green";
			$final="both";
		}elsif($types[0] eq "up"){
			$col="red";
		}elsif($types[0] eq "down"){
			$col="blue";
		}
		print DIFF "$RNA\t$deg\t$final\n";
		next if($RNA eq "miRNA");

		if($RNA eq "circRNA"){
			my $deg_tmp=$deg;
			$deg_tmp=~s/\|/:/;
			my @tmp=split(/:/,$deg_tmp);
			print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t1\tcolor=$col\n";
		}else{
			print OUT "$locs{$deg}{chr}\t$locs{$deg}{tss}\t$locs{$deg}{tts}\t1\tcolor=$col\n";
		}
	}
	close(OUT);
}
close(DIFF);
if(@diff>1){
	my $command="perl $Bin/../tools/circos_v1.pl --chr $in/ideogram.info --od $od/Total ";
	$command .="--circle $od/Total/mRNA.txt --type heatmap "        if(-e "$od/Total/mRNA.txt");
	$command .="--circle $od/Total/lncRNA.txt --type heatmap "        if(-e "$od/Total/lncRNA.txt");
	$command .="--circle $od/Total/circRNA.txt --type heatmap "        if(-e "$od/Total/circRNA.txt");
	$command .="--circle $od/Total/miRNA.txt --type heatmap "        if(-e "$od/Total/miRNA.txt");
	&run_or_die($command);
	&run_or_die("$Rscript $Bin/plot/geom_bar.r --input $od/Total/diff_type_num.xls --outpath $od/Total --key total_diff");
}
=cut
###############################
#
#	Self defined function
###############################
sub relate{     ##获取不同RNA样品的对应关系
        my $file=shift;
        open(REL,$file)||die $!;
        my %sample=();
        while(<REL>){
                chomp;
                next if($_ !~/^SampleID/);
                my @tmp=split(/\s+/,$_);
                $sample{$tmp[-1]}{sRNA}=$tmp[2];
                $sample{$tmp[-1]}{lncRNA}=$tmp[3];
                $sample{$tmp[-1]}{circRNA}=$tmp[3];
                $sample{$tmp[-1]}{mRNA}=$tmp[3];
        }
        close(REL);
        return %sample;
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

sub trans{  ###transfer diff group for different RNA type
        my $v=shift;
        my ($g1,$g2)=split(/_vs_/,$v,2);
	my @group1=split(/_/,$g1);
        my @group2=split(/_/,$g2);
        foreach my $RNA(@RNAs){
                my @ng1=();
		my @ng2=();
                foreach my $g1(@group1){
                        push @ng1,$sample{$g1}{$RNA};
                }
                foreach my $g2(@group2){
                        push @ng2,$sample{$g2}{$RNA};
                }
		my @nng1=&unique(\@ng1);
		my @nng2=&unique(\@ng2);
                my $ngs1=join("_",@nng1);
                my $ngs2=join("_",@nng2);
                $vs{$v}{$RNA}="$ngs1"."_vs_"."$ngs2";
		print "$v:$RNA:$vs{$v}{$RNA}\n";
        }
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
sub mean{
	my $s=shift;
	my @set=@{$s};
	my $sum=0;
	foreach my $s(@set){$sum=$sum+$s;}
	return $sum/scalar(@set);
}
sub readGFF{####get lncRNA/mRNA locs based on GFF files
	my $gff=shift;
	open(GFF,$gff)||die $!;
	while(<GFF>){
		chomp;
		my($chr,$source,$type,$tss,$tts,$tmp1,$strand,$tmp2,$info)=split(/\t/,$_);
#		next if($type ne "mRNA" && $type ne "gene" && $type ne "lincRNA_gene");
		if($type eq "gene" ){
			$info=~/ID=(.*?)$/;
			my $id=(split(/;/,$1))[0];
			$locs{$id}{chr}=$chr;
			$locs{$id}{tss}=$tss;
			$locs{$id}{tts}=$tts;
		}
		if($type eq "lincRNA_gene"){
			if($info=~/ID=gene:(.*);Name=(.*)/){
				$locs{$1}{chr}=$chr;
				$locs{$1}{tss}=$tss;
				$locs{$1}{tts}=$tts;
			}
		}
		if($type eq "mRNA"){
			$info=~/ID=(.*?)$/;
			my $id=(split(/;/,$1))[0];
                        $locs{$id}{chr}=$chr;
                        $locs{$id}{tss}=$tss;
                        $locs{$id}{tts}=$tts;
		}
		if($type eq "transcript" || $type eq "lincRNA"){
			if($info=~/ID=transcript:(.*);Parent=(.*)/){
				$locs{$1}{chr}=$chr;
				$locs{$1}{tss}=$tss;
				$locs{$1}{tts}=$tts;
			}
		}
	}
	close(GFF);		

}
sub getExp{
        my ($file,$rna)=@_;
	return if(!-e $file);
        open(FILE,$file)||die $!;
        my $header=<FILE>;
        chomp($header);my @heads=split(/\t/,$header);shift @heads;
        while(<FILE>){
                chomp;
                my @tmp=split(/\t/,$_);
                for(my $i=1;$i<@tmp;$i++){
                        ${$rna}{$heads[$i-1]}{$tmp[0]}=$tmp[$i];
                }
        }
        close(FILE);
}

sub diff_loc{	#get circos format based on deg.files
	my ($in,$out,$type)=@_;  #$in is deg.file  $out deg.loc $type mRNA/lncRNA/circRNA
	open(IN,$in)||die $!;
	open(OUT,">$out")||die $!;
	my $head=<IN>;
	chomp($head);my @tmp=split(/\t/,$head);
	my %hash=();
	for(my $i=0;$i<@tmp;$i++){$hash{$tmp[$i]}=$i;}	
	while(<IN>){
		chomp;
		next if($_=~/^#/);
		my @tmp=split(/\t/,$_);
		my ($col,$logFC);
		if($tmp[-1] eq "up"){$col="red";}elsif($tmp[-1] eq "down"){$col="blue";}
		my ($logp,$order);
		if(exists $hash{PValue}){$order=$hash{PValue};}elsif(exists $hash{FDR}){$order=$hash{FDR};}
		$logp=&log10($tmp[$order]+0.0001);
		$logp=-$logp;
		###not all DEG.xls would return pvalue, so taken 1-FDR as significant
		if($type eq "circRNA"){
			my $circID=$tmp[0];
			$circID=~s/\|/:/;
			my ($id,$start,$end)=split(/:/,$circID);
			##chr start end -log10pvalue color=red
			print OUT "chr$id\t$start\t$end\t",sprintf("%.1f",$logp),"\tcolor=$col\n";
		}else{
			print OUT "chr$locs{$tmp[0]}{chr}\t$locs{$tmp[0]}{tss}\t$locs{$tmp[0]}{tts}\t",sprintf("%.1f",$logp),"\tcolor=$col\n";
		}
		$total{$type}{$tmp[0]}++;##整合所有差异分组的差异情况
		push @{$tmp[0]},$tmp[-1];
	}
	close(IN);
	close(OUT);
}
sub diff_miRNA{##get miRNA diff info
	my $in=shift;
	open(IN,$in)||die $!;
	while(<IN>){
		chomp;
		next if($_=~/^#/);
		my @tmp=split(/\t/,$_);
		$total{miRNA}{$tmp[0]}++;
		push @{$tmp[0]},$tmp[-1];
	}
	
	close(IN);
}



sub log10 {
	my $n = shift;
	return log($n)/log(10);
}
sub which{
	my ($s,$value)=@_;
	my @set=@{$s};
	my $loc=-1;
	for(my $i=0;$i<@set;$i++){
		if($set[$i] eq $value){
			$loc=$i;
		}
	}
	return $loc;
}
sub getDeg{
        my $deg=shift;
        my %hash=();
        open(DEG,$deg)||die $!;
        while(<DEG>){
                chomp;
                next if($_=~/^#/);
                my @tmp=split;
                $hash{$tmp[0]}=$tmp[-1];
        }
        close(DEG);
        return %hash;
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
	-od		<path>	output path, forced
	-cfg1		<file>	data config file
	-cfg2		<file>	detail config
	-lnc_gff	<file>	lncRNA gff
	-gene_gff	<file>	gene gff
	-mi_loc		<file>	miRNA loc
	-DEG		<path>	DEG Analysis path
        -ideogram    	<file>	ideogram file for circos plot

	-h	Help

Example:

USAGE
	print $usage;
	exit;
}


