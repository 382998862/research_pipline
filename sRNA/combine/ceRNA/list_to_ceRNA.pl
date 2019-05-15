use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path getcwd);
my ($mi2l,$mi2c,$mi2g,$od,$num,$pvalue,$fdr,$coexp,$diff,$queue);
GetOptions(
        "h|?"           =>\&USAGE,
	"mi2l:s"	=>\$mi2l,
	"mi2c:s"	=>\$mi2c,
	"mi2g:s"	=>\$mi2g,
	"od:s"     	=>\$od,
	"diff:s"	=>\$diff,
	"num:s"		=>\$num,
	"pvalue:s"	=>\$pvalue,
	"fdr:s"		=>\$fdr,
	"coexp:s"	=>\$coexp,
	"queue:s"	=>\$queue,
)or &USAGE;
&USAGE unless (($mi2l || $mi2c || $mi2g) && $od);

`mkdir -p $od`	unless(-d $od);
$od = abs_path($od);
$mi2l = abs_path($mi2l) if(defined $mi2l);
$mi2c = abs_path($mi2c) if(defined $mi2c);
$mi2g = abs_path($mi2g) if(defined $mi2g);
my $coexp_dir;
if(defined $coexp){
	$coexp = abs_path($coexp);
	$coexp_dir = dirname $coexp;
}
print "$coexp_dir\n";
$num||=5;
$pvalue||=0.01;
$fdr||=0.05;
$queue||="medical.q";

my %totalRNA=();
my %miRNA=();
my %targets=();

my %type_hash=();
$type_hash{"lncRNA"}=$mi2l if(defined $mi2l);
$type_hash{"circRNA"}=$mi2c if(defined $mi2c);
$type_hash{"gene"}=$mi2g if(defined $mi2g);

foreach my $type(sort keys %type_hash){
	print "$type\t$type_hash{$type}\n";
	open(IN,"$type_hash{$type}") or die $!;
	while(<IN>){
		chomp;
		my @info = split;
		my @rnas = split /,|;|\//,$info[1];
		$miRNA{$info[0]}++;
		foreach my $r(@rnas){
			push @{"$type:$r"},$info[0];
			push @{$targets{$r}},$info[0];
			$totalRNA{$type}{$r}++;
		}
	}
	close(IN);
}

my @RNAs= sort keys %totalRNA;
print "$RNAs[0]\t$RNAs[1]\t$RNAs[2]\n";
if(!-e "$od/ceRNA_pair_adjust_p.txt"){
	open(OUT,">$od/tmp_ceRNA_pairs.txt")||die $!;
	for(my $i=0;$i<@RNAs-1;$i++){
		for(my $j=$i+1;$j<@RNAs;$j++){
			my @RNA1=keys %{$totalRNA{$RNAs[$i]}};
			my @RNA2=keys %{$totalRNA{$RNAs[$j]}};
			foreach my $r1(@RNA1){
				my @cerna1=&unique(\@{"$RNAs[$i]:$r1"});
				my %hash_a = map{$_=>1} @cerna1;
				foreach my $r2(@RNA2){
					my @cerna2=&unique(\@{"$RNAs[$j]:$r2"});
					my @overlap=grep {$hash_a{$_}} @cerna2;
					if(@overlap>0){
						print OUT "$RNAs[$i]:$r1\t$RNAs[$j]:$r2\t".scalar(@cerna1)."\t".scalar(@cerna2)."\t".scalar(@overlap)."\n";
					}
				}
			}
		}
	}
	close(OUT);
}

`mkdir -p $od/tmp` if(!-e "$od/tmp");
chdir ("$od/tmp");
#my $cmd= "awk -F \$\"\\t\" '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5}' $od/tmp_ceRNA_pairs.txt >$od/tmp/tmp.txt && split -l 100000 $od/tmp/tmp.txt -d ceRNA -a 4\n";
my $cmd = "split -l 100000 $od/tmp_ceRNA_pairs.txt -d ceRNA -a 4\n";
&run_or_die($cmd);


my @files=glob("$od/tmp/ceRNA*");
my $totalnum=scalar(keys %miRNA);
print "Total miRNA num is $totalnum\n";

open(SH,">$od/pvalue.sh")||die $!;
my $command="cat ";
for(my $i=0;$i<@files;$i++){
        my $id=sprintf("%04d",$i);
        print SH "Rscript $Bin/hyper_pvalue.r -i $od/tmp/ceRNA$id -o $od/tmp/p.ceRNA$id.txt -t $totalnum\n";
        $command .="$od/tmp/p.ceRNA$id.txt ";
}
$command .=">$od/tmp/pvalue.txt";
close(SH);

&qsub("$od/pvalue.sh");

&run_or_die("cat $od/tmp/p.ceRNA*.txt >$od/tmp/pvalue.txt ");
&run_or_die("Rscript $Bin/p_adjust.r -i $od/tmp/pvalue.txt -o $od/tmp/adjust_pvalue.txt");
&run_or_die("paste $od/tmp_ceRNA_pairs.txt $od/tmp/adjust_pvalue.txt > $od/ceRNA_pair_adjust_p.txt && sed -i '1i\#ceRNA1\\tceRNA2\\tceRNA1_num\\tceRNA2_num\\toverlap_num\\tpvalue\\tFDR' $od/ceRNA_pair_adjust_p.txt");

############filter the result based on the pvalue and fdr 
my ($ceRNA_result,$diff_ceRNA,$diff_node,%hash_coexp);
if(defined $coexp){
	&show_log("Start Check if the coexpression result exists or not...");
	$ceRNA_result = "$od/coexp_ceRNA_pair_adjust_p_Sig.txt";
	$diff_ceRNA = "$od/coexp_ceRNA_pair_adjust_p_Sig_diff.txt";
	$diff_node = "$od/coexp_ceRNA_pair_adjust_p_Sig_diff.node";
	my $times =10;
	for(my $t=1;$t<$times;$t++){
		if(!-e "$coexp_dir/coexpression_finish"){
			for(my $i=1;$i<$times;$i++){
				sleep(60);
				$i++;
				if(-e "$coexp_dir/coexpression_finish"){
					&show_log("I have waited for 60*$i seconds.You finished finally.Cheer up!");last;
				}
			}
		}
		last if(-e "$coexp_dir/coexpression_finish");
	}
	if(!-e "$coexp_dir/coexpression_finish"){
		&show_log("I have waited for 1h,but you coexpression result is still null.I will be killed,please rerun this program after your coexpression!");
		die;
	}else{
		&show_log("coexpression analysis done!this program will generates coexp_ceRNA file.");
	}

	open(CO,"$coexp")or die $!;
	while(<CO>){
		chomp;
		next if(/^#|^RNA1/);
		my @tmp = split /\t/,$_;
		$hash_coexp{$tmp[0]}{$tmp[1]}=join("\t",$tmp[2],$tmp[3]);
		$hash_coexp{$tmp[1]}{$tmp[0]}=join("\t",$tmp[2],$tmp[3]);
	}
	close(CO);
	&show_log("read coexpression finish!");
}else{
	&show_log("you have no coexpression analysis!this pregram will generates ceRNA file.");
	$ceRNA_result = "$od/ceRNA_pair_adjust_p_Sig.txt";
	$diff_ceRNA = "$od/ceRNA_pair_adjust_p_Sig_diff.txt";
	$diff_node = "$od/ceRNA_pair_adjust_p_Sig_diff.node";
}


open(CERNA,"$od/ceRNA_pair_adjust_p.txt")||die $!;
open(SIG,">$ceRNA_result")||die $!;
my $head="\#ceRNA1\tceRNA2\tceRNA1_num\tceRNA2_num\toverlap_num\tceRNA1_miRNA_uniq\tceRNA2_miRNA_uniq\toverlap_miRNA\tpvalue\tFDR";
$head .= "\tcoexp_cor\tcoexp_pvalue" if(defined $coexp);
print SIG "$head\n";
while(<CERNA>){
	chomp;
	next if(/^#|^ceRNA/);
	my ($c1,$c2,$cn1,$cn2,$over,$p,$f)=split /\s+/,$_ ;
	next if($over<$num || $p>$pvalue ||$f>$fdr);
#	print "$_\n";
	my ($type1,$id1)=split /:/,$c1,2;
	my ($type2,$id2)=split /:/,$c2,2;
	my @mi1 = @{$targets{$id1}};
	my @mi2 = @{$targets{$id2}};
	my %hash_a = map{$_=>1} @mi1;
	my %hash_b = map{$_=>1} @mi2;
	my @uniq1 = grep {!$hash_b{$_}} @mi1;
	my @uniq2 = grep {!$hash_a{$_}} @mi2;
	my @miover = grep {$hash_b{$_}} @mi1;
	if(-e $coexp){
		if(exists $hash_coexp{$id1}{$id2}){
			print SIG "$c1\t$c2\t$cn1\t$cn2\t$over\t",join(",",@uniq1),"\t",join(",",@uniq2),"\t",join(",",@miover),"\t$p\t$f\t$hash_coexp{$id1}{$id2}\n";
		}
	}else{
		print SIG "$c1\t$c2\t$cn1\t$cn2\t$over\t",join(",",@uniq1),"\t",join(",",@uniq2),"\t",join(",",@miover),"\t$p\t$f\n";
	}
}
close(CERNA);
close(SIG);

############
if(defined $diff){
        my %diffList=();
        open(DIFF,$diff)||die $!;
        while(<DIFF>){
                chomp;next if($_=~/^#/);
		my@tmp=split(/\t/,$_);
                $diffList{$tmp[0]}++;
        }
        close(DIFF);
        open(CE,"$ceRNA_result")||die $!;
        open(DIFFCE,">$diff_ceRNA")||die $!;
	
	my $head1 = "#ceRNA1\tceRNA2\tceRNA1_num\tceRNA2_num\toverlap_num\tceRNA1_miRNA_uniq\tceRNA2_miRNA_uniq\toverlap_miRNA\tpvalue\tFDR";
	$head1 .= "\tcoexp_cor\tcoexp_pvalue" if(defined $coexp);
	print DIFFCE "$head1\n";
        my %types=();
        while(<CE>){
                chomp;
                my @tmp=split(/\t/,$_);
                my ($t1,$id1)=split(/:/,shift @tmp,2);
                my ($t2,$id2)=split(/:/,shift @tmp,2);
                if(exists $diffList{$id1} || exists $diffList{$id2}){
                        print DIFFCE "$id1\t$id2\t",join("\t",@tmp),"\n";
                        $types{$id1}=$t1;
                        $types{$id2}=$t2;
                }
        }

        open(NODE,">$diff_node")||die $!;
        foreach my $rna(keys %types){
                print NODE "$rna\t$types{$rna}\n";
        }
        close(NODE);
        close(DIFFCE);
        close(CE);
}

################SUB FUNCTION
sub qsub()
{
        my $shfile= shift;
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $shfile --queue $queue";
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


sub readConfig{
	my $configFile=shift;
	my $d=Config::General->new(-ConfigFile => "$configFile");
	my %config=$d->getall;	
	return %config;
}
sub intersect{
        my ($s1,$s2)=@_;
        my @a=@{$s1};
        my @b=@{$s2};
        my %a=map{$_ => 1} @a;
        my %b=map{$_ => 1} @b;
        my @inter=grep{$a{$_}} @b;
        return @inter;
}
sub unique{
        my $s=shift;
        my @set=@{$s};
        my %hash=();
        foreach (@set){ $hash{$_}++;}
        my @uniq=sort keys %hash;
        return @uniq;
}

sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE{
	my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
	Options:
	-mi2l	<file>		miRNA target lncRNA file, #miRNA_id lncRNA1,lncRNA2,lncRNA3...
	-mi2c	<file>		miRNA target circRNA file, #miRNA_id circRNA1,circRNA2,circRNA3...
	-mi2g	<file>		miRNA target mRNA file, #miRNA_id mRNA1,mRNA2,mRNA3...
	-od	<path>		output path, forced
	-diff	<file>		Diff file, if exists this file will give additional output files, not must. 
	-coexp	<file>		coexp file, if exists, get coexp_ceRNA result,else ceRNA result
				if your coexp result is not exists,this program will wait for 10*10*60s until "coexp_dir/coexpression_finish" exists,
				or the program die and you have to rerun this script.
	-num	<int>		the min number of common miRNA for ceRNA pairs, default 5
	-pvalue	<float>		the min pvalue, default 0.01
	-fdr	<float>		the min fdr, default 0.05
	
	-h	Help

Example: perl $0 -mi2l miRNA2lncRNA.xls -mi2c miRNA2circRNA.xls -mi2g miRNA2mRNA.xls -od ceRNA_dir -coexp coexp_result

USAGE
	print $usage;
	exit;
}


