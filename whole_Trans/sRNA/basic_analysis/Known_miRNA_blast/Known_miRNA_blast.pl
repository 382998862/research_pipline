#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.1";

my($miRBase_mature,$idir,$main_cfg,$reads,$key,$species,$od,$step,$sbs);
GetOptions(
	"help|?"	=>\&USAGE,
	"mature:s"	=>\$miRBase_mature,
	"idir:s"	=>\$idir,
	"maincfg:s"	=>\$main_cfg,
	"reads:s"	=>\$reads,
	"miR:s"		=>\$key,
	"species:s"	=>\$species,
	"od:s"		=>\$od,
	"step:s"	=>\$step,
	"onestep"	=>\$sbs,
) or &USAGE;
&USAGE unless ($miRBase_mature and $key and $od);
&USAGE unless (($idir and $main_cfg) or $reads);
&USAGE if (($idir and $main_cfg) and $reads);

$miRBase_mature = &ABSOLUTE_DIR($miRBase_mature);
$od = &MAKE_DIR($od);

$step ||= 0;
$species ||= "NONE";

################################################################################################################
#main code
################################################################################################################

#read main.cfg
my %Raw_FQ;
unless (defined $reads and -f $reads) {
	open(IN,$main_cfg) or die $!;
	while(<IN>){
		chomp;
		next if(/^#/||/^$/);
		if(/^NAME/){
			my $name=(split/\s+/,$_)[1];
			my $fq=<IN>;
			$Raw_FQ{$name}=(split/\s+/,$fq)[1];
			unless (-e $Raw_FQ{$name}) {
				print "Check Your config file:$main_cfg!Expecially $name\n\n";die;
			}
		}
	}
	close IN;
}
else {
	$reads = &ABSOLUTE_DIR($reads);
}


################## Get Known Mature miRNA
if($step==0){
        stepStart(0,"Get Known Mature miRNA");
		my $miRBase_dir = dirname($miRBase_mature);
		my $miRBase_pre = "$miRBase_dir/hairpin.fa";
		runOrDie("perl $Bin/bin/select_fa.pl -i $key -fa $miRBase_mature -o $od/$key.mature.fa >>$od/select_mature.out");
		runOrDie("perl $Bin/bin/select_fa.pl -i $key -fa $miRBase_pre -o $od/$key.hairpin.fa >>$od/select_mature.out");
		mkdir "$od/Known_miRNA_result" unless (-d "$od/Known_miRNA_result");
		system("cp $od/$key.hairpin.fa $od/Known_miRNA_result/Known_Pre_miRNA.fa");
		system("cp $od/$key.mature.fa $od/Known_miRNA_result/Known_Mature_miRNA.fa");
        stepTime(0);
        $step++ unless ($sbs);
}


################# Get all collapsed reads to identify Known miRNA
if($step==1){
        stepStart(1,"Get All Collapsed Reads");
        unless(defined $reads and -f $reads){
		foreach my $name (keys %Raw_FQ){
			my $fa = "$idir/$name/genome_map/$name.genome_map.fa";
			if (-f $fa){
				system("perl $Bin/bin/collapse_reads_md.pl $fa $name >$od/$name.collapsed.fa");
				system("cat $od/$name.collapsed.fa >> $od/Total_reads_collapsed.fa");
			}
			else {
				die "There is no file: $fa\nPlease Check Your Data.\n";
			}
		}
		$reads = "$od/Total_reads_collapsed.fa";
        }
        else {
		`ln -s $reads $od/Total_reads_collapsed.fa`;
		$reads = "$od/Total_reads_collapsed.fa";
        }

        stepTime(1);
        $step++ unless($sbs);
}

################## miRNA Expression Quantify
if($step==2){
        stepStart(2,"miRNA Expression Quantify");
	runOrDie("cd $od && perl $Bin/bin/quantifier.pl -p $od/Known_miRNA_result/Known_Pre_miRNA.fa -m $od/Known_miRNA_result/Known_Mature_miRNA.fa -r $reads -d -g 1 -e 2 -f 5 -W ");
	my $exp_file = glob "$od/miRNAs_expressed_all_samples_*.csv";
	unless (defined $exp_file and -f $exp_file) {
		die "There is no file: $od/miRNAs_expressed_all_samples_*.csv\nPlease Check Your Data.\n";
	}
	runOrDie("cd $od && perl $Bin/bin/expression_extra_v0.pl -exp $exp_file -fa $od/Known_miRNA_result/Known_Mature_miRNA.fa -pre $od/Known_miRNA_result/Known_Pre_miRNA.fa -oexp $od/Known_miRNA_result/All_Known_miRNA.count.list -ofa $od/Known_miRNA_result/All_Known_miRNA.expressed.fa -opre $od/Known_miRNA_result/All_Known_miRNA_Pre.expressed.fa ");
        stepTime(2);
        $step++ unless ($sbs);
}


################## Get Unquantified Reads
if($step==3){
        stepStart(3,"Get Unquantified Reads");
        my $quantified_file = glob "$od/expression_analyses/expression_analyses_*/Total_reads_collapsed_mapped.bwt";
        unless (defined $quantified_file and -f $quantified_file) {
		die "There is no file: $od/expression_analyses/expression_analyses_*/Total_reads_collapsed_mapped.bwt\nPlease Check Your Data.\n";
        }
	runOrDie("cd $od && perl $Bin/bin/abstract_fasta_seq_by_id.pl -i $quantified_file -fa $reads -o $od/unquantified_reads.fa -v 1 ");
        stepTime(3);
        $step++ unless ($sbs);
}


################################################################################################################
#sub functions
################################################################################################################

sub OutFileCheck {#检查输出文件路径是否为文件夹，如果是，则添加 _数字 后缀; 并会创建输出目录
	my $out_file = shift;
	my $out_file_check = $out_file;
	while (-d $out_file_check) {
		my $num ++;
		$out_file_check = $out_file."_$num";
	}
	my $out_dir = dirname($out_file_check);
	&MAKE_DIR($out_dir);
	return $out_file_check;
}

################################################################################################################

sub GetTime {
        my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
        return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################

sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
        my $cur_dir=`pwd`;chomp($cur_dir);
        my ($in)=@_;
        my $return="";
        if(-f $in){
                my $dir=dirname($in);
                my $file=basename($in);
                chdir $dir;$dir=`pwd`;chomp $dir;
                $return="$dir/$file";
        }elsif(-d $in){
                chdir $in;$return=`pwd`;chomp $return;
        }else{
                warn "Warning just for file and dir $in\n";
                exit;
        }
        chdir $cur_dir;
        return $return;
}

################################################################################################################

sub cmd_call {
        print "@_\n";
        system(@_) == 0 or die "system @_ failed: $?";
}

################################################################################################################

sub LOG_FORMAT {
        my $info = shift;
        my $star = '*'x80;
        my $time = &GetTime;
        my $error = $star."\n$time\n$info\n".$star."\n";
        return $error;
}

################################################################################################################

sub ERROR_AND_DIE {
        my $info = shift;
        my $error = &LOG_FORMAT($info);
        die "$error";
}

################################################################################################################

sub MAKE_DIR {
        my $directory = shift;
        if (-f $directory) {
                &ERROR_AND_DIE("$directory is a file!");
        }
        elsif (-d $directory) {
#               &ERROR_AND_DIE("Directory $directory exists!");
        }
        else {
                &cmd_call("mkdir -p $directory");
        }
        $directory = &ABSOLUTE_DIR($directory);
        return $directory;
}

################################################################################################################

sub check_array_no_same {#&check_array_no_same("array_name",@array_name); 检查数组中是否有相同的元素，有则die。
	my $array_name = shift;
	my @array = @_;
	foreach my $i (0..$#array-1) {
		foreach my $j ($i+1..$#array) {
			if ($array[$i] eq $array[$j]) {
				print "ERROR: \"$array[$i]\" appears twice in \@$array_name at least. Please check your input.\n";
				die "Illegal parameter input: -$array_name $array[$i]\n";
#				print "$i,$j:\tSame\n";
			}
		}
	}
}

################################################################################################################

sub USAGE {
	my $usage=<<"USAGE";
Program: $Script
Version: $version
Usage:
  Options:
    forced:
	-mature	 <file>	miRBase mature miRNA fa file, e.g /share/nas2/database/miRBase/v21/mature.fa
	-miR	 <str>	species or class to identify Known miRNA. This can be Latin abbreviations,e.g "vvi",
			OR species class: All,Viridiplantae,Metazoa,Viruses,Chromalveolata,Mycetozoa
	-od		output directory
	
	
    Must choose one:
	-idir	 <dir>	sRNA project Alignment dir
	    -maincfg	<file>	main.cfg of sRNA pipeline, must give this if you choose -idir
	 OR
	-reads	 <file>	All collapsed reads to identify Known miRNA, e.g Total_reads_collapsed.fa
	
	
    options:
	-species <str>	The species priority to retain, e.g  "hsa,mmu,rno"  or  "ahy"
	
	-step           <num> default 0
                0               Get Known Mature miRNA
                1               Get All Collapsed Reads
                2               miRNA Expression Quantify
                3               Get Unquantified Reads
                
        -onestep        only one step
	
	-h		Help

USAGE
	print $usage;
	exit;
}

