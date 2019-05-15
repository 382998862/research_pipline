#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0";

my($miRBase_mature,$idir,$main_cfg,$reads,$key,$species,$od,$step,$sbs,$min,$max);
GetOptions(
	"help|?"	=>\&USAGE,
	"mature:s"	=>\$miRBase_mature,
	"idir:s"	=>\$idir,
	"maincfg:s"	=>\$main_cfg,
	"reads:s"	=>\$reads,
	"miR:s"		=>\$key,
	"species:s"	=>\$species,
	"od:s"		=>\$od,
	"min:s"		=>\$min,
	"max:s"		=>\$max,
	"step:s"	=>\$step,
	"onestep"	=>\$sbs,
) or &USAGE;
&USAGE unless ($miRBase_mature and $key and $od and $main_cfg);
&USAGE unless ($idir or $reads);
&USAGE if ($idir and $reads);


$miRBase_mature = &ABSOLUTE_DIR($miRBase_mature);
$od = &MAKE_DIR($od);

$step ||= 0;
$species ||= "NONE";
$min ||= 18;
$max ||= 30;

############### read CFG
my %CFG=%{readconf("$Bin/../CFG")};
$miRBase_mature ||= $CFG{mature};

################################################################################################################
#main code
################################################################################################################

################## Creat miRBase file
if($step==0){
	stepStart(0,"Creat miRBase file");
	runOrDie("perl $Bin/bin/select_fa.pl -i $key -fa $miRBase_mature -o $idir/Ref_Database/$key.mature.fa >>$idir/Ref_Database/select_mature.out");
	stepTime(0);
	$step++ unless ($sbs);
}



################## Known miRNA Identify
if($step==1){
        stepStart(1,"Known miRNA Identify");
	if (defined $idir and -d $idir) {
		runOrDie("perl $Bin/Known_miRNA_blast/Known_miRNA_blast.pl -miR $key -mature $miRBase_mature -od $od/Known_miRNA_identify -idir $idir -maincfg $main_cfg -species $species ");
	}elsif (defined $reads and -f $reads) {
		runOrDie("perl $Bin/Known_miRNA_blast/Known_miRNA_blast.pl -miR $key -mature $miRBase_mature -od $od/Known_miRNA_identify -reads $reads -species $species ");
	}
        stepTime(1);
        $step++ unless ($sbs);
}


################# Novel miRNA prediction
if($step==2){
        stepStart(2,"Novel miRNA prediction");
	my $uq_reads = (glob "$od/Known_miRNA_identify/unquantified_reads.fa")[0];
	if (defined $uq_reads && -f $uq_reads) {
		runOrDie("perl $Bin/Novel_miRNA_prediction.pl -uq $uq_reads -idir  $od -main_cfg $main_cfg -min $min -max $max ");
	}else {
		if (defined $idir and -d $idir) {
			runOrDie("perl $Bin/basic_analysis_without_miRBase.pl -idir $idir -main_cfg $main_cfg -min $min -max $max ");
		}
	}
	
        stepTime(2);
        $step++ unless($sbs);
}


if($step==3){
	my @sums=glob("$od/miRDeep2/miRNA_Quantify/sum_*.txt");
	my $sum=join(" ",@sums);
	stepStart(3,"All miRNA location");
	`mkdir $od/miRNA_loc`	unless(-d "$od/miRNA_loc");
	runOrDie("perl $Bin/get_miRNA_position/extract_miRNA_position.pl -i $sum -k $Bin/get_miRNA_position/know_miRNA_pos.list -e $od/miRDeep2/miRNA_Quantify/All_miRNA_expression.list -o $od/miRNA_loc");
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
	-od		output directory, Analysis of sRNA project
	-maincfg <file> main.cfg of sRNA pipeline, must give this if you choose -idir
	-genome	 <file>	reference genome file, for novel miRNA prediction
	-plant	 <int>	animal 0; plant 1;
	-max_step<int>	default 5	
	
    Must choose one:
	-idir	 <dir>	sRNA project Analysis directory, contains Alignment dir
	 OR
	-reads	 <file>	All collapsed reads to identify Known miRNA, e.g Total_reads_collapsed.fa
	
	
    options:
	-species <str>	The species priority to retain, e.g  "hsa,mmu,rno"  or  "ahy"
	-min		the min length of predicted miRNA, default 18
	-max		the max length of predicted miRNA, default 30
	
	-step           <num> default 0
                0               Creat miRBase file
                1               Known miRNA Identify
                2               Novel miRNA prediction
                
        -onestep        only one step
	
	-h		Help

USAGE
	print $usage;
	exit;
}

