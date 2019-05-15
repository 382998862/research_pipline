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

my($miRBase,$blast_parsed_result,$od);
GetOptions(
	"help|?"	=>\&USAGE,
	"miRBase:s"	=>\$miRBase,
	"bp:s"		=>\$blast_parsed_result,
	"od:s"		=>\$od,
) or &USAGE;
&USAGE unless ($blast_parsed_result and $od);
$miRBase ||= "$Bin/miRNA.txt";
$miRBase = &ABSOLUTE_DIR($miRBase);
$blast_parsed_result = &ABSOLUTE_DIR($blast_parsed_result);
$od = &MAKE_DIR($od);

#####
# Note: Take all pre of one mature miRNA.
#####

################################################################################################################
#main code
################################################################################################################
# get all miRNA id to be taken
my %miRNA;
open IN,"$blast_parsed_result";
while (<IN>) {
	chomp;
	next if (/^#/ or /^\s*$/);
	my $miRNA_id = (split "\t",$_)[1];
	$miRNA{$miRNA_id} = 1;
}
close IN;


# take all pre and mature. One mature to N pre. Maybe one pre to N mature.
open MIR,$miRBase;
open PRE,">$od/Known_Pre_miRNA.fa";
while (<MIR>) {
	chomp;
	next if (/^#/ or /^\s*$/);
	my @line = split "\t";
	my $miRNA_num = (scalar(@line) - 4)/3;
	my $check = 0;
	for my $i (1..$miRNA_num) {
		my $id_index = $i*3+2;
		my $miRNA_seq_index = $id_index + 1;
		if (defined $miRNA{$line[$id_index]}) {
			$line[$miRNA_seq_index] =~ tr/VDBHWSKMYRvdbhwskmyrUu/AATAACTACAaataactacaTt/;
			$miRNA{$line[$id_index]} = $line[$miRNA_seq_index];
			$check = 1;
		}
		
	}
	if ($check == 1) {
		$line[3] =~ tr/VDBHWSKMYRvdbhwskmyrUu/AATAACTACAaataactacaTt/;
		print PRE ">$line[1]\n$line[3]\n";
	}
}
close PRE;


open MATURE,">$od/Known_Mature_miRNA.fa";
foreach my $miRNA_id (sort keys %miRNA) {
	print MATURE ">$miRNA_id\n$miRNA{$miRNA_id}\n";
}
close MATURE;



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
	-bp	 <file>	blast parsed result file, col 2 is miRNA ID
	-od	 <dir>	output directory
	
    options:
	-miRBase <file>	pre-miRNA and mature miRNA relationship, tabular by tab, see miRNA.xls of miRBase
			default: $Bin/miRNA.txt
	
	-h		Help

USAGE
	print $usage;
	exit;
}

