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

my($unquantified_reads,$od);
GetOptions(
	"help|?"	=>\&USAGE,
	"uq:s"		=>\$unquantified_reads,
	"od:s"		=>\$od,
) or &USAGE;
&USAGE unless ($unquantified_reads and $od);
$unquantified_reads = &ABSOLUTE_DIR($unquantified_reads);
$od = &MAKE_DIR($od);

################################################################################################################
#main code
################################################################################################################
my %sample_reads;
$/ = ">";
open FA,$unquantified_reads;
while (<FA>) {
	chomp;
	next if (/^#/ || /^\s*$/);
	my ($id,$seq) = split "\n",$_,2;
	$seq =~ s/\n//g;
	print STDERR "ID of ypur reads should be this format: >SampleName_\\d+_x\\d+ .\nThe read with ID $id will be abandoned." unless $id =~ /^([^\s_]+)_(\d+)_x(\d+)$/;
	my ($sample_name,$reads_number) = ($1,$3);
	$sample_reads{$sample_name}{$seq} = $reads_number;
}
close FA;


foreach my $sample (sort keys %sample_reads) {
	my $index = 0;
	open OUT,">$od/$sample.mapper.fa";
	foreach my $seq (keys %{$sample_reads{$sample}}) {
		for (1..$sample_reads{$sample}{$seq}) {
			$index++;
			print OUT ">seq_$index\n$seq\n";
		}
	}
	close OUT;
}


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

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
	-uq	 <file>	"unquantified_reads.fa", fasta file. The ID should be this format: >SampleName_\\d+_x\\d+ 
	-od	 <dir>	output directory, e.g "miRDeep2"
	
	-h		Help

USAGE
	print $usage;
	exit;
}

