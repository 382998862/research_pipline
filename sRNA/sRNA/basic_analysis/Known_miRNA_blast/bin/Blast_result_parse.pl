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

my($fa,$db,$blast_result,$species,$out);
GetOptions(
	"help|?"	=>\&USAGE,
	"fa:s"		=>\$fa,
	"db:s"		=>\$db,
	"bo:s"		=>\$blast_result,
	"species:s"	=>\$species,
	"out:s"		=>\$out,
) or &USAGE;
&USAGE unless ($fa and $db and $blast_result and $out);
$fa = &ABSOLUTE_DIR($fa);
$db = &ABSOLUTE_DIR($db);
$blast_result = &ABSOLUTE_DIR($blast_result);

################################################################################################################
#main code
################################################################################################################
#The species priority to retain 
$species ||= "NONE";
my @sp = reverse(split ",",$species);


#The seq length of reads
my %len;
$/ = ">";
open FA,$fa;
while (<FA>) {
	chomp;
	next if (/^#/ || /^\s*$/);
	my ($id,$seq) = split "\n",$_,2;
	$seq =~ s/\n//g;
	my $len = length($seq);
	$len{$id} = $len;
}
close FA;

#The seq length of known miRNAs
my %dblen;
open DB,$db;
while (<DB>) {
	chomp;
	next if (/^#/ || /^\s*$/);
	my ($id,$seq) = split "\n",$_,2;
	$seq =~ s/\n//g;
	my $len = length($seq);
	$dblen{$id} = $len;
}
close DB;


#filter blast results
my %blast;
$/ = "\n";
open IN,$blast_result;
while (<IN>) {
	chomp;
	my @line = split "\t";
	next unless $line[3] == $len{$line[0]};
	next unless $line[3] == $dblen{$line[1]};
	next unless $line[2] == 100;
	next unless $line[6] == 1;
	next unless $line[8] == 1;
	$blast{$line[0]}{$line[1]} = 1;
#	print "$_\n";
}
close IN;


open OUT,">$out";
print OUT "#Read_id\tKnown_miRNA\n";
#my %final_reads_2_db; #only one to one.
foreach my $reads_id (sort keys %blast) {
	my $remain_id = "";
	foreach my $db_id (sort keys %{$blast{$reads_id}}) {
		$remain_id ||= $db_id;
		foreach my $sp (@sp) {
			$remain_id = $db_id if $db_id =~ /^$sp/;
		}
	}
#	$final_reads_2_db{$reads_id} = $remain_id;
	print OUT "$reads_id\t$remain_id\n";
}
close OUT;



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
	-fa	 <file>	Total reads collapsed fa file, Total_reads_collapsed.fa
	-db	 <file>	blast database fa file
	-bo	 <file>	blast result file
	-out		output file
	
    options:
	-species <str>	The species priority to retain, e.g  "hsa,mmu,rno"  or  "ahy"
	
	-h		Help

USAGE
	print $usage;
	exit;
}

