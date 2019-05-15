#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$samples);
GetOptions(
				"help|?" =>\&USAGE,
				"id:s"=>\$fIn,
				"samples:s"=>\$samples,
				) or &USAGE;
&USAGE unless ($fIn and $samples);

$fIn=&ABSOLUTE_DIR($fIn);

foreach my $sam (split/,/,$samples) {
	my $dir="$fIn/$sam";
	mkdir $dir unless (-d $dir);
	mkdir "$dir/Len_stat" unless -d "$dir/Len_stat";
	#################################################################clean data
	runOrDie("perl $Bin/length_stat.pl -i $dir/$sam.clean.fa  -od $dir/Len_stat -prefix $sam.clean_len");
	#################################################################maped genome data
	runOrDie("perl $Bin/length_stat.pl -i $dir/genome_map/$sam.genome_map.fa  -od $dir/Len_stat -prefix $sam.map_genome");
	my %Len;
	my %Total;
	open (IN,"$dir/Len_stat/$sam.clean_len.stat") or die $!;
	<IN>;
	while (<IN>) {
		chomp;
		next if /^\s*$/;
		my ($len,$num)=split/\s+/,$_;
		$Len{$len}{Clean}=$num;
		$Total{Clean}+=$num;
	}
	close IN;
	open (IN,"$dir/Len_stat/$sam.map_genome.stat") or die $!;
	<IN>;
	while (<IN>) {
		chomp;
		next if /^\s*$/;
		my ($len,$num)=split/\s+/,$_;
		$Len{$len}{Genome}=$num;
		$Total{Genome}+=$num;
	}
	close IN;
	open (OUT,">$dir/Len_stat/$sam.Total.stat") or die $!;
	print OUT "#Length\tClean_num\tGenome_num\tGenome_per\n";
	foreach my $key (sort {$a<=>$b} keys %Len) {
		$Len{$key}{Genome}||=0;
		my $per_Genome=sprintf "%.2f",$Len{$key}{Genome}/$Len{$key}{Clean}*100;
		print OUT "$key\t$Len{$key}{Clean}\t$Len{$key}{Genome}\t$per_Genome%\n";
	}
	$Total{Genome}||=0;
	my $per_Genome=sprintf "%.2f",$Total{Genome}/$Total{Clean}*100;
	print OUT "Total\t$Total{Clean}\t$Total{Genome}\t$per_Genome%\n";
	close OUT;
}




#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
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
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

################################################################################################################

sub max{#&max(lists or arry);
	#求列表中的最大值
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}

################################################################################################################

sub min{#&min(lists or arry);
	#求列表中的最小值
	my $min=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$min=$min<$temp?$min:$temp;
	}
	return $min;
}

################################################################################################################

sub revcom(){#&revcom($ref_seq);
	#获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return uc $seq;			  
}

################################################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	Zhang XueChuan <zhangxc\@biomarker.com.cn> 
Usage:
  Options:
  -id         <dir>  input dir,forced 
  
  -samples    <str>  samples,"S01,S02",forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
