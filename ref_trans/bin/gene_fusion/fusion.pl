#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($gff,$indir,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"gff:s"=>\$gff,
				"fusion:s"=>\$indir,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($gff and $indir and $od);


mkdir $od unless -d $od;
$od=&ABSOLUTE_DIR($od);
###############################################################################################################
`cp $Bin/readme.txt $od`;

#--------------------------------------------gff   file   -----------------------------------------------------
my %GTF;
my ($chr,$symbol,$start,$end,$name);
open (IN,$gff) or die $!;
while(<IN>){
	chomp;
	next if (/^\s*$/ || /^\#/);
	my @tmp=split("\t",$_);
	next if($tmp[2] !~ /gene/);
	($chr,$symbol,$start,$end,$name)=($tmp[0],$tmp[2],$tmp[3],$tmp[4],$tmp[8]);
	#($chr,$symbol,$start,$end,$name)=($tmp[0],$tmp[2],$tmp[4],$tmp[3],$tmp[8]) if($tmp[3]<$tmp[4]);
	$name=~ s/ID=//;
	push@{$GTF{$chr}{$name}},$start;
	push@{$GTF{$chr}{$name}},$end;
}
close(IN);
#print Dumper %GTF;
#die;
###############################################################################################################


#--------------------------------------------------------------------------------------------------------------
my ($fusion,$pos1,$pos2,$orient,$depth,$chr1,$chr2,@sub1,@sub2,$len1,$len2);
my @file=glob "$indir/*/fusions.out";
my $base;
foreach my $in (@file) {
    $base=dirname($in);
    $base=(split("/",$base))[-1];
    open IN,$in or die $!;
    open OUT,">$od/$base.txt" or die $!;
    print OUT "Fusion\torientations\tGeneID1\tChromosome\tFusion_position\tGeneID2\tChromosome\tFusion_position\treads_span_fusion\n";
    while(<IN>){
        chomp;
        my @tmp=split("\t",$_);
        ($fusion,$pos1,$pos2,$orient,$depth)=($tmp[0],$tmp[1],$tmp[2],$tmp[3],$tmp[4]);
        next if($depth<100);
        ($chr1,$chr2)=split("-",$fusion);
        next if(!exists $GTF{$chr1});
        next if(!exists $GTF{$chr2});
        @sub1=();
        @sub2=();
        foreach my $gene1 (keys %{$GTF{$chr1}}) {
	  my $a=(@{$GTF{$chr1}{$gene1}})[0];
	  my $b=(@{$GTF{$chr1}{$gene1}})[1];
	  if($pos1 >=$a && $pos1<=$b){
	      push @sub1,$gene1;
	  }
        }
        $len1=@sub1;
        if($len1 !=0){
	  print OUT "$fusion\t$orient\t",join(";",@sub1),"\t",$chr1,"\t",$pos1,"\t";
        }else{
	  print OUT "$fusion\t$orient\t","-","\t",$chr1,"\t",$pos1,"\t";
        }
        foreach my $gene2 (keys %{$GTF{$chr2}}) {
	  my $c=(@{$GTF{$chr2}{$gene2}})[0];
	  my $d=(@{$GTF{$chr2}{$gene2}})[1];
	  if($pos2 >= $c && $pos2 <= $d){
	      push @sub2,$gene2;
	  }
        }
        $len2=@sub2;
        if($len2 !=0){
	  print OUT join(";",@sub2),"\t",$chr2,"\t",$pos2,"\t",$depth,"\n";
        }else{
	  print OUT "-","\t",$chr2,"\t",$pos2,"\t",$depth,"\n";
        }
    }
    close IN;
    close OUT;
}





###############################################################################################################
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
sub QSUB {
	my $file = shift;
	chomp(my $host = `hostname`);
	if ($host =~ /cluster/) {
		&cmd_call("qsub-sge.pl $file --maxproc 50  --resource vf=15G --reqsub --independent");
	}else{
		&cmd_call("ssh cluster qsub-sge.pl $file  --maxproc 50  --resource vf=15G --reqsub --independent");
	}
}

sub cmd_call {
	print "@_\n";
	system(@_) == 0 or die "system @_ failed: $?";
}

sub PARA_CONFIG {
	my $hash = shift;
	my $fcfg = shift;
	open IN,"$fcfg" or die "cann open $fcfg, $!";
	while (<IN>) {
		chomp;
		next if (/^\#/||/^\s*$/);
		$_ =~ s/\#.*$//;
		chomp $_;
		my @tmp = split /\s+/,$_;
		$$hash{$tmp[0]} = $tmp[1];
	}
	close IN;
}
sub check_time{
	my $dir=shift;
	my $check;
	my @Old_time=(0,0,0,0,0,0);
	foreach my $d (glob "$dir/run*") {
		next unless -d $d;
		my $run_name=basename $d;
		my @Time=(split/_/,$run_name)[3,2,1,5,6,7];
		foreach my $i (0..5) {
			if ($Time[$i]>$Old_time[$i]) {
				@Old_time=@Time;last;
			}
		}
	}
	return "run_$Old_time[2]_$Old_time[1]_$Old_time[0]_t_$Old_time[3]_$Old_time[4]_$Old_time[5]";
}

#sub output_expression{
#	my $hash=shift;
#	my $type=shift;
#	foreach my $key (keys %{$$hash{$type}}) {
#		my $index=
#	}
#}

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	Wang Yajing <wangyj\@biomarker.com.cn> 
Usage:
  Options:
  -fusion <dir>    Tophat dir 
  -gff <file>      Genome GFF
  -od  <dir>    output dir,forced 
  

USAGE
	print $usage;
	exit;
}
