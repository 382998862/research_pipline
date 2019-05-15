#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my ($outDir,$data_cfg);
GetOptions(
                "help|?" =>\&USAGE,
                "o:s"=>\$outDir,
				"data_cfg:s"=>\$data_cfg
                ) or &USAGE;
&USAGE unless ($data_cfg and $outDir);
$data_cfg=&ABSOLUTE_DIR($data_cfg);
my %Samples=();
&parse_sample_config ($data_cfg,\%Samples);

open (OUT,">$outDir/All.fq") or die $!;
my $flag = 1;
foreach my $sam (sort keys %Samples) 
{
	my $fq = $Samples{$sam}{FQ};
	open (IN,"$fq")or die $!;
	open (COUNT,">$outDir/$sam.count")or die $!;
	my $count = 0;
	while (<IN>) {
		next if(/^\s*$/);
		if($flag%2!=0 && /^@/)
		{
			$count++;
			my $line = $_;
			$line =~s/^(@)/\@$sam:/;
			print OUT $line;
		
		}
		else
		{
			print OUT $_;
		}
		$flag ++;
	}
	
	print COUNT "$count";
	close COUNT;
	close IN;
}
close OUT;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
####################subs
sub parse_sample_config
{ #load config file
	my $config_file= shift;
	my $SamData= shift;
	
	my $error_status = 0;
	open IN,$config_file || die "fail open: $config_file";
	while (<IN>) {
		chomp;
		s/\s+$//;s/\r$//;
		next if(/^$/ or /^\#/);
		if ( /^Sample/ ) {
			my $fq_info = <IN> ; chomp $fq_info;
			my $sample_name = (split /\s+/, $_)[1];
			my $fq = (split /\s+/, $fq_info)[1];

			if (! -e $fq) {
				warn "Non-exist: Reads file dir for $sample_name isn't exists\n" ;
				exit;
			}
			$SamData->{$sample_name}->{FQ} = $fq;
		}
	}
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

################################################################################################################
sub USAGE {
    my $usage=<<"USAGE";
 ProgramName:
     Version:   $version
     Contact:   Simon Young <yangxh\@biomarker.com.cn> 
Program Date:   2012.07.02
      Modify:   
 Description:   This program is used to ......
       Usage:
        Options:
        -data_cfg <file>   data config file, forced
        -o <file>   output dir,optional

        -h      help

USAGE
    print $usage;
    exit;
}

