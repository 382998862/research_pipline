#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($i,$fa,$o);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$i,
				"fa:s"=>\$fa,
				"o:s"=>\$o,
				) or &USAGE;
&USAGE unless ($i and $fa and $o);
#my %baseCompRegPartern=(
#	"V"=>"[ACG]",
#	"D"=>"[ATG]",
#	"B"=>"[TCG]",
#	"H"=>"[ATC]",
#	"W"=>"[AT]",
#	"S"=>"[CG]",
#	"K"=>"[TG]",
#	"M"=>"[AC]",
#	"Y"=>"[CT]",
#	"R"=>"[AG]",
#);

$/=">";
open (OUT,">$o") or die $!;
open (IN,$fa) or die $!;
<IN>;
while (<IN>) {
	chomp;
	my $id=(split/\s+/,$_)[0];
	my $info=(split/\n/,$_,2)[1];
	$info=~tr/VDBHWSKMYRvdbhwskmyr/AATAACTACAaataactaca/;
	$info=~tr/atcguU/ATCGTT/;
	$info =~ s/\n//g;
	if ($i!~/All/) {
		if ($i=~m/Chromalveolata|Metazoa|Mycetozoa|Viridiplantae|Viruses/){
			my $system_LATIN=&system_Latin($i);
			next if ($id!~/^($system_LATIN)/);
		}else{
			next if ($id!~/^$i/) ;
		}
	}
	print OUT ">$id\n$info\n";
}
close IN;
close OUT;
$/="\n";





sub system_Latin{
	my ($system)=@_;
	my $Litan;
	if ($system=~m/Chromalveolata/i) {
		$Litan="esi|pti|pin|pra|psj";
	}elsif($system=~m/Metazoa/i) {
		$Litan="aqu|nve|hma|sko|spu|cin|csa|odi|bfl|xla|xtr|gga|cfa|mdo|age|lla|sla|mml|mne|pbi|ggo|hsa|ppa|ppy|ptr|ssy|lca|oan|cgr|mmu|rno|bta|oar|ssc|dre|fru|tni|aga|ame|bmo|dan|der|dgr|dme|dmo|dpe|dps|dse|dsi|dvi|dwi|dya|lmi|tca|dpu|isc|cbr|cel|sja|sma|sme|cte|cla|hru|lgi|crm|ppc|eca|ola|bma|api|aae|cqu|tgu|nvi|smr|ngi|nlo|pma|xbo|egr|emu|hme|gpy|tre|rmi|asu|aca|meu|sha|ocu|aja|efu|chi|tch|ccr|hhi|ipu|pol|ssa|oha|lva|pmi|tur|mja|mse|pxy|cbn|hco|prd|str|gsa|lco|sci";
	}elsif($system=~m/Mycetozoa/i) {
		$Litan="ddi";
	}elsif($system=~m/Viridiplantae/i) {
		$Litan="cre|pta|ppt|smo|ath|bna|bol|bra|cpa|gma|lja|mtr|vun|ghb|ghr|gra|ptc|sly|vvi|osa|sbi|sof|tae|zma|pvu|mdm|bdi|aqc|peu|csi|ccl|crt|ctr|rco|gar|aly|ahy|gso|pab|ttu|ata|hvu|far|bgy|bcy|tcc|rgl|cme|ssp|amg|aau|ssl|stu|pde|atr|pgi|cca|han|har|hci|hex|hpa|hpe|htu|hbr|mes|ama|dpr|lus|ppe|nta|egu|cln";
	}elsif($system=~m/Viruses/i) {
		$Litan="bkv|ebv|hcmv|hiv1|hsv1|jcv|kshv|mcmv|mcv|mdv1|mdv2|mghv|rlcv|rrv|sv40|hsv2|hbv|iltv|hvt|bhv1|bpcv1|bpcv2|hvsa|bfv|bhv5|blv|dev|hhv6b|prv";
	}
#	my @Litans = split ('|',$Litan);
#	$Litan = join ('^|',@Litans);
	return $Litan;
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
Program Date:   2013.12.26
Usage:
  Options:
  -i   <str>   keyword,special choice 'All','Viridiplantae','Metazoa','Viruses','Chromalveolata','Mycetozoa',forced 
  
  -fa  <file>  input file,fasta format,forced 
  
  -o   <file>  output file,fasta format,forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
