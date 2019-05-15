use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);
my ($genegff,$lncgff,$type,$output,$cis);
GetOptions(
        "h|?"           =>\&USAGE,
	"lncgff:s"	=>\$lncgff,
	"genegff:s"	=>\$genegff,
	"type:s"     	=>\$type,
	"output:s"	=>\$output,
	"cis:s"		=>\$cis,
)or &USAGE;
&USAGE unless ($genegff and $lncgff  and $output);

$cis||="100000;10000;1000";
$cis  ="100000;10000;1000"	if($cis eq "100000");
my @cis=split(/\;/,$cis);
push @cis,0	if($cis[-1] !=0);

open(GFF,$genegff)||die $!;
my %genes=();
while(<GFF>){
	chomp;
	my ($chr,$source,$typ,$start,$end,$t1,$strand,$t2,$info)=split(/\t/,$_);
	next if($typ ne "gene");
	my $id=(split(/\;/,$info))[0];
	$id=~s/ID=//;
	$genes{$chr}{$start}{$end}=$id;
}
close(GFF);

my %code=();
if(defined $type){
	open(TYPE,$type)||die $!;
	while(<TYPE>){
		chomp;
		my ($id,$type)=split(/\t/,$_,2);
		$code{$id}=$type;
	}
	close(TYPE);
}


my @set=&set(\@cis);
print join("\t",@set),"\n";
open(GFF,$lncgff)||die $!;
open(OUT,">$output")||die $!;
my $header="#ID";
$header.="\tType"	if(defined $type);
print OUT "$header\tTotal\tOverlap\t",join("\t",map{"$_(upstream)\t$_(downstream)"} @set),"\n";
my %test=();
while(<GFF>){
	chomp;
	my ($chr,$source,$typ,$start,$end,$t1,$strand,$t2,$info)=split(/\t/,$_);
	if($typ =~/mRNA|^lincRNA$|^transcript|lncRNA$/){
	my $id=(split(/\;/,$info))[0];
	$id=~s/ID=//;
	$id=~s/transcript://;
	my @total=();my @over=();
	foreach my $s (sort{$a <=>$b} keys %{$genes{$chr}}){
		next if($s>$end+$cis[0]);
		foreach my $e(sort{$a<=>$b} keys %{$genes{$chr}{$s}}){
			next if($e<$start-$cis[0]);

			push @total,$genes{$chr}{$s}{$e};
			if(($s>=$start && $s<=$end)||($e>=$start && $e<=$end) ||($s<$start && $e>$end)){
				push @over,$genes{$chr}{$s}{$e};
			}else{
				my $loc=&relativeLoc($start,$end,$s,$e);
				my $f=0;
				for(my $i=0;$i<@cis-1;$i++){
					next if($f==1);
					if(abs($loc)>$cis[$i+1] && abs($loc)<=$cis[$i]){
						my $flag=($loc>0)?"(downstream)":"(upstream)";
						my $t=($cis[$i+1]/1000)."-".($cis[$i]/1000)."Kb".$flag;
						$test{$t}++;
						${$id}{$t}{$genes{$chr}{$s}{$e}}=abs($loc);
						$f=1;
					}
				}
				if($f==0){
					print "$id\t$genes{$chr}{$s}{$e}\n";
				}
								
			}
		}
	}
	if(@total>0){
		my $out=$id;
		$out .="\t$code{$id}" if(defined $type);
		$out .="\t".join(";",@total)."\t".join(";",@over);
		foreach my $s(@set){
		#	my %up		=keys %{${$id}{$s."(upstream)"}};
		#	my %down	=keys %{${$id}{$s."(downstream)"}};
		#	my @up_rnas	=sort{$up{$a} <=> $up{$b}} keys %up;
		#	my @down_rnas	=sort{$down{$a} <=> $down{$b}} keys %down;
			my @up_rnas	=sort {$a<=>$b} keys %{${$id}{$s."(upstream)"}};
			my @down_rnas	=sort {$a<=>$b} keys %{${$id}{$s."(downstream)"}};
			$out .="\t".join(";",@up_rnas)."\t".join(";",@down_rnas);
		}
		print OUT "$out\n";
	}
	}
	
}
close(OUT);
close(GFF);
print join("\t",keys %test),"\n";
sub set{
	my $cis=shift;
	my @cis=@{$cis};
	my @out=();
	for(my $i=0;$i<@cis-1;$i++){
		my ($t1,$t2)=($cis[$i+1]/1000,$cis[$i]/1000);
		push @out,"$t1-$t2"."Kb";
	}
	return reverse @out;
}


sub relativeLoc{
	my ($s1,$e1,$s2,$e2)=@_;	##s1,e1 =>lncRNA     s1,e2 =>mRNA
	my $loc;
	if($e2<$s1){
		$loc=$e2-$s1;
	}elsif($s2>$e1){
		$loc=$s2-$e1;
	}
	return $loc;
}


sub readConfig{
	my $configFile=shift;
	my $d=Config::General->new(-ConfigFile => "$configFile");
	my %config=$d->getall;	
	return %config;
}
sub USAGE{
	my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
	Options:
	-lncgff		<file>	input lncRNA gff file, forced
	-genegff	<file>	input gene gff file, forced
	-type		<file>	lncRNA type file, contained two columns eg:#ID	type, not must
	-output		<file>	output file, forced
	-cis		<str>	default 100000,10000,1000
				means the Cis file would be divided into 3 conditions
				eg: 0-1K,1-10K,10-100k


	-h	Help

Example:

USAGE
	print $usage;
	exit;
}


