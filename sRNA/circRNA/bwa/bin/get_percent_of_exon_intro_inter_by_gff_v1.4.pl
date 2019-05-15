#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
#my $version="1.0.0";

my %config=%{readconf("$Bin/../../project.cfg")};
my $Title=$config{Title};												#流程的名称，必填
my $version=$config{version};
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$gff,$index,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"gff:s"=>\$gff,
				"index:s"=>\$index,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($fIn and $gff and $index and $od);

mkdir $od unless -d $od;
$od=&ABSOLUTE_DIR($od);

my %ST;
my %Gene;
my %Exon;
my %Len;
my $geneids;
open (IN,$gff) or die $!;
while (<IN>) {
	chomp;
	next if (/\#/||/^$/) ;
	my ($chr,$type,$start,$end,$strand,$ID)=(split /\s+/,$_)[0,2,3,4,6,8];
	if ($type eq 'gene') {
		$geneids=(split/\;|=/,$ID)[1];
		$Gene{$chr}{$start}{$end}=$geneids;
		$ST{$geneids}=1 if $strand eq '-';
	}
	if ($type eq 'CDS'||$type eq 'exon') {
		if (exists $Exon{$geneids}{$start}) {
			$Exon{$geneids}{$start}=$end if $end>$Exon{$geneids}{$start};
		}
		else{$Exon{$geneids}{$start}=$end}
	}
}
close (IN) ;

foreach my $chr (keys %Gene) {
	foreach my $g_s (sort {$a<=>$b} keys %{$Gene{$chr}}) {
		foreach my $g_e (keys %{$Gene{$chr}{$g_s}}) {
			my $gene_name=$Gene{$chr}{$g_s}{$g_e};
			my $limit=keys %{$Exon{$gene_name}};
			while (1) {
				my $o_l=$limit;
				my $before_s=-2;
				my $before_e=-1;
				foreach my $e_s (sort {$a<=>$b} keys %{$Exon{$gene_name}}) {
					if ($before_e<$e_s) {
						$before_s=$e_s;
						$before_e=$Exon{$gene_name}{$e_s};
					}
					else{
						my $new_s=($before_s>$e_s)?$e_s:$before_s;
						my $new_e=($before_e>$Exon{$gene_name}{$e_s})?$before_e:$Exon{$gene_name}{$e_s};
						delete $Exon{$gene_name}{$before_s};
						delete $Exon{$gene_name}{$e_s};
						$Exon{$gene_name}{$new_s}=$new_e;
						$limit--;
						last;
					}
				}
				last if $o_l==$limit;
			}
			if ($limit>1) {
				my $length;
				foreach my $e_s (sort {$a<=>$b} keys %{$Exon{$gene_name}}) {
					$length+=$Exon{$gene_name}{$e_s}+1-$e_s;
				}
				$Len{$gene_name}=$length;
			}
			if ($limit==1) {
				foreach my $e_s (keys %{$Exon{$gene_name}}) {
					$Len{$gene_name}=$Exon{$gene_name}{$e_s}+1-$e_s;
				}
			}
		}
	}
}




#########先无视最后一个基因后面的基因间区
my %T;
foreach my $chr (keys %Gene) {
	my $b_e=0;
	foreach my $g_s (sort {$a<=>$b} keys %{$Gene{$chr}}) {
		foreach my $g_e (keys %{$Gene{$chr}{$g_s}}) {
			my $gene=$Gene{$chr}{$g_s}{$g_e};
			my $inter_s=$b_e+1;
			my $inter_e=$g_s-1;
			if ($inter_s<$inter_e) {
				$T{$chr}{$inter_s}{end}=$inter_e;
				$T{$chr}{$inter_s}{type}="Intergenic";
				$T{$chr}{$inter_s}{value}=0;
			}
			$b_e=$g_e;
			my $b_e_e=$g_s-1;
			foreach my $s (sort {$a<=>$b} keys %{$Exon{$gene}}) {
				my $i_s=$b_e_e+1;
				my $i_e=$s-1;
				if ($i_s<$i_e) {
					$T{$chr}{$i_s}{gene}=$gene;
					$T{$chr}{$i_s}{end}=$i_e;
					$T{$chr}{$i_s}{type}="Intron";
					$T{$chr}{$i_s}{value}=0;
				}
				$b_e_e=$Exon{$gene}{$s};
				$T{$chr}{$s}{gene}=$gene;
				$T{$chr}{$s}{end}=$b_e_e;
				$T{$chr}{$s}{type}="Exon";
				$T{$chr}{$s}{value}=0;
			}
			$b_e_e++;
			if ($b_e_e<$g_e) {
				$T{$chr}{$b_e_e}{gene}=$gene;
				$T{$chr}{$b_e_e}{end}=$g_e;
				$T{$chr}{$b_e_e}{type}="Intron";
				$T{$chr}{$b_e_e}{value}=0;
			}
		}
	}
}

my %F;
foreach my $chr (keys %T) {
	foreach my $s (sort {$a<=>$b} keys %{$T{$chr}}) {
		push @{$F{$chr}},$s;
	}
}

my %Original;
my $chr_old="";
my $n=0;
open (IN,"$fIn") or die $!;
while (<IN>) {
	chomp;
	next if /^\s*$/;
	my @A=split /\s+/,$_;
	if ($A[0] ne $chr_old) {
		$n=0;
		$chr_old=$A[0];
	}
	if ($A[0] eq $chr_old) {
		next unless exists $F{$A[0]};
		my $limit=@{$F{$A[0]}};
		next if $n>=$limit;
		my $start=$F{$A[0]}->[$n];
		my $end=$T{$A[0]}{$start}{end};
		while ($A[1]>$end) {
			$n++;
			last if $n>=$limit;
			$start=$F{$A[0]}->[$n];
			$end=$T{$A[0]}{$start}{end};
			last if $A[1]<=$end;
		}
		next if $n>=$limit;
		$T{$A[0]}{$start}{value}+=$A[2];
		if ($T{$A[0]}{$start}{type} eq 'Exon') {
			$Original{$T{$A[0]}{$start}{gene}}{site}{$A[1]}+=$A[2];
			$Original{$T{$A[0]}{$start}{gene}}{total}+=$A[2];
		}
	}
}
close IN;



my ($e_n,$i_n,$b_n)=(0,0,0);
open (OUT,">$od/$index.type.stat") or die $!;
open (OUT1,">$od/$index.Exon.xls") or die $!;
open (OUT2,">$od/$index.Intron.xls") or die $!;
open (OUT3,">$od/$index.Intergenic.xls") or die $!;
foreach my $chr (keys %T) {
	foreach my $s (keys %{$T{$chr}}) {
		next if $T{$chr}{$s}{value}==0;
		my $len=$T{$chr}{$s}{end}-$s+1;
		my $reads=sprintf "%.2f",$T{$chr}{$s}{value}/$len;
		if ($T{$chr}{$s}{type} eq 'Exon') {
			print OUT1 "$reads\t$len\n";
			$e_n+=$T{$chr}{$s}{value};
		}
		elsif ($T{$chr}{$s}{type} eq 'Intron') {
			print OUT2 "$reads\t$len\n";
			$i_n+=$T{$chr}{$s}{value};
		}
		elsif ($T{$chr}{$s}{type} eq 'Intergenic') {
			print OUT3 "$reads\t$len\n";
			$b_n+=$T{$chr}{$s}{value};
		}
	}
}

print OUT "#Type\t$index\n";
print OUT "Exon\t$e_n\n";
print OUT "Intron\t$i_n\n";
print OUT "Intergenic\t$b_n\n";
close OUT;
close OUT1;
close OUT2;
close OUT3;

my $Rscript=$config{Rscript};
my $cmd="$Rscript $Bin/pie.R --infile $od/$index.type.stat --outfile $od/$index.type --bg white ";
print $cmd,"\n";
system"$cmd";
#############################做深度随机性分布图

=pod
my %num;
my %Total_num;
foreach my $gene (keys %Original) {
	my $now_site=0;
	foreach my $s (sort {$a<=>$b} keys %{$Exon{$gene}}) {
		foreach my $ss ($s..$Exon{$gene}{$s}) {
			$now_site++;
			$Original{$gene}{site}{$ss}||=0;
			my $site=(exists $ST{$gene})?(0.01*int(100-100*$now_site/$Len{$gene})):(0.01*int(100*$now_site/$Len{$gene}));
			$Total_num{$site}+=$Original{$gene}{site}{$ss};
			$num{$gene}{site}{$site}+=$Original{$gene}{site}{$ss};
			$num{$gene}{total}+=$Original{$gene}{site}{$ss};
		}
	}
}

chdir $od;
##############################数量
my $max_y=0;
foreach (values %Total_num) {
	if ($_>$max_y) {
		$max_y=$_;
	}
}
my ($y,$ystep)=&define_axis ($max_y);
open (OUT,">$od/$index.randcheck_num.list")||die "$!";
print OUT << "Usage End.";
Type:Line
Width:600
Height:400
WholeScale:0.9
FontSize:25
X:Relative Position in Genes(5'-3')
Y:Number of Depth
XStep:0.1
YStep:$ystep
XStart:0
YStart:0
XEnd:1
YEnd:$y

Color:red
Usage End.
foreach  (sort{$a<=>$b}keys %Total_num) {
	print OUT "$_:$Total_num{$_}\n";
}
close OUT;
#`perl /share/nas2/genome/bmksoft/tool/distributing_svg/v4.75/distributing_svg.pl $od/$index.randcheck_num.list $od/$index.randcheck_num.svg`;
#`/share/nas2/genome/bmksoft/tool/svg2xxx/v1.0/svg2xxx $index.randcheck_num.svg`;


######################百分比

my $gene_num=keys %num;
my %per;
foreach my $key1 (keys %num) {
	if ($num{$key1}{total}<=400) {
		$gene_num--;
		next;
	}
	foreach my $key2 (keys %{$num{$key1}{site}}) {
		$per{$key2}+=$num{$key1}{site}{$key2}/$num{$key1}{total};
	}
}
foreach my $key (keys %per) {
	$per{$key}=sprintf "%.2f",$per{$key}/$gene_num*100;
}
$max_y=0;
foreach (values %per) {
	if ($_>$max_y) {
		$max_y=$_;
	}
}

open (OUT,">$od/$index.randcheck_per.list")||die "$!";
($y,$ystep)=&define_axis ($max_y);
print OUT << "Usage End.";
Type:Line
Width:600
Height:400
WholeScale:0.9
FontSize:25
X:Relative Position in Genes(5'-3')
Y:Percent of Depth
XStep:0.1
YStep:$ystep
XStart:0
YStart:0
XEnd:1
YEnd:$y

Color:red
Usage End.
foreach  (sort{$a<=>$b}keys %per) {
	print OUT "$_:$per{$_}\n";
}
close OUT;
#`perl /share/nas2/genome/bmksoft/tool/distributing_svg/v4.75/distributing_svg.pl $od/$index.randcheck_per.list $od/$index.randcheck_per.svg`;
#`/share/nas2/genome/bmksoft/tool/svg2xxx/v1.0/svg2xxx $index.randcheck_per.svg`;
=cut








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
sub define_axis () {
        my $i=0;
        my ($max)=@_;
        my $time=1;
        my @ret=();
        my @limit=(1,2,3,4,5,6,8,10,12,15,16,20,24,25,30,40,50,60,80,100,120);
        my @unlim=(0,1,2,3,4,5,6,8 ,10,11,14,15,18.5,21,23,29,37,47,56,76 ,92 ,110);

        while ($max >$unlim[21]) {
                 $max=$max/10;
                 $time=$time*10;
        }
        for ($i=0;$i<=20 ;$i++) {
                 if ($max>$unlim[$i] && $max<=$unlim[$i+1]) {
                         $ret[0]=$limit[$i]*$time;

                         if ($i==2 || $i==5 || $i==9 || $i==14) {
                                 $ret[1]=$ret[0]/3;
                         }
                         elsif ($i==4 || $i==7 || $i==13 || $i==16){
                                 $ret[1]=$ret[0]/5;
                         }
                         else {
                                 $ret[1]=$ret[0]/4;
                         }

                 }
        }
        return @ret;
}


sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	Zhang XueChuan <zhangxc\@biomarker.com.cn> 
Usage:
  Options:
  -i     <file>  tophat result *.bam.depth file,forced 
  
  -gff   <file>  gff file,forced 
  
  -index <str>   index of file,forced 
  
  -od    <dir>   output dir,forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
