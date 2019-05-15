#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(max);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
#GetOptions
# ------------------------------------------------------------------
my ($chrom,$gen,$dep,$index,$out,$NUM_0,$NUM_1,$NUM_2,$LOG);
GetOptions(
				"help|?" =>\&USAGE,
				"c:s"=>\$chrom,
				"g:s"=>\$gen,
				"d:s"=>\$dep,
				"o:s"=>\$out,
				"i:s"=>\$index,
				"n0:s"=>\$NUM_0,
				"n1:s"=>\$NUM_1,
				"n2:s"=>\$NUM_2,
				) or &USAGE;
&USAGE unless ($chrom and $gen and $dep and $out and $index);
######################################################################################


mkdir "$out" if (!-d $out);
$chrom = ABSOLUTE_DIR ($chrom);
$gen = ABSOLUTE_DIR ($gen);
$dep = ABSOLUTE_DIR ($dep);
$out = ABSOLUTE_DIR ($out);
mkdir "$out/maps";
$NUM_0 ||= 1000000;
$NUM_1 ||= 500;
$NUM_2 ||= 100;

# ------------------------------------------------------------------
# DataStructure of chromosome
# ------------------------------------------------------------------
open IN, "<", "$chrom" or die "Can't open $chrom:$!\n";
my %chrdata;
$/='>';
my $fa;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp=split/\n+/,$_,2;
	my $id=(split/\s+/,$tmp[0])[0];
	$tmp[1]=~s/\s+//g;
	$chrdata{$id}=length($tmp[1]);
}
$/="\n";
close IN;
for( keys %chrdata){
	if ($chrdata{$_} < $NUM_0){
		delete $chrdata{$_};
	}
}


# ------------------------------------------------------------------
# DataStructure of gff for chromosome maps and random map,random map contains scaffold data
# ------------------------------------------------------------------
my %genedata;			#for chromosome maps
my %ref;				#for random map
my $chr    = '';		#染色体名
my $mRNA   = '';		#mRNA名
my $factor = 2000;		#mRNA 分组参数
my $block;				#mRNA 索引号

my $exon_type=0;
my $gene_type=0;
my $mRNA_type=0;
my $index2_type=0;
my $CDS_type=0;
my $index1_type=0;
my %flag;
open IN, "<", "$gen" or die "Can't open $gen:$!\n";
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/) ;
	if ($_=~/exon/ && $index1_type==0) {
		$exon_type=1;
		$index1_type=1;
	}
	elsif ($_=~/CDS/ && $index1_type==0) {
		$CDS_type=1;
		$index1_type=1;
	}
	if ($_=~/gene/ && $index2_type==0) {
		$gene_type=1;
		$index2_type=1;
	}
	elsif ($_=~/mRNA/ && $index2_type==0) {
		$mRNA_type=1;
		$index2_type=1;
	}
	my @tmp=split /\t+/,$_;
	#print Dumper @tmp;
	if ($tmp[2] eq "mRNA") {
		if ($tmp[6] eq "-"){ 
			$tmp[8] =~m /ID(.?)([^;]+)/;
			$flag{$2}=1;
		}
	}

	if ($tmp[2] eq "gene" && exists $chrdata{$tmp[0]} && $gene_type==1) {
		$genedata{$tmp[0]}->[int($tmp[3] * $NUM_1 / $chrdata{$tmp[0]})]++;
	}
	if ($tmp[2] eq "mRNA" && exists $chrdata{$tmp[0]} && $mRNA_type==1) {
		$genedata{$tmp[0]}->[int($tmp[3] * $NUM_1 / $chrdata{$tmp[0]})]++;
	}
	if ($tmp[2] eq "mRNA") {
		$chr = $tmp[0] if ($tmp[0] ne $chr);
		$tmp[8] =~m /ID(.?)([^;]+)/;
		$ref{$chr}->{int($tmp[3]/$factor)}->{$2}->[0] = $tmp[3];
		$ref{$chr}->{int($tmp[3]/$factor)}->{$2}->[1] = $tmp[4];
		$mRNA = $2;
		$block = int($tmp[3]/$factor);
	}
	if ($tmp[2] eq "exon" && $exon_type==1) {
		push @{$ref{$chr}->{$block}->{$mRNA}->[2]}, [$tmp[3],$tmp[4],defined $ref{$chr}->{$block}->{$mRNA}->[3] ? $ref{$chr}->{$block}->{$mRNA}->[3] : 0];
		$ref{$chr}->{$block}->{$mRNA}->[3]+=($tmp[4]-$tmp[3]+1);
	}
	if ($tmp[2] eq "CDS" && $CDS_type==1) {
		push @{$ref{$chr}->{$block}->{$mRNA}->[2]}, [$tmp[3],$tmp[4],defined $ref{$chr}->{$block}->{$mRNA}->[3] ? $ref{$chr}->{$block}->{$mRNA}->[3] : 0];
		$ref{$chr}->{$block}->{$mRNA}->[3]+=($tmp[4]-$tmp[3]+1);
	}
}
close IN;
#print Dumper %ref;die;

# ------------------------------------------------------------------
# DataStructure of depth for chromosome maps
# ------------------------------------------------------------------
my @depth;
my %depdata;
open IN, "<", "$dep" or die "Can't open $dep:$!\n";
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my @tmp=split/\t+/,$_;
	if (exists $chrdata{$tmp[0]}) {
		my $sub = int ($tmp[1] * $NUM_1 / $chrdata{$tmp[0]}); #subscript is autovivification
		$depdata{$tmp[0]}->[$sub][0] += $tmp[2];
		$depdata{$tmp[0]}->[$sub][1] ++ if ($tmp[2] != 0);
	}
	if ($tmp[2] != 0 && exists $ref{$tmp[0]}->{int($tmp[1]/$factor)}) {
		for my $m (keys %{$ref{$tmp[0]}->{int($tmp[1]/$factor)}}) {
			if ($tmp[1] >= $ref{$tmp[0]}->{int($tmp[1]/$factor)}->{$m}->[0] && $tmp[1] <= $ref{$tmp[0]}->{int($tmp[1]/$factor)}->{$m}->[1]) {
				for (@{$ref{$tmp[0]}->{int($tmp[1]/$factor)}->{$m}->[2]}) {
					if ($tmp[1] >= $_->[0] && $tmp[1] <= $_->[1]) {
						if (defined $flag{$m}) {
							$depth[int($NUM_2 * ($ref{$tmp[0]}->{int($tmp[1]/$factor)}->{$m}->[3]-($tmp[1] - $_->[0] + $_->[2])-1)/$ref{$tmp[0]}->{int($tmp[1]/$factor)}->{$m}->[3])] += $tmp[2];
						}else{
							$depth[int($NUM_2 * ($tmp[1] - $_->[0] + $_->[2])/$ref{$tmp[0]}->{int($tmp[1]/$factor)}->{$m}->[3])] += $tmp[2];
						}
					}
				}
			}
		}
	}
}
close IN;
#print Dumper @depth;die;
# ------------------------------------------------------------------
# SVG
# ------------------------------------------------------------------

my @chrs = sort {$a cmp $b} keys %depdata;
#chr map

for my $chr (@chrs) {
	open OUT, ">", "$out/maps/${chr}.svg" or die "Can't creat $out/maps/${chr}.svg:$!\n";
	print OUT &svg_paper('650','450');

	print OUT &svg_line('100','120','600','120','#EE0000');
	print OUT &svg_line('100','120','100','20','#EE0000');
	print OUT &svg_line('100','260','600','260','#00EE00');
	print OUT &svg_line('100','260','100','160','#00EE00');
	print OUT &svg_line('100','400','600','400','#0000EE');
	print OUT &svg_line('100','400','100','300','#0000EE');

	print OUT &svg_txt('75','110','20','#000000','log2(depth)',3);
	print OUT &svg_txt('75','250','20','#000000','Coverage',3);
	print OUT &svg_txt('75','390','20','#000000','GeneNumber',3);

	my (@Y1,@Y2,@Y3);
	for (@{$depdata{$chr}}) {
		if (defined $_) {
			push @Y1,$_->[0];
			if (defined $_->[1]) {
				push @Y2,$_->[1];
			}
		}
	}
	for (@{$genedata{$chr}}) {
		if (defined $_) {
			push @Y3,$_;
		}
	}
	my $Y1_max = max @Y1;
	my $Y2_max = max @Y2;
	my $Y3_max = max @Y3;

	my $Y1_top = int(log($Y1_max)/log(2)) + 1;
	my $Y2_top = (int(10 * $Y2_max * 500 / $chrdata{$chr} ) + 1 ) / 10;
	my $Y3_top = $Y3_max; 
	
	print OUT &svg_txt('90','130','20','#000000','0');
	print OUT &svg_txt('80','30','20','#000000',"$Y1_top");
	print OUT &svg_txt('90','270','20','#000000','0');
	print OUT &svg_txt('80','170','20','#000000',"$Y2_top");
	print OUT &svg_txt('90','410','20','#000000','0');
	print OUT &svg_txt('80','310','20','#000000',"$Y3_top");
	my $region = int($chrdata{$chr}/500);
	print OUT &svg_txt('100','430','25','#000000',"chr $chr (${region}nt/windows,500windows)");
	my ($Y1_last,$Y2_last,$Y3_last) = (0,0,0);
	my $X = 0;

	for (@{$depdata{$chr}}) {
		my ($Y1,$Y2) = defined($_->[1]) ? ($_->[0],$_->[1]) : (0,0);
		$Y1 = ($Y1 == 0) ? 0 : (100 * log($Y1)/log(2) / $Y1_top);
		$Y2 = 100 * ($Y2 * 500 / $chrdata{$chr} ) / $Y2_top;
		print OUT &svg_line(100+$X,120-$Y1_last,101+$X,120-$Y1,"#EE0000");
		print OUT &svg_line(100+$X,260-$Y2_last,101+$X,260-$Y2,"#00EE00");
		$X++;
		($Y1_last,$Y2_last) = ($Y1,$Y2);
	}
	#Need another loop since the keys's number are not equal
	$X = 0;

	for (@{$genedata{$chr}}) {
		my $Y3 = defined $_ ? $_ : 0;
		$Y3 = 100 * $Y3 / $Y3_top;
		print OUT &svg_line(100+$X,400-$Y3_last,101+$X,400-$Y3,"#0000EE");
		$X++;
		$Y3_last = $Y3;
	}
	print OUT &svg_end();
	close OUT;
}

#depth map
my @x = ('0.0','0.2','0.4','0.6','0.8','1.0');
my $y_max = max @depth;
open OUT, ">", "$out/maps/$index.randCheck.svg" or die "Can't creat $out/maps/$index.randCheck.svg:$!\n";
open OUT1, ">", "$out/maps/$index.randCheck" or die "Can't creat $out/maps/$index.randCheck:$!\n";
print OUT &svg_paper(650,450);
print OUT &svg_line(100,370,600,370,'black');
print OUT &svg_line(80,350,80,50,'black');
my $wide_step = 500/$NUM_2;

for ( @x ) {
	print OUT &svg_line(100+500*$_,370,100+500*$_,378,'black');
	print OUT &svg_line(80,350-300*$_,72,350-300*$_,'black');
	print OUT &svg_txt(90+500*$_,400,12,'black',$_);
	print OUT &svg_txt(60,375-300*$_,12,'black',sprintf ("%.1e",$y_max * $_),3);
}

for (0..$#depth) {
	my @rec = (100 + $wide_step * $_, 350 - 300 * $depth[$_] / $y_max, $wide_step, 300 * $depth[$_] / $y_max,'red','black', 1);
	print OUT1 "$_:$depth[$_]\n";
	print OUT "<rect x=\"$rec[0]\" y=\"$rec[1]\" width=\"$rec[2]\" height=\"$rec[3]\" fill=\"$rec[4]\" stroke=\"$rec[5]\" stroke-width=\"$rec[6]\"/>\n";
}

print OUT &svg_txt (260,420,18,'black',"Relative Position");
print OUT &svg_txt (20,250,18,'black',"Total Depth",3);
print OUT &svg_end();
close OUT;
close OUT1;

chdir "$out/maps";
my $svg2png = "/share/nas2/genome/biosoft/distributing_svg_4.74/svg2xxx_release/svg2xxx";
for my $svg_map (glob "*.svg"){
	`$svg2png $svg_map`;
}


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################



# ------------------------------------------------------------------

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

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub svg_txt (){#&svg_txt(x,y,size,color,text,[vertical,0/1/2/3]);
	my @svg_x=@_;
	if (!defined $svg_x[5]) {
		$svg_x[5]=0;
	}
	my $svg_matrix='';
	if ($svg_x[5]==0) {
		$svg_matrix="1 0 0 1";
	}
	if ($svg_x[5]==1) {
		$svg_matrix="0 1 -1 0";
	}
	if ($svg_x[5]==2) {
		$svg_matrix="-1 0 0 -1";
	}
	if ($svg_x[5]==3) {
		$svg_matrix="0 -1 1 0";
	}
	my $line="<text fill=\"$svg_x[3]\" transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" font-family=\"ArialNarrow-Bold\" font-size=\"$svg_x[2]\">$svg_x[4]</text>\n";
	return $line;
}

sub svg_line (){#&svg_line(x1,y1,x2,y2,color,[width])
	my @svg_x=@_;
	my $line="<line fill=\"$svg_x[4]\" stroke=\"$svg_x[4]\" x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\"/>\n";
	if (defined $svg_x[5]) {
		$line="<line fill=\"$svg_x[4]\" stroke=\"$svg_x[4]\" stroke-width=\"$svg_x[5]\" x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\"/>\n";
	}
	return $line;
}

sub svg_paper (){#&svg_paper(width,height,[color])
	my $svg_drawer = (getlogin() || getpwuid($<))."@".(`hostname`);
	chomp $svg_drawer;
	my @svg_x=@_;
	my $line="";
	$line.="<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n";
	$line.="<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 20001102//EN\" \"http://www.w3.org/TR/2000/CR-SVG-20001102/DTD/svg-20001102.dtd\">\n\n";
	$line.="<svg width=\"$svg_x[0]\" height=\"$svg_x[1]\" viewBox=\"0 0 $svg_x[0] $svg_x[1]\">\n";
	$line.="<Drawer>$svg_drawer</Drawer>\n";
	$line.="<Date>".(localtime())."</Date>\n";
	if (defined $svg_x[2]) {
		$line.="<rect x=\"0\" y=\"0\" width=\"$svg_x[0]\" height=\"$svg_x[1]\" fill=\"$svg_x[2]\"/>\n";
	}
	return $line;
}

sub svg_end (){#
	return "</svg>\n";
}

# ------------------------------------------------------------------

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Zhang Yuheng <zhangyh\@biomarker.com.cn>

  It is used to Creat depth maps and random map

  example
  -----------------------------------------------------------------------------
  perl $0 -c Gmax1.01.RM.N.fasta -g Glyma1_highConfidence.gff -d YZF1.genome.bam.depth -o gffmaps
  -----------------------------------------------------------------------------
  -c         chromosome
  -g         gff
  -d         depth
  -o         outdir
  -i         index of out file

  -n0        default is 1000000      [optional]
  -n1        default is 500          [optional]
  -n2        default is 100           [optional]
  -h         Help

USAGE
	print $usage;
	exit;
}
