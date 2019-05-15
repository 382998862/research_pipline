#!usr/bin/perl -w
use strict;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Getopt::Long;

my ($inList,$nameList,$out_dir);
GetOptions(
	"help|?"=>\&USAGE,
	"inList:s" =>\$inList,
	"nameList:s" =>\$nameList,
	"od:s"=>\$out_dir,
)|| die &USAGE;
&USAGE unless ($inList,$nameList,$out_dir);
$file[2]||="None";
$file[3]||="None";
$file[4]||="None";
$out_dir||=".";
mkdir $out_dir unless -d $out_dir;
if (@file >5) {
	&Usage;
}
#my @name=split/\:/,$ID_name;
my %list;my %dir;
my $file_dir = dirname($file[0]);
for (my $i=0;$i <@file ;$i++) {
	if ($file[$i] ne "None") {
		$dir{$name[$i]} = dirname($file[$i]);
		$file[$i]= basename($file[$i]);
#		my ($name)=basename($file[$i])=~/(.*?)\./;
		$list{$name[$i]}=$file[$i];
	}
}
#print Dumper\%list;die;
my @sample_name=sort keys %list;
my @dir = sort keys %dir;
my $sample_num=scalar(@sample_name);
my (%sample_ids,%total_ids);
for(my $i=0;$i<$sample_num;$i++){
	%{$sample_ids{$i}}=&read_list("$dir{$dir[$i]}/$list{$sample_name[$i]}",\%total_ids);
}
#print Dumper\%sample_ids;die;
my ($out_list,$venn_file);
$out_list="$out_dir/list.txt";
open (OUT,">",$out_list)||die "Check your output file!\n";
my $list_header="#ID\t".(join "\t",@sample_name);
print OUT "$list_header\n";

my %com_num;
foreach my $id(keys %total_ids){
	#print $id;die;
	my $mark=10**$sample_num;
		for(my $i=0;$i<$sample_num;$i++){
			if(exists $sample_ids{$i}{$id}){
					$mark+=10**$i;
			}
		}
		my @marks=split //,(reverse $mark);
		my (@com_samples,@infos);
		for(my $i=0;$i<$sample_num;$i++){
			if($marks[$i]){
				push @com_samples,$sample_name[$i];
			}
			my $info=($marks[$i])?$id:"NA";       push @infos,$info;
		}
	my $com_samples=join ",",@com_samples;    $com_num{$com_samples}++;
	my $infos=join "\t",@infos;
	my $out_info=join "\t",($id,$infos);      print OUT $out_info."\n";
}
close OUT;

foreach(sort keys %com_num){
	my $com_samples=$_;
	my $num=$com_num{$com_samples};
	my $out_info=join "\t",($com_samples,$num);   print $out_info."\n";
}

my $R_file=($venn_file)?$venn_file:$out_list;
open (IN,$R_file)||die $!;
my $head=<IN>;chomp $head;
#print $head;die;
close IN;
my @heads=split /\t+/,$head;
shift @heads;
#print $heads[0];die;
my @samples=@heads;

## draw venn picture
	chdir($out_dir);
	$R_file=basename($R_file);
	my $library_grid='library(grid)'."\n";
	my $library_VD='library(VennDiagram)'."\n";
	my $read_table="data<-read.table(\"$R_file\",header=F)"."\n";
	$sample_num=scalar(@samples);
	
	#---------------------------add category.names by huangls--------------------------------------#
	my @category_names=&add_com(@samples);
	my $category_name=join(",",@category_names);
	$category_name="category.name=c(".$category_name.")";		
	my@sa=qw(aa bb cc dd ee);
	@samples=@sa[0..($sample_num-1)];
	#-------------------------------------------------------------------------------------------#
	
	
	my @read_data;
	for(my $i=0;$i<$sample_num;$i++)
	{
		my $column=$i+2;
		push @read_data,"temp<-data[,$column]\n$samples[$i]<-temp[!is.na(temp)]";
	}
	my $read_data=(join "\n",@read_data)."\n";
	my @temp;
	foreach(@samples)
	{
		push @temp,"$_=$_";
		#print $_;die;
	}
	my $sample_name=join ",",@temp;
	my @colors=("'cornflowerblue'","'green'","'yellow'","'darkorchid1'","'red'");
	my $colors=join ",",@colors[0..($sample_num-1)];
	my $list="list($sample_name)";
	my $margin='margin=0.25';
	my $height='height=800';
	my $width='width=800';
	my $resolution='resolution=200';
	my $units="units='px'";
	my $lwd='lwd=1';
	my $fill="fill=c($colors)";
	my $cex=($sample_num==5)?'cex=0.4':'cex=0.7';
	my $cat_cex='cat.cex=0.9';
	my $scaled='scaled=0';
	my $cat_pos=($sample_num==2)?'cat.pos=0':'';
	my $out_pic="filename=\"venn.png\"";
	my $cat_dist=($sample_num==3)?'cat.dist=c(0.12,0.12,0.07)':'';
	my $venn_command="venn.diagram($list,$margin,$height,$width,$resolution,$units,$lwd,$fill,$cex,$cat_cex,$scaled,$cat_pos,$category_name,$out_pic,$cat_dist,imagetype=\"png\")"."\n";
	open OUT,">","venn.R";
	print OUT $library_grid;
	print OUT $library_VD;
	print OUT $read_table;
	print OUT $read_data;
	print OUT $venn_command;
	close OUT;
	chdir($out_dir);
	`/share/nas2/genome/biosoft/R/2.15.1/bin/R -f venn.R`;
## sub functions
## read_file
sub read_list
{
	my ($file,$in_hash)=@_;
	my %hash;
	open IN,"$file";
	print $file,"\n";
	while (<IN>)
	{
		chomp;
		next if /^\s*$/ || /^#/;
		my $id=(split /\t+/)[0];
		$hash{$id}++;
		$in_hash->{$id}++;
	}
	close IN;
	return %hash;
}

sub add_com(){
	my@aa=@_;
	for(my$i=0;$i<@aa;$i++){
		$aa[$i]="\"$aa[$i]\"";
	}
	return @aa;
}
## help information
sub USAGE
{
	my $usage = <<"USAGE";
\nDescription: this script is used to process veen statitic && draw venn pictures;\n
Usage:     perl veen2.pl -T1 C1.txt -T2 M1.txt -od veen_dir -n "sampleID1:sampleID2"
Options:
     -T1   <file>  inputlist1 the file record geneid in the first column.
     -T2   <file>  inputlist2 the format is the same as above.
     -T3   <file>  inputlist3 the format is the same as above.
     -T4   <file>  inputlist4 the format is the same as above.
     -T5   <file>  inputlist5 the format is the same as above.
     -ID1   <file>  give a name for input file one,which will be shown in the veen plot.Note that the name should not be numbers.
     -ID2   <file>  give a name for input file two,which will be shown in the veen plot.Note that the name should not be numbers.
     -ID3   <file>  give a name for input file three,which will be shown in the veen plot.Note that the name should not be numbers.
     -ID4   <file>  give a name for input file four,which will be shown in the veen plot.Note that the name should not be numbers.
     -ID5   <file>  give a name for input file five,which will be shown in the veen plot.Note that the name should not be numbers.
      -od   <str>  output directiory.
Attention:(1) parameter T1 and T2 must be given
          (2) input list num shouldn't be more than 5

USAGE
	print $usage;
	exit;
}
