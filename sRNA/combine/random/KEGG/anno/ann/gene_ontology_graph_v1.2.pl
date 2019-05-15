#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd qw(realpath);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my %config=%{readconf("$Bin/../database.config")}; 

my $BEGIN_TIME=time();
my $version="1.0.0";
my $author="Sunhy";
my $date="2013-01-15";
my @Times = localtime();
my $year=$Times[5]+1900;
my $month=$Times[4]+1;
my $day=$Times[3];
#######################################################################################
my $Time_Start = sub_format_datetime(localtime(time()));
print STDOUT "$Script start at:[$Time_Start]\n";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
#---------------------------------options----------------------------------------------
my ($component_file,$function_file,$process_file,$foutdir,$help,@opt_i,@opt_mark,$fkey);
GetOptions(
			"h|help"=>\$help,
			"i:s"=>\@opt_i,
			"mark:s"=>\@opt_mark,
			"k:s"=>\$fkey,
			"c:s"=>\$component_file,
			"f:s"=>\$function_file,
			"p:s"=>\$process_file,
			"od:s"=>\$foutdir,
);

#$component_file   ||=  "/share/nas2/database/GO/releases20120916/GO_Classify/cellular_component.txt"; # 2012-09-16 ~ 2014-12-12
#$function_file    ||=  "/share/nas2/database/GO/releases20120916/GO_Classify/molecular_function.txt"; # 2012-09-16 ~ 2014-12-12
#$process_file     ||=  "/share/nas2/database/GO/releases20120916/GO_Classify/biological_process.txt"; # 2012-09-16 ~ 2014-12-12
#$component_file   ||=  "/share/nas2/database/GO/releases20141211/GO_Classify/cellular_component.txt"; # 2014-12-12 ~ 
#$function_file    ||=  "/share/nas2/database/GO/releases20141211/GO_Classify/molecular_function.txt"; # 2014-12-12 ~ 
#$process_file     ||=  "/share/nas2/database/GO/releases20141211/GO_Classify/biological_process.txt"; # 2014-12-12 ~ 
$foutdir          ||=  './';
$help             ||=  "";
#-------------------------------------------------------------------------------------
&help () if($help);
&usage() if(@opt_i <=0 || !defined $foutdir || @opt_i != @opt_mark || !defined $fkey || !defined $component_file ||!defined $function_file ||!defined $process_file );

foreach (@opt_mark){
	s/\_/ /g;
	$_ = ucfirst $_;	
}
$foutdir = Cwd::realpath($foutdir);
mkdir("$foutdir") if (!-d "$foutdir");
my $nowdir = Cwd::realpath();
die"\n\tERROR:The number of input file should not bigger than 3!" if (@opt_i>4);
# --------------------------------------------------------------------------------------


my $dissvg_pl=$config{distributing_svg2};
my $svgxxx=$config{svg2xxx2};

my $note_type = "GO Classification";
my $graph_move=0.1;
my @color_set=("#800000","#191970","#556b2f","#F96611");
my @rect_width=("0.6","0.35","0.25","0.2");
my @moveper=("0.2","0.15","0.125","0.1");
my @three_go_tree_file=();
my @percentage_out=();
my @tree_root = ();
my @real=();


my %percentage_out=();
my %three_go_tree_key=();
my %three_go_tree_hash=();
my %file_provide=();
my %group_count=();
my %correlation=();
my %out_hash=();
my %out_key=();
my %clean_key=();
my %total_gene;
@three_go_tree_file = ($component_file, $function_file, $process_file);
&get_GO_tree_hash(\@three_go_tree_file,\%three_go_tree_hash,\%three_go_tree_key,\@tree_root);

&get_file_info(\@opt_i, \@opt_mark, \%file_provide,\@real);


#------------------find correlation--------------------------------------------
foreach my $k (keys %three_go_tree_hash) {
	foreach my $k1 (keys %{$three_go_tree_hash{$k}}) {
		foreach my $k2 (keys %file_provide) {
			if (defined $file_provide{$k2}{$k1}) {
				foreach my $k3 (keys %{$file_provide{$k2}{$k1}}) {
					$correlation{$k}{$k2}{$k3}++;
				}
			}
		}
	}
}

foreach my $k (keys %correlation) {
	foreach my $k1 (keys %{$correlation{$k}}) {
		my $num=0;
		$num = (keys %{$correlation{$k}{$k1}});
		$out_hash{$k}{$k1} = $num;
	}
}

#----------------------------------clean key-------------------------------------------
foreach my $k (@tree_root) {
	my @new_cla=();
	foreach my $k1 (@{$three_go_tree_key{$k}}) {
		if (defined $out_hash{$k1}) {
			push @new_cla,$k1;
		}
	}
	$clean_key{$k} = \@new_cla;
}
#----------------------------------output statistics------------------------------------
open STAT,">$foutdir/$fkey.GO_enrichment.stat.xls" || die;
print STAT "#GO_classify1\tGO_classify2\t",join("\t",@opt_mark),"\n";
print STAT "#Total_gene\t";
foreach(@opt_mark){
	print STAT "\t$total_gene{$_}";
}
print STAT "\n";

foreach my $k (@tree_root) {
	foreach my $k1 (@{$clean_key{$k}}) {
		print STAT "$k\t$k1";
		foreach my $k2 (@opt_mark) {
			$out_hash{$k1}{$k2} ||= 0;
			print STAT "\t$out_hash{$k1}{$k2}";
		}
		print STAT "\n";
	}

}
close STAT;


my $cmd="$config{Rscript} $Bin/GOClassificationMap2.r --infile $foutdir/$fkey.GO_enrichment.stat.xls --outpath $foutdir --key $fkey ";
print $cmd,"\n";
system $cmd;

#---------------------------------graph--------------------------------------------------#
my $Time_End = sub_format_datetime(localtime(time()));
print STDOUT "\n$Script Done at: [$Time_End]\t\tTotal elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub get_GO_tree_hash ($$$$){  ## from three go tree file to a tree hash and a array
	my ($go_file,$go_tree_hash,$go_tree_key,$root) = @_;
	my $key;
	foreach my $file (@$go_file) {
		open IN,$file or die "cannot open file $file, $!\n";
		my @array = ();
		my $name;
		my $key;
		while (my $line = <IN>) {
			chomp($line);
			next if ($line =~ /^$|^\#/);
			next if ($line !~ /(^\s+>)|(^\s+%)/);
			if ($line =~ /^\s>(.*)\s;/) {
				$name = $1;
				$name =~ s/[_\-]/ /;
				next;
			}
			if ($line=~/^\s\s%/) {
				my @key_in=split(/\s+;\s+/,$line);
				$key_in[0]=~s/^\s\s%//;
				$key = $key_in[0];
				push @array, $key;
				my ($go_id) = $key_in[1] =~ /(GO:\d+)/;
				$$go_tree_hash{$key}{$go_id}++;
			}else{
				my ($go_id) = $line =~ /(GO:\d+)/;
				next if (!defined $go_id);
				$$go_tree_hash{$key}{$go_id}++;
			}
		}
		$$go_tree_key{$name} = \@array;
		push @$root,$name;
		close IN;
	}
}
sub get_file_info ($$$$) {
	my ($files, $marks, $hash, $re) = @_;
	for (my $i = 0;$i < @$files ;$i++) {
		open IN,"$$files[$i]" or die "cannot open file $$files[$i], $!\n";
		my $line_num=0;
		while (<IN>) {
			chomp;
			next if (/^$/||/\#/);
			s/\t+$//;
			my @lines = split /\t/, $_;
			next if (@lines < 2);
			my $gene_id = shift @lines;
			$line_num++;
			foreach (@lines) {
				$$hash{$$marks[$i]}{$_}{$gene_id}++;
			}
			
		}
		push @$re, $line_num;
		$total_gene{$$marks[$i]}=$line_num;
		close IN;
	}
}
sub tuning_svg {
	my $raw=shift;
	open (IN,$raw) ||die;
	open (OUT,">$raw.raw")||die;
	while(<IN>)
	{
		chomp;
		if(/^\<svg width=\"([.\d]+)\"/){
			s/($1)/$1+150/e;
		}
		print OUT "$_\n";
	}
	close IN;
	close OUT;
}

sub graph_out {
	my ($x_end,$width,$moveper)=@_;
	$moveper+=$graph_move;
	my $dis =<<"OUT";
Type:Simple
Width:3000
Height:800
BothYAxis:1
MultiRY:1
#MultiY:1
ScaleLen:8
WholeScale:1
OffsetPer:$width
UnitPer:$width
MovePer:$moveper
XStep:1
YStep:33.3
RYStep:6
XScalePos:0.5
XScaleRoate:80
XStart:0
YStart:0.1
RYStart:0.1
XEnd:$x_end
YEnd:100
RYEnd:100
YNeedLog:10
MarkNoBorder:1
WholeScale:0.95
Note:$note_type 
X:
Y:Percent of genes
RY:Number of genes
XUnit:1
Scale:
OUT

return $dis;
}

sub graph_part_out ($$$$) {
	my ($x_end,$width,$moveper,$part_note)=@_;
	$moveper+=$graph_move;
	my $part_width = 75 * ($x_end - 1);
	my $dis =<<"OUT";	
Type:Simple
Width:$part_width
Height:800
BothYAxis:1
MultiRY:1
#MultiY:1
ScaleLen:8
WholeScale:1
OffsetPer:$width
UnitPer:$width
MovePer:$moveper
XStep:1
YStep:33.3
RYStep:6
XScalePos:0.5
XScaleRoate:80
XStart:0
YStart:0.1
RYStart:0.1
XEnd:$x_end
YEnd:100
RYEnd:100
YNeedLog:10
MarkNoBorder:1
WholeScale:0.95
Note: 
X:
Y:Percent of genes
RY:Number of genes
XUnit:1
Scale:
OUT

return $dis;
}

sub get_graph_list ($$$$){
	my ($r_percentage,$root,$tree_key,$mark)=@_;
	my $i=0;
	my $list_line;
	foreach my $k (@$root) {
		foreach my $k1 (@{$$tree_key{$k}}) {
			my $nu = $$r_percentage{$k1}{$mark};
			$nu = 0.1 if ($nu < 0.1);
			$list_line .= "$i:$nu\n";
			$i++;
		}
	}
	return $list_line;
}

sub get_part_graph_list ($$$$$) {
	my ($r_percentage,$tree_key,$mark,$part)=@_;
	my $i=0;
	my $list_line;
	foreach my $k1 (@{$$tree_key{$part}}) {
		my $nu = $$r_percentage{$k1}{$mark};
		$nu = 0.1 if ($nu < 0.1);
		$list_line .= "$i:$nu\n";
		$i++;
	}
	return $list_line;
}

sub get_percentage ($$$$$$) {  # get ordinate
	my ($outhash, $mark,$root,$tree_key, $real ,$percent) = @_;
	foreach my $k (@$root) {
		foreach my $k1 (@{$$tree_key{$k}}) {
			for(my $i = 0;$i< @$mark;$i++) {
				$$percent{$k1}{$$mark[$i]} = sprintf "%.4f",($$outhash{$k1}{$$mark[$i]}/$$real[$i]*100000+0.5)/100000*100;
			}
		}
	}
}

sub sub_format_datetime {   #Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub usage {
    print STDERR <<"    _EOT_";
	 +---------------------------------------------------------------------+
	 | According to the GO id that relate with Gene, we classify the genes |
	 | to GO Classification. and generate a table for output. then get a   |
	 | intuitionistic SVG&PNG plot.                                        |
	 +---------------------------------------------------------------------+
    For more help type -h or -help;

    
    Usage: program <options> <specification file> [-i ...[-i ...]]
    
-c        Component file name                 "component_ontology.txt"   Optional
           
     -f        Function file name                  "funtion_ontology.txt"     Optional
           
     -p        Process file name                   "process_ontology.txt"     Optional
                       
     -i        Input nm go file name               "*.list.txt"               Required
           
     -mark     Mark to distinguish multi-column    "RICE" e.g                 Required
           
     -od       Output dir name                     default[./]                Optional
      
	 -k        the prefix                          "Pumpkin.Unigene.fa.GO"    Required
           
     -h       Help                                
           
    _EOT_
    exit(1);
}

sub help () {
	print STDERR <<"    _EOT_";
	 +---------------------------------------------------------------------+
	 | According to the GO id that relate with Gene, we classify the genes |
	 | to GO Classification. and generate a table for output. then get a   |
	 | intuitionistic SVG&PNG plot.                                        |
	 +---------------------------------------------------------------------+
    
     -c        Component file name                 "component_ontology.txt"   Optional
           
     -f        Function file name                  "funtion_ontology.txt"     Optional
           
     -p        Process file name                   "process_ontology.txt"     Optional
                       
     -i        Input nm go file name               "*.list.txt"               Required
           
     -mark     Mark to distinguish multi-column    "RICE" e.g                 Required
           
     -od       Output dir name                     default[./]                Optional
      
	 -k        the prefix                          "Pumpkin.Unigene.fa.GO"    Required
           
           
     -h       Help                                
           
    -c, -f, -p the three files are the three Gene Ontology files.
        
        [The three files maybe in TXT format,not HTML,please do this 
         translate after you get the three GO tree file from the web.]
         
    
    -i is the Gene with GO number you want to classify.
    
        [Gene number followed with GO numbers, they must split by 
         <TABLE>, the GO number can be split by <TABLE> .
         You are allowed to give multiple input files.]
         
         
    -mark is to distinguish the data in the plots in different color.
        
        [the number of marks must equal that of input files. 
			and the marks must be different]
        
    -o is the output file name you want.
    
    	[It must be with ".svg" as extend name.]
    
        
    REMEMBER to update your GO tree when new version produced!

     
    _EOT_
    exit;
}
