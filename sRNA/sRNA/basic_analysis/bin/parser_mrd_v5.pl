#!usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
##***********
#This script is used to parser .mrd file.
#4 files will be generated.
#1.mature.mrd
#2.novel_conservative.mrd
#3.novel_unconservative.mrd
#4.miRNA_expression.list
#************   
my($mrd_file,$prefix,$out_dir,$fa,$Tscore,$Tmfe,$Tstar,$Mcount,$coord_file,$str_file)	;
GetOptions(
				"help|?" 	=>\&USAGE,
				"mrd:s"		=>\$mrd_file,
				"samples:s"	=>\$prefix,
				"od:s"		=>\$out_dir,
				"known:s"	=>\$fa,
				"Tscore:i"	=>\$Tscore,
				"Tmfe:i"	=>\$Tmfe,
				"Tstar:i"	=>\$Tstar,
				"Mcount:i"	=>\$Mcount,
				"coord:s"	=>\$coord_file,
				"str:s"		=>\$str_file,
				) or &USAGE;
&USAGE unless ($mrd_file and $prefix and $out_dir and $coord_file and $str_file);
#&USAGE unless ($mrd_file and $prefix and $coord_file and $str_file);
#print "$mrd_file\n";
#print "$prefix\n";
#print "$out_dir\n";
#print "$coord_file\n";
#print "$str_file\n";

$Tscore||=0 ;
$Tmfe||=0.4;
$Tstar||=-2;
$Mcount||=5;
$mrd_file=&ABSOLUTE_DIR($mrd_file);
$coord_file=&ABSOLUTE_DIR($coord_file);
mkdir $out_dir unless -d $out_dir;
$out_dir=&ABSOLUTE_DIR($out_dir);

my @prefixs=split /,/,$prefix;

my %FA;
if(defined $fa)
{
	extract_fa($fa,%FA);
}

my @header;
push @header,"<tag id>";
push @header,"<mature MiRNA id>";
push @header,"<potential ids>";
push @header,"<conserved miRNAs>";
push @header,"<mature length>";#改名字
push @header,"<mature sequence>";#改名字
push @header,"<precursor sequence>";#增加一列
push @header,"<total score >";
push @header,"<score for star>";
push @header,"<mfe>";
push @header,"<mismatch>";
push @header,"<total read count>";
push @header,"<mature read count>";

push @header,"<star sequence>";
push @header,"<strand>";
push @header,"<start>";
push @header,"<end>";
push @header,"<hairpin structure>";
push @header,"<hairpin energy>";
#foreach(sort @prefixs)
#{
#push @header,"<mature read count:$_>";
#}

open(IN,"$coord_file")||die "check whether $coord_file exist!\n";
my (%strand,%start,%end);
while(<IN>){
	my @line=split(/\s+/,$_);
	$line[0]=~s/>//g;
	$strand{$line[0]}=$line[1];
	$start{$line[0]}=$line[2];
	$end{$line[0]}=$line[3];
}
close IN;

open(IN,"$str_file")||die "check whether $str_file exist!\n";
$/=">";<IN>; my %energy;
while(<IN>){
	chomp;
	my @line=split(/\n/,$_);
	#my $len=length($line[2]);
	if($line[2]=~/.+\((.+)\)/){
		#print "$1\n";
        $energy{$line[0]}=$1;
	}
}
close IN;
$/="\n";

parser_mrd($mrd_file,$prefix,$out_dir);


my $novel_conserve="$out_dir/novel_conservative.mrd";
my $novel_conserve_len=`less -S $novel_conserve|wc -l`;
chomp $novel_conserve_len; 

my %CONSERV;
if($novel_conserve_len!=0){
    open (IN,$novel_conserve) or die $!;
	$/=">";
	<IN>;
	while(<IN>){
		chomp;
		my ($id,$seed)=(split /\n/,$_)[0,7];
		$seed=(split /\t/,$seed)[1];
		$CONSERV{$id}=$seed;
	}
	close IN;
	$/="\n";
	
	open (OUT, ">$out_dir/novel_conservative.species") or die $!;
	foreach my $name(keys %CONSERV){
		print OUT "conservative_$name\t$CONSERV{$name}\n";
	}
	close OUT;
}


sub parser_mrd
{
	my ($file,$pre,$dir)=@_;

	chdir($dir);
	
	open (MRD,$file) || die "Check your file $file whether exsisted!\n";
	$/=">";
	open MRD1,">","mature.mrd";
	open MRD2,">","novel_conservative.mrd";
	open MRD3,">","novel_unconservative.mrd";
	open OUT,">","miRNA_expression.list";
	
	print OUT "#";
	print OUT (join "\t",@header);
	print OUT "\n";

	while(<MRD>)
	{
		chomp;
		my $origion_info=$_;
		next if /^\s*$/;

		my ($chr)         =$_=~/^(\S+)/;
		my ($score)       =$_=~/score\s+total\s+(\S+)/;
		my ($total_count) =$_=~/total\s+read\s+count\s+(\S+)/;
		my ($read_count)  =$_=~/\nmature\s+read\s+count\s+(\S+)/;
		my ($score_star)  =$_=~/\nscore\s+for\s+star\s+read\(s\)+\s+(\S+)/;
		my ($mfe)         =$_=~/\nscore\s+for\s+mfe\s+(\S+)/;
		my ($pri_seq)     =$_=~/\npri_seq\s+(\S+)/;
		my ($obs_seq)     =$_=~/\nobs\s+(\S+)/;
		my ($exp_seq)     =$_=~/\nexp\s+(\S+)/;
		my ($conserve_miR)=$_=~/\nmiRNA\s+with\s+same\s+seed\s*(\S+)\n/;
		my ($pri_struct)    =$_=~/\npri_struct\s+(\S+)/;
		$conserve_miR||="--";
		
		next if $score<$Tscore;
		#next if $mfe<$Tmfe;
		#next if $score_star<$Tstar;
		next if $read_count<$Mcount;
		
		my @lines=split /\n+/,$_;
		#print "@lines\n";
		my %muti_count;
		foreach(@lines)
		{
			next unless /^(\S{3})\s+mature\s+read\s+count\s+(\d+)/;
			print "line:$_\n";
			my $p_sam=$1;
			my $p_count=$2;
			$muti_count{$p_sam}=$p_count;
		}
		my $muti_count;
		foreach my $k (sort @prefixs)
		{
			$muti_count{$k}||=0;
			$muti_count.="$muti_count{$k}\t";
		}
		$muti_count=~s/\t$//;
		print "$muti_count\n";


		while()	
		{
			last unless @lines;
			unless($lines[0]=~/\s+#MM/){
				shift @lines;   next;
			}
			shift @lines;   last;
		}
		#print "@lines\n";
		my @lines_=@lines;
		foreach(@prefixs)
		{
			my $p_=$_;
			@lines_=grep (!/^$p_/,@lines_);
			#print "round:@lines_\n";
		}
		#print "@lines_\n";
		if(scalar(@lines_)==0){## Consider the novel miRNAs
			my $map_seq=($obs_seq?$obs_seq:$exp_seq);
			my $M_num  =($map_seq=~s/M/M/g);
			my $MMMMM  ="M"x$M_num;
			my $index  =index($map_seq,$MMMMM);
			my $seed   =substr($pri_seq,$index,$M_num);
			
			my $S_num  =($map_seq=~s/S/S/g);
			my $SSSSS  ="S"x$S_num;
			   $index  =index($map_seq,$SSSSS);
			my $star_seq=substr($pri_seq,$index,$S_num);
			
			print OUT "$chr\t";
			print OUT "$chr\t";
			print OUT "novel\t";
			print OUT "$conserve_miR\t";
			print OUT length($seed)."\t";
			print OUT "$seed\t";
			print OUT "$pri_seq\t";
			print OUT "$score\t";
			print OUT "$score_star\t";
			print OUT "$mfe\t";
			print OUT "0\t";
			print OUT "$total_count\t";
			print OUT "$read_count\t";
#			print OUT "$muti_count\n";
			
			
			print OUT "$star_seq\t";
			print OUT "$strand{$chr}\t";
			print OUT "$start{$chr}\t";
			print OUT "$end{$chr}\t";
			print OUT "$pri_struct\t";
			print OUT "$energy{$chr}\n";
			
			
			if($conserve_miR ne '--')
			{
				print MRD2 ">$origion_info";
			}else
			{
				print MRD3 ">$origion_info";
			}

			next;

		}
		print MRD1 ">$origion_info";
		
		## Consider known miRNAs
		#print "Origin:start\n";
		#print "@lines_\n";
		#print "Origin:end\n";
		@lines_    =map {[split /\s+/]} @lines_;
=c		
		foreach (@lines_){
			print "start:\n";
			print "@{$_}" ;
			print "\nend\n";
		}	
=cut		
		@lines_    =sort {$a->[-1]<=>$b->[-1]} @lines_;
=c
		foreach (@lines_){
			print "start:\n";
			print "@{$_}" ;
			print "\nend\n";
		}		
=cut		
		#print "$lines_[0][-1]";
		my $min_mis=$lines_[0]->[-1];
		#print "$min_mis\n";
		@lines_    =map {useful($_,$min_mis)} @lines_;
		@lines_    =grep (/[\S+]/,@lines_);
=c		
		foreach (@lines_){
			print "start:\n";
			print "@{$_}" ;
			print "\nend\n";
		}		
=cut		
		foreach (@lines_)
		{
			$_->[1]=~s/\.//g;
			$_->[1]=lc($_->[1]);
		}
		@lines_=sort{length($a->[1])<=>length($b->[1])} @lines_;
		my @id =map{$_->[0]} @lines_;
		print "@id\n";
		my ($id_use,$seq_use,$mis_use)=@{$lines_[0]};
		my $ids=join ",",@id;
		my $mature_seq=(exists $FA{$ids})?$FA{$ids}:$seq_use;
		print OUT "$chr\t";
		print OUT "$id_use\t";
		print OUT "$ids\t";
		print OUT "$conserve_miR\t";
		print OUT length($mature_seq)."\t";
		print OUT "$mature_seq\t";
		print OUT "$pri_seq\t";
		print OUT "$score\t";
		print OUT "$score_star\t";
		print OUT "$mfe\t";
		print OUT "$mis_use\t";
		print OUT "$total_count\t";
		print OUT "$read_count\t";
#		print OUT "$muti_count\n";
		
		my $map_seq=($obs_seq?$obs_seq:$exp_seq);
		my $M_num  =($map_seq=~s/M/M/g);
		my $MMMMM  ="M"x$M_num;
		my $index  =index($map_seq,$MMMMM);
		my $seed   =substr($pri_seq,$index,$M_num);
			
		my $S_num  =($map_seq=~s/S/S/g);
		my $SSSSS  ="S"x$S_num;
		$index  =index($map_seq,$SSSSS);
		my $star_seq=substr($pri_seq,$index,$S_num);
			
		#my $map_seq=($obs_seq?$obs_seq:$exp_seq);
		$index  =index($map_seq,$seq_use);
		my $seq   =substr($pri_seq,$index,20);
		my $star_seq_m="";
		if($seq=~/S/){
			$star_seq_m=$seed;
		}
		else{
			$star_seq_m=$star_seq;
		}
		
		print OUT "$star_seq_m\t";
		print OUT "$strand{$chr}\t";
		print OUT "$start{$chr}\t";
		print OUT "$end{$chr}\t";
		print OUT "$pri_struct\t";
		print OUT "$energy{$chr}\n";

	}
	close MRD3;
	close MRD2;
	close MRD1;
	$/="\n";
	close MRD;

}

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
		print "$in\n";
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub extract_fa
{
	my ($fa_file,%fa)=@_;
	open (IN,$fa_file) or die $!;  
	$/=">";
	while(<IN>)
	{
		chomp;          
		next if /^\s*$/;
		my ($id,@seqs)=split /\n+/;
		($id)         =$id=~/^(\S+)/;
		my $seq       =join "",@seqs;
		$fa{$id}      =$seq;
	}
	close IN;           $/="\n";
}

sub useful
{## For map return value
	my ($in,$target)=@_;
	if($in->[-1] eq $target)
	{
		return $in;
	}
}



sub USAGE {#
	my $usage=<<"USAGE";
Program:	$0
Version:	$version
Usage:
  Options:

	-mrd      <str>	input mrd file 
	-samples  <str>	samples,"S01,S02",forced 
	-od       <str>	output dir
	-known    <str>	known miRNA fa file
	-Tscore         Total score of miRDeep2, default 0
	-Tmfe           the score of mfe(minimal free energy), default 0.4
	-Tstar          score for star reads, default -2
	-Mcount   <num>	threshold of mature reads count, default 5
	-coord          precursor.coord
	-str            precursor.str
	-h              Help

USAGE
	print $usage;
	exit;
}
