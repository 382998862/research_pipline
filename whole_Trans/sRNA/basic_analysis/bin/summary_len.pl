use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0";
my  ($od,$fyes,$min,$max,$seq);
GetOptions(
	"help|?"=>\&USAGE,
	"dir:s"=>\$od,
	"fyes:s"=>\$fyes,
	"min:s"=>\$min,
	"max:s"=>\$max,
	"seq:s"=>\$seq,
	) or &USAGE;
&USAGE unless ($od);
&USAGE unless (defined $fyes and $seq);

$od=&ABSOLUTE_DIR($od);


$min||=18;
$max||=30;

my %Raw_FQ;
my @raw=split/,/,$seq;
foreach(@raw){
	$Raw_FQ{$_}=1;
}
	my $total_mirna_known=0;my $total_mirna_predicted=0;
	if($fyes==1){
		open(IN1,"$od/miRDeep2/miRNA_Quantify/sum_known_miRNA.txt")||die $!;
		
		my $first=<IN1>;
		my @first=split/\s+/,$first;
		splice(@first,0,11);
		#pop @first;   ##### do not have target genes
		my %sample_len;
		while(<IN1>){
			$total_mirna_known++;
			my @temp=split/\s+/,$_;
			my $len=length($temp[2]);
			splice(@temp,0,11);
			#pop @temp;
			#print "@temp\n";
			for(my $i=0;$i<=$#temp;$i++){
				if($temp[$i]>0){
					$sample_len{$first[$i]}{$len}++;
				}
			}
		}
		close IN1;

		foreach(keys %sample_len){
			my $sample=$_;
			open(OUT1,">$od/miRDeep2/miRNA_Quantify/$_.known_miRNA_len.stat")||die $!;
			print OUT1 "Length\tNumber\n";
			for(my $i=$min;$i<=$max;$i++){
				if(exists $sample_len{$sample}{$i}){
					print OUT1 "$i\t$sample_len{$sample}{$i}\n";
				}
				else{
					print OUT1 "$i\t0\n";
				}
			}
			close OUT1;
		}
	}
	
	
	open(IN1,"$od/miRDeep2/miRNA_Quantify/sum_predicted_miRNA.txt")||die $!;
	my $first=<IN1>;
	my @first=split/\s+/,$first;
	splice(@first,0,11);
	#pop @first;
	my %sample_len;
	while(<IN1>){
		$total_mirna_predicted++;
		my @temp=split/\s+/,$_;
		my $len=length($temp[2]);
		splice(@temp,0,11);
		#pop @temp;
		for(my $i=0;$i<=$#temp;$i++){
			if($temp[$i]>0){
				$sample_len{$first[$i]}{$len}++;
			}
		}
	}
	close IN1;
	foreach(keys %sample_len){
			my $sample=$_;
			print "$sample\n";
			open(OUT1,">$od/miRDeep2/miRNA_Quantify/$_.predicted_miRNA_len.stat")||die $!;
			print OUT1 "Length\tNumber\n";
			for(my $i=$min;$i<=$max;$i++){
				if(exists $sample_len{$sample}{$i}){
					print OUT1 "$i\t$sample_len{$sample}{$i}\n";
				}
				else{
					print OUT1 "$i\t0\n";
				}
			}
			close OUT1;
	}
	
	
	open(TOTAL,">$od/miRDeep2/miRNA_Quantify/Total_miRNA.stat")||die $!;
	
	print TOTAL "BMK-ID\tKnown-miRNAs\tNovel-miRNAs\tTotal\n";

	foreach my $name (sort keys %Raw_FQ) {
		my $known_mirna=0;my $novel_mirna=0;
#		open(IN1,"$od/Alignment/$name/Len_stat/Total.stat")||die $!;
		print "$od/miRDeep2/miRNA_Quantify/$name.predicted_miRNA_len.stat\n";
		open(IN3,"$od/miRDeep2/miRNA_Quantify/$name.predicted_miRNA_len.stat")||die "$od/miRDeep2/miRNA_Quantify/$name.predicted_miRNA_len.stat do not exist";
		my %total_length;my %total_length_known;my %total_length_predicted;
		if($fyes==1 and -f "$od/miRDeep2/miRNA_Quantify/$name.known_miRNA_len.stat"){
			open(IN2,"$od/miRDeep2/miRNA_Quantify/$name.known_miRNA_len.stat")||die $!;
			<IN2>;
			while(<IN2>){
				chomp;
				my @temp=split/\s+/,$_;
				$total_length{$temp[0]}+=$temp[1];
				$total_length_known{$temp[0]}=$temp[1];
			}
			close IN2;
			foreach(keys %total_length_known){
				$known_mirna+=$total_length_known{$_};
			}
		}
		<IN3>;
		while(<IN3>){
			chomp;
			my @temp=split/\s+/,$_;
			$total_length{$temp[0]}+=$temp[1];
			$total_length_predicted{$temp[0]}=$temp[1];
		}
		foreach(keys %total_length_predicted){
			$novel_mirna+=$total_length_predicted{$_};
		}
		
		foreach(keys %total_length){
			$total_length{Total}+=$total_length{$_};
		}
		close IN3;#close OUT;
		my $k_n=$known_mirna+$novel_mirna;
		print TOTAL "$name\t$known_mirna\t$novel_mirna\t$k_n\n";
	}
	my $k_n=$total_mirna_known+$total_mirna_predicted;
	print TOTAL "Total\t$total_mirna_known\t$total_mirna_predicted\t$k_n\n";
	close TOTAL;
	
	
######################## sub
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
sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	Liu xiaoshuang <liuxs\@biomarker.com.cn> 
Usage:
  Options:
  -dir			<file>   file dir,forced 
  -fyes   		<num>    have known miRNA,1 ; orthwise ,0
  -min			<num>    the min of sRNA length
  -max         <max>	the max length of sRNA
  -seq          <str>   the name of sample
  -h			Help

USAGE
	print $usage;
	exit;
}
