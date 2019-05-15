#!/usr/bin/perl
use strict;
use warnings;
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME = time();
my $version    = "1.0.0";
my ( $fIn, $dir, $sam );
## input is CIRI output file;
GetOptions(
	"help|?"   => \&USAGE,
	"dir:s"    => \$dir,
	"i:s"      => \$fIn,
	"sample:s" => \$sam
) or &USAGE;
&USAGE unless ( $fIn and $dir and $sam );
open( IN, $fIn ) or die $!;

if ( !-d $dir ) {
	`mkdir $dir`;
}
open( OUT,  ">$dir/All_gene_counts.list" )       or die $!;
open( OUT1, ">$dir/All_gene_counts_detail.xls" ) or die $!;

my @sam = split /,/, $sam;
$sam = join( "\t", @sam );

print OUT "circRNA_ID\t$sam\tgeneLength\n";
print OUT1 "circRNA_ID\t$sam\tgeneLength\tgene_id\n";

while (<IN>) {
	next if (/\#/);
	next if (/^intergenic/);
	my @line        = split /\t/, $_;
	my $length      = $line[3] - $line[2] + 1;
	my $junction_id = $line[-1];
	my @junction_id = split /,/, $junction_id;
	my %expression;
	## count every sample junctions reads
	foreach my $id (@junction_id) {
		if ( $id =~ /\s+/ ) {
			next;
		}
		$id = ( split( /:/, $id ) )[0];
		if ( exists $expression{$id} ) {
			$expression{$id} = $expression{$id} + 1;
		}
		else {
			$expression{$id} = 1;
		}

	}
	
	my $flag=1;
	my $gene_list = $line[9];
	$gene_list=~s/,$//;
	if ( $gene_list =~ /,/ ) {
		my @gene_id = split /,/, $gene_list;
		for ( my $i = 0 ; $i < @gene_id ; $i++ ) {
			if($flag ==1)
			{
				print OUT "$line[0]\t";
			}
			print OUT1 "$line[0]\t";
			foreach (@sam) {
				if ( !exists $expression{$_} ) {
					if($flag ==1)
					{
						print OUT "0\t";
					}
					print OUT1 "0\t";
				}
				else {
					my $value = $expression{$_};
					if($flag ==1)
					{
						print OUT "$value\t";
					}
					print OUT1 "$value\t";
				}
			}
			if($flag ==1)
			{
				print OUT "$length\n";
			}
			print OUT1 "$length\t$gene_id[$i]\n";
			$flag++;
		}
	}
	else {
		print OUT "$line[0]\t";
		print OUT1 "$line[0]\t";
		foreach (@sam) {
			if ( !exists $expression{$_} ) {
				print OUT "0\t";
				print OUT1 "0\t";
			}
			else {
				my $value = $expression{$_};
				print OUT "$value\t";
				print OUT1 "$value\t";
			}
		}
		print OUT "$length\n";
		print OUT1 "$length\t$gene_list\n";

	}
}

close IN;
close OUT;
close OUT1;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ", time() - $BEGIN_TIME, "s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

sub ABSOLUTE_DIR {    #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir = `pwd`;
	chomp($cur_dir);
	my ($in) = @_;
	my $return = "";
	if ( -f $in ) {
		my $dir  = dirname($in);
		my $file = basename($in);
		chdir $dir;
		$dir = `pwd`;
		chomp $dir;
		$return = "$dir/$file";
	}
	elsif ( -d $in ) {
		chdir $in;
		$return = `pwd`;
		chomp $return;
	}
	else {
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

################################################################################################################

sub max {    #&max(lists or arry);
	         #求列表中的最大值
	my $max = shift;
	my $temp;
	while (@_) {
		$temp = shift;
		$max = $max > $temp ? $max : $temp;
	}
	return $max;
}

################################################################################################################

sub min {    #&min(lists or arry);
	         #求列表中的最小值
	my $min = shift;
	my $temp;
	while (@_) {
		$temp = shift;
		$min = $min < $temp ? $min : $temp;
	}
	return $min;
}

################################################################################################################

sub revcom() {    #&revcom($ref_seq);
	 #获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
	my $seq = shift;
	$seq =~ tr/ATCGatcg/TAGCtagc/;
	$seq = reverse $seq;
	return uc $seq;
}

################################################################################################################

sub GetTime {
	my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) =
	  localtime( time() );
	return sprintf(
		"%4d-%02d-%02d %02d:%02d:%02d",
		$year + 1900,
		$mon + 1, $day, $hour, $min, $sec
	);
}

################################################################################################################
sub USAGE {
	my $usage = <<"USAGE";
 ProgramName:
     Version:   $version
     Contact:   Simon Young <yangxh\@biomarker.com.cn> 
Program Date:   2012.07.02
      Modify:   
 Description:   This program is used to ......
       Usage:
        Options:
        -i <file>   input file,xxx format,forced
        -dir <file>   output dir,forced
	-sample <str> sample name list split by comma
        -h      help

USAGE
	print $usage;
	exit;
}
