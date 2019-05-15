#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Spreadsheet::WriteExcel;
use Spreadsheet::ParseExcel;  
use Encode;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($id,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"id:s"=>\$id,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($id and $od);


#############################转换成excel表格
#my $stat="$od/Project.stat.xls";
#$stat = Spreadsheet::WriteExcel->new("$stat");
#my @filter=glob ("$od/*.xls");
#foreach my $f (@filter) 
#{
#	my $name=basename($f);
#	next if($name=~/^Project/);
#	$name=~s/\.xls//;
#	$name=~s/.*\.fa\.//;
#	Beta($f,$stat,$name,0);
#}
#
#


my $parser   = Spreadsheet::ParseExcel->new();
my $workbook = $parser->parse("Project.stat.xls");

if ( !defined $workbook ) {
	die $parser->error(), ".\n";
}

for my $worksheet ( $workbook->worksheets() ) {
	my $True_sheet = $worksheet->get_name();
	print "$True_sheet\n";
	my ( $row_min, $row_max ) = $worksheet->row_range();
	my ( $col_min, $col_max ) = $worksheet->col_range();

	for my $row ( $row_min .. $row_max ) {
		for my $col ( $col_min .. $col_max ) {

			my $cell = $worksheet->get_cell( $row, $col );
			next unless $cell;

			print "Row, Col    = ($row, $col)\n";
			print "Value       = ", $cell->value(),       "\n";
			print "Unformatted = ", $cell->unformatted(), "\n";
			print "\n";die;
		}
	}
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
sub num_format{
	my $number=shift @_;
	1 while ($number=~s/^(-?\d+)(\d\d\d)/$1,$2/);
	$number;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub Beta #将linux下的xls文件转换成windows下可读的excel文件
{
	my ($file,$xls,$sheet,$f)=@_;
	open (IN,"$file")||die "can't open file $file\n";
	my $worksheet = $xls->add_worksheet(decode('GB2312',"$sheet"));
	my $format = $xls->add_format();
	$format->set_border($f);
	$format->set_align('left');
#	$format->set_font('Times New Roman');
	$format->set_size('11') if($f ==0);
	$format->set_font('宋体') if($f ==0);
	my $num=0;
	while(<IN>)
	{
		chomp;
		my @line=split/\t/,$_;
		for(my $i=0;$i<@line;$i++)
		{
			$format->set_num_format('#,###') if($f==1 && $line[$i]=~/\d+/ && $line[$i]!~/[\.\%]/);
			$worksheet->write($num,$i,decode('GB2312',"$line[$i]"),$format);
		}
		$num++;
	}
	close IN;
}


sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	Zhang XueChuan <zhangxc\@biomarker.com.cn> 
Usage:
  Options:
  -id    <dir>  input dir,forced 
  
  -od    <dir>  output dir,forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
