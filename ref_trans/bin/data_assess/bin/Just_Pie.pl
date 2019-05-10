#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use SVG::TT::Graph::Pie; 
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my %config=%{readconf("$Bin/../../../config/db_file.cfg")}; 
my $version="1.0.0";
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($inf,$fout,$fheight,$fwidth,$fcss,$fnote);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$inf,
				"o:s"=>\$fout,
				"h:i"=>\$fheight,
				"w:i"=>\$fwidth,
				"css:s"=>\$fcss,
				"note:s"=>\$fnote,
				) or &USAGE;
&USAGE unless ($inf and $fout);
$fcss ||= "$Bin/pie12.css";
my @fields;
my @datas;
my @Rate;
my $data_sum;
my $svgxxx = $config{svg2xxx2};
open IN,$inf or die "cannot open file $inf, $!\n";
while (<IN>) {
	chomp;
	next if (/^$/||/\#/);
	my @lines = split /\t/,$_;
	push @fields,$lines[0];
	push @datas,$lines[1];
	$data_sum += $lines[1];
}  
close IN;
if ($data_sum == 0) {
	die "please check file ,some wrong cat data volumn\n";
}
$fheight ||= 300;
$fwidth  ||= 500;
my $graph;
if (defined $fnote) {
	$graph = SVG::TT::Graph::Pie->new( 
	{ 
		'height'              => $fheight, 
		'width'               => $fwidth, 
		'fields'              => \@fields, 
		'show_graph_title'    => 1,
		'graph_title'         => "$fnote",
		'show_data_labels'    => 0,
		'show_actual_values'  => 0,
		'show_percent'        => 1,
		'rollover_values'     => 1,
		'key'                 => 1,
		'key_placement'       => 'R',
		'rollover_values'     => 1,
		'style_sheet'         => "$fcss", # internal stylesheet
		'random_colors'       => 0,
	  } 
	);
}else{
	$graph = SVG::TT::Graph::Pie->new( 
	{ 
		'height'              => $fheight, 
		'width'               => $fwidth, 
		'fields'              => \@fields, 
		'show_data_labels'    => 0,
		'show_actual_values'  => 0,
		'show_percent'        => 1,
		'rollover_values'     => 1,
		'key'                 => 1,
		'key_placement'       => 'R',
		'rollover_values'     => 1,
		'style_sheet'         => "$fcss", # internal stylesheet
		'random_colors'       =>0,
	  } 
	);
}
 
$graph->add_data( 
  { 
      'data'  => \@datas, 
      'title' => 'datas', 
  } 
); 
  
open( my $fh, '>', "$fout" ) or die "cannot open file $fout, $!\n"; 
select $fh; 
binmode $fh; 
print $graph->burn(); 
close($fh); 
if ($fout=~/Raw_data_ratio/) {
	open IN,"$fout" ||die $!;
	open OUT,">$fout.svg" ||die $!;
	while (<IN>) {

		if ($_=~/keyText/) {
			$_=~s/\]/%\]/;
		}
		print OUT $_;
	}
	close(IN);
	close(OUT);
my $cmd=`mv $fout.svg $fout `;
}
my $nowdir = Cwd::realpath();
my $newdir = dirname($fout);
my $filename = basename($fout);
chdir($newdir);
` $svgxxx ./$filename -t png`;
chdir($nowdir);

sub USAGE{
	print <<"	Usage End.";
	Description:
		Writer by sunhy <sunhy\@biomarker.com.cn>
	Usage:
		-i    <file>     input file                              [required]
		-o    <file>     output file (svg file)                  [required]

		-w    <int>      the graph width  default[500]           [optional]
		-h    <int>      the graph height default[300]           [optional]
		-css  <file>     the css file that you want to use       [optional]
				 default [\$Bin/pie.css]
		-note <note>	 a main title of graph                   [optional]
		-h               Help document
	Usage End.
	exit;
}
  

  
