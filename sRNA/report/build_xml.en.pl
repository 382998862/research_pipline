#!/usr/bin/perl
use Getopt::Long;
use Getopt::Std;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);
use File::Basename qw(basename dirname);
use XML::Writer;
use IO::File;
use Encode;
my ($o,$in);
GetOptions(
        "h|?"           =>\&USAGE,
	"in:s"		=>\$in,
        "o:s"           =>\$o,
)or &USAGE;
&USAGE unless ($in);

`mkdir -p $in/HTML/dataassess`	unless(-d "$in/HTML/dataassess");
$o ||= "$in/dataassess.xml";

$in=abs_path($in);
$o =abs_path($o);
my $od=dirname($o);

my ($pid,$tid)=(0,0);

my $report_time=&gaintime();
my $writer = XML::Writer->new(OUTPUT => 'self');
$writer->xmlDecl('UTF-8');
$writer->startTag('report');
$writer->emptyTag('report_version','value'=>$config{report_version});
my $report_name;
$report_name=$config{report_name}	if($config{report_name});
$report_name ||= "Material and Methods";
$writer->emptyTag('report_name','value'=>$report_name);
$writer->emptyTag('report_code','value'=>$config{report_code});
$writer->emptyTag('report_user','value'=>$config{report_user});
$writer->emptyTag('report_user_addr','value'=>$config{report_user_addr});
$writer->emptyTag('report_time','value'=>"XXXX/XX/XX;XXXX/XX/XX;XXXX/XX/XX;$report_time");
$writer->emptyTag('report_abstract','value'=>"");
###############################摘要
$writer->emptyTag('h1','name'=>"1 Sample collection and preparation",'type'=>'type1','desc'=>"1 Sample collection and preparation");
$writer->emptyTag('h1','name'=>"",'type'=>'type1','desc'=>"");
$writer->emptyTag('p','desc'=>"");
$writer->emptyTag('p','desc'=>"",'type'=>"type1");
$writer->emptyTag('p','desc'=>"",'type'=>"type1");




$writer->endTag('report');

open OUT,">:utf8", "$o";
my $xmlstr=&decorate($writer->to_string);
$xmlstr = Encode::decode("utf-8", $xmlstr);
print OUT $xmlstr;
close(OUT);
$writer->end();
print "/share/nas2/genome/biosoft/Python/2.7.8/bin/python $Bin/bin/xml2HtmlConverter.py -i $o -o $in/HTML/dataassess -n dataassess\n";
`/share/nas2/genome/biosoft/Python/2.7.8/bin/python $Bin/bin/xml2HtmlConverter.py -i $o -o $in/HTML/dataassess -n dataassess`;

################################################################################################
#
#		                        sub function
################################################################################################

sub decorate{
	my $xmlstr=shift;
	$xmlstr=~s/(<[^>]*)\>/$1\>\n/mgo;
	return $xmlstr;
}

sub reference  {
    my ($writer,$list)=@_;
    my ($data,$value,$unit,$item);
    
    foreach $data (@$list)    {		
	$writer->emptyTag('ref','name'=>$$data[0],'link'=>$$data[1],'id'=>$$data[2]);
    }    
}

sub piclist{
	my ($name,$desc,$pics,$path)=@_;
	$pid++;
	$writer->startTag('pic_list','name'=>"图$pid. $name",'desc'=>"$desc",'type'=>"type1");
	my @images=glob("$pics");
	foreach my $s(@images){
		my $base=basename $s;
		my $dir=dirname $s;
		my $tmp=(split(/$path/,$dir))[1];
		my $new =$path.$tmp;
        	$writer->emptyTag('pic','desc'=>"",'name'=>"$base",'type'=>"type1",'path'=>"$new/$base");
	}
	$writer->endTag('pic_list');
}

sub filelist{
	my ($desc,$pics,$path)=@_;
	$writer->startTag('file_list','name'=>"",'desc'=>"$desc",'type'=>"type1");
	my @images=glob("$pics");
	foreach my $s(@images){
		my $base=basename $s;
		my $dir=dirname $s;
		my $tmp=(split(/$path/,$dir))[1];
		$path .=$tmp;
        	$writer->emptyTag('file','desc'=>"",'name'=>"$base",'action'=>"xls",'path'=>"$path/$base",'type'=>"xls");
	}
	$writer->endTag('file_list');
}
sub relate{    
        my $file=shift;
        open(REL,$file)||die $!;
        while(<REL>){
                chomp;
                next if($_ !~/^SampleID/);
		$_=~s/\s+$//;
                my @tmp=split(/\s+/,$_);
                $sample{$tmp[-1]}{miRNA}=$tmp[2];
                $sample{$tmp[-1]}{lncRNA}=$tmp[3];
                $sample{$tmp[-1]}{circRNA}=$tmp[3];
                $sample{$tmp[-1]}{mRNA}=$tmp[3];
		$sample{$tmp[-1]}{user}=$tmp[1];
        }
        close(REL);
}

sub gaintime{
	my $timestamp=time(); 
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($timestamp); 
	my $y = $year + 1900; 
	my $m = $mon + 1; 
	$timestamp=sprintf("%4d-%02d-%02d",$y,$m,$mday);
	return $timestamp
}




sub USAGE{
        my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
        Options:
	-in	<path>	input BMK_Result path, forced
        -o      <path>	output xml name, default inpath/assess.xml

        -h      Help

Example: perl $0 -in BMK_Result -o HTML/assess.xml

USAGE
        print $usage;
        exit;
}

