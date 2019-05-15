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
my $template_dir = "template";

my ($pid,$tid)=(0,0);

my $report_time=&gaintime();
my $writer = XML::Writer->new(OUTPUT => 'self');
$writer->xmlDecl('UTF-8');
$writer->startTag('report');
$writer->emptyTag('report_version','value'=>$config{report_version});
my $report_name;
$report_name=$config{report_name}	if($config{report_name});
$report_name ||= "测序数据评估与统计";
$writer->emptyTag('report_name','value'=>$report_name);
$writer->emptyTag('report_code','value'=>$config{report_code});
$writer->emptyTag('report_user','value'=>$config{report_user});
$writer->emptyTag('report_user_addr','value'=>$config{report_user_addr});
$writer->emptyTag('report_time','value'=>"XXXX/XX/XX;XXXX/XX/XX;XXXX/XX/XX;$report_time");
$writer->emptyTag('report_abstract','value'=>"");
###############################摘要
$writer->emptyTag('h1','name'=>"原始数据介绍",'type'=>'type1','desc'=>"原始数据介绍");

$writer->emptyTag('p','desc'=>"基于边合成边测序（Sequencing By Synthesis，SBS）技术，Illumina HiSeq高通量测序平台对cDNA文库进行测序，产出大量的高质量Data，称为原始数据（Raw Data）。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"Raw Data通常以FASTQ格式提供，每个测序样品的Raw Data包括两个FASTQ文件，分别包含所有cDNA片段两端测定的Reads。",'type'=>"type1");
$writer->emptyTag('p','desc'=>"FASTQ格式文件示意图如下：",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"注：FASTQ格式文件中每个Read由四行描述，其中第一行以“@”开头，随后为Illumina测序识别符（Sequence Identifiers）和描述文字（选择性部分）；第二行是碱基序列；第三行以“+”开头，随后为Illumina测序识别符（选择性部分）；第四行是对应序列的测序质量。",'name'=>"图$pid. FASTQ格式文件示意图",'type'=>"type1",'path'=>"$template_dir/fastq.png");
$writer->emptyTag('h1','name'=>"1 测序碱基质量值",'type'=>'type1','desc'=>"1 测序碱基质量值");
$writer->emptyTag('p','desc'=>"碱基质量值（Quality Score或Q-score）是碱基识别（Base Calling）出错的概率的整数映射。通常使用的碱基质量值Q公式<a href=\"#ref1\">[1]</a>为：",'type'=>"type1");
$pid++;
$writer->emptyTag('pic','desc'=>"",'name'=>"图$pid. 质量值计算公式",'type'=>"type1",'path'=>"$template_dir/Qscore.png");
$writer->emptyTag('p','desc'=>"其中P为碱基识别出错的概率。下表给出了碱基质量值与碱基识别出错的概率的对应关系：",'type'=>"type1");
$tid++;
$writer->emptyTag('table','desc'=>'注：碱基质量值越高表明碱基识别越可靠，准确度越高。比如，对于碱基质量值为Q20的碱基识别，100个碱基中有1个会识别出错，以此类推。','type'=>"4",'name'=>"表$tid. 碱基质量值与碱基识别出错的概率的对应关系表",'path'=>"$template_dir/error.txt");

my ($lncRNA,$sRNA,@sRNA);
if (-d "$in/BMK_1_rawData/BMK_1_Data_Assess/miRNA_PNG"){
	$lncRNA = "$in/BMK_1_rawData/BMK_1_Data_Assess/nonRib_PNG";
	$sRNA = "$in/BMK_1_rawData/BMK_1_Data_Assess/miRNA_PNG";
}else{
	$lncRNA = "$in/BMK_1_rawData/BMK_1_Data_Assess/PNG";
	$sRNA = "$in/BMK_1_rawData/BMK_1_Data_Assess/PNG";
	@sRNA = glob "$in/BMK_1_rawData/BMK_1_Data_Assess/PNG/*error.png";
}

if(-d "$in/BMK_1_rawData/BMK_1_Data_Assess/miRNA_PNG"){
	$writer->emptyTag('p','desc'=>"lncRNA样品碱基测序错误率分布见下图：",'type'=>"type1");
	&piclist("碱基测序错误率分布图","注：横坐标为Reads的碱基位置， 纵坐标为单碱基错误率。测序错误率受测序仪本身、测序试剂、样品等多个因素共同影响。由于测序过程中化学试剂的消耗，测序错误率会随着测序序列(Sequenced Reads)长度的增加而升高。此特征是Illumina高通量测序平台所具有的。","$lncRNA/*.quality.png","BMK_1_rawData");
	$writer->emptyTag('p','desc'=>"sRNA样品碱基测序错误率分布见下图：",'type'=>"type1");
	&piclist("碱基测序错误率分布图","注：横坐标为reads的碱基位置，纵坐标为每个测序反应的碱基错误率的平均值。","$sRNA/*.error.png","BMK_1_rawData");
	$writer->emptyTag('h1','name'=>"2 测序碱基分布含量",'type'=>'type1','desc'=>"2 测序碱基分布含量");
	$writer->emptyTag('p','desc'=>"碱基类型分布检查用于检测有无AT、GC分离现象,由于RNA-Seq所测的序列为随机打断的cDNA片段，因随机性打断及碱基互补配对原则，理论上，G和C、A和T的含量每个测序循环上应分别相等，且整个测序过程稳定不变，呈水平线。由于Reads 5’端的前几个碱基为随机引物序列存在一定的偏好性，因此会在碱基分布图中出现前端波动较大的现象。",'type'=>"type1");
	&piclist("ATCG含量分布图","注：横坐标为Reads的碱基位置，纵坐标为单碱基所占比例。","$lncRNA/*.acgtn.png","BMK_1_rawData");
}else{
	$writer->emptyTag('p','desc'=>"样品碱基测序错误率分布见下图：",'type'=>"type1");
	if(@sRNA>0){
		&piclist("碱基测序错误率分布图","注：横坐标为reads的碱基位置，纵坐标为每个测序反应的碱基错误率的平均值。","$sRNA/*.error.png","BMK_1_rawData");
	}else{
		&piclist("碱基测序错误率分布图","注：横坐标为Reads的碱基位置， 纵坐标为单碱基错误率。测序错误率受测序仪本身、测序试剂、样品等多个因素共同影响。由于测序过程中化学试剂的消耗，测序错误率会随着测序序列(Sequenced Reads)长度的增加而升高。此特征是Illumina高通量测序平台所具有的。","$lncRNA/*.quality.png","BMK_1_rawData");
		$writer->emptyTag('h1','name'=>"2 测序碱基分布含量",'type'=>'type1','desc'=>"2 测序碱基分布含量");
		$writer->emptyTag('p','desc'=>"碱基类型分布检查用于检测有无AT、GC分离现象,由于RNA-Seq所测的序列为随机打断的cDNA片段，因随机性打断及碱基互补配对原则，理论上，G和C、A和T的含量每个测序循环上应分别相等，且整个测序过程稳定不变，呈水平线。由于Reads 5’端的前几个碱基为随机引物序列存在一定的偏好性，因此会在碱基分布图中出现前端波动较大的现象。",'type'=>"type1");
		&piclist("ATCG含量分布图","注：横坐标为Reads的碱基位置，纵坐标为单碱基所占比例。","$lncRNA/*.acgtn.png","BMK_1_rawData");
	}
}
############################
$writer->startTag('ref_list',desc=>"参考文献",type=>"type1",name=>"参考文献");
&reference ($writer,[
        ["Kim D, Langmead B, Salzberg S L. HISAT: a fast spliced aligner with low memory requirements.[J]. Nature Methods, 2015, 12(4):357-360.","https://www.nature.com/articles/nmeth.3317","1"],
]);
$writer->endTag('ref_list');


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
        	$writer->emptyTag('pic','desc'=>"",'name'=>"$base",'type'=>"type1",'path'=>"../$new/$base");
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
		my $new = $path.$tmp;
        	$writer->emptyTag('file','desc'=>"",'name'=>"$base",'action'=>"xls",'path'=>"../$new/$base",'type'=>"xls");
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

