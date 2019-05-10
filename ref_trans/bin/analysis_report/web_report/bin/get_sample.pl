use strict;
use warnings;
use Cwd qw(abs_path); #获取工作路径，即当前目标所在的路径（函数括号里的东西）
use Getopt::Long; #获取选项
use Data::Dumper; #可打印引用的东西。eg：print Dumper(\%hash \@array);
use FindBin qw($Bin $Script); #$Bin  调用脚本的bin目录的路径，$Script  脚本名称  $RealBin 调用脚本的绝对路径  $RealScript  与脚本相关的脚本（用不着）
use File::Basename qw(basename dirname); #basename函数获取文件名  dirname函数获取路径  fileparse函数获取扩展名
my $BEGIN_TIME = time(); #获取系统时间，可两者做差得到运行时间，单位为秒（s）
my $version = "1.0.0";
use newPerlBase;
my %config=%{readconf("$Bin/../../../../config/db_file.cfg")}; 
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#GetOptions
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
my ( $input, $list, $species, $output,$cloud,%species_for_database )
  ;    #选项的hash表
GetOptions(
    "input=s" => \$input,   
    "list=s" => \$list,  
    "output:s"   => \$output,
    "species:s"   => \$species,
    "cloud:s"   => \$cloud,
    "help|?" => \&USAGE,
) or &USAGE;
&USAGE unless ( defined $input and defined $list and defined $species );    #判断选项的合法性
&log_current_time("$Script start.....");

open READ_2, "$list" || die "can't open the file:$_\n";
while (<READ_2>) {
    chomp;
    my @arr = split /\s+/, $_;
    $species_for_database{ $arr[0] } = [ @arr[ 1 .. 3 ] ];
}
###############2015/11/18 by luml
unless ( ${ $species_for_database{$species} }[0] ) {
    print "$species is not exist in $list!";
    #	next;
}
else {
    my $database = ${ $species_for_database{$species} }[0];
    my $dataset  = ${ $species_for_database{$species} }[1];
    my $host     = ${ $species_for_database{$species} }[2]
      ;    ##biomaRt change 2015-11-10 by linhj
    if (defined $cloud) {
        system "$config{Rscript} $Bin/id_conversion.r --infile $input --outfile $output --id_name ensembl_gene_id --id.col 1 --database $database --dataset $dataset --host $host";    
    } else {
        system "ssh $config{net_node}  $config{Rscript} $Bin/id_conversion.r --infile $input --outfile $output --id_name ensembl_gene_id --id.col 1 --database $database --dataset $dataset --host $host";        

    }
}


###################################################################################################################################################################################
&log_current_time("$Script end……");    #调用时间函数
my $run_time = time() - $BEGIN_TIME;
print "the program run time is:$run_time.s\n";

sub log_current_time {

    # get parameter    #获取参数
    my ($info) = @_;

    # get current time with string   #获取当前串的时间
    my $curr_time = &date_time_format( localtime( time() ) )
      ;                                    #格式化获取时间表示方法

    # print info with time
    print "[$curr_time] $info\n";
}

sub date_time_format {
    my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) =
      localtime( time() );
    return sprintf(
        "%4d-%02d-%02d %02d:%02d:%02d",
        $year + 1900,
        $mon + 1, $day, $hour, $min, $sec
    );
}

###################################################################################################################################################################################
sub USAGE {
    my $usage = <<"USAGE";
----------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: huangxy <huangxy\@biomarker.com.cn>
      Date: 2015-08-13

     Usage:
            --input     <FILE>   input file which need add gene name 
            --list      <FILE>   Must(local_and_biomaRt_database.txt)
            --output    <FILE>   result file which include gene name 
            --species      <FILE>   species name for get gene name 
                          some name from  $Bin/bin/local_and_biomaRt_database.txt  force
            --cloud     <FILE>   analysis at biocloud
     Example:
            perl $Script --input ref_trans_full_table.xls  --list local_and_biomaRt_database.txt

----------------------------------------------------------------------------------------------------
USAGE
    print $usage;
    exit;
}
