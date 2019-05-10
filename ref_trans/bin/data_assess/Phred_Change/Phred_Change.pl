#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
use newPerlBase;
#my %config=%{readconf("$Bin/../../../config/db_file.cfg")}; 
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($config,$qu,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"cfg:s"=>\$config,
                "od:s"=>\$od
				) or &USAGE;
&USAGE unless ($config);
my (%data_cfg);
my $cmd;
mkdir "$od/data" unless -d "$od/data";
&data_cfg_read($config,\%data_cfg);
system("cp $config $config.bak");
open(CONF,">$config")||die $!;
my %All_Q;
foreach my $sample(keys %{$data_cfg{rawdata}}){
    foreach my $f(keys %{$data_cfg{rawdata}{$sample}}){
        my $file1=$data_cfg{rawdata}{$sample}{$f};
        #my $file2=$data_cfg{rawdata}{$sample}{fq2};
        $cmd = "perl $Bin/fastq_phred.pl -n 1000 -a -od $od  $file1";
        &run_or_die($cmd);
        my $phred_file=(glob("$od/*.phred"))[0];
        my ($fq,$Q); 
        open(IN,"$phred_file")||die $!;
        while(<IN>)
        {
            next if(/^\s*$/);
            chomp;
            $Q=$_ ;			
        }
        close IN;
        push @{$All_Q{$Q}},$file1;
        
        system("rm $od/*.phred");		   
    }
}

my @all_q=keys %All_Q;
my $Final_Q;
if ($#all_q==0) {
    $Final_Q=$all_q[0];
    
}
elsif($#all_q==1){
    $Final_Q=33;
    my @file=@{$All_Q{64}};
    my @file_33=@{$All_Q{33}};
    if ($#file<0 or $#file_33<0){
        my $all_q=join(";",@all_q);
        die "Unknow file type $all_q\n";
    }
    foreach my $file1(@file){
            my $filename=basename $file1;
            my $newfile="$od/data/$filename";
            my $sample_id=$data_cfg{filename}{sample}{$file1};
            my $fq=$data_cfg{filename}{fq}{$file1};
            $data_cfg{rawdata}{$sample_id}{$fq}=$newfile;
            print"######111######raw:$file1\n#######111#######fq:$fq\n";
            my $cmd="/share/nas2/genome/bmksoft/tool/fastq_quality_convert/v1.0/fastq_quality_convert -i $file1 -o $newfile -f 1 -t 0\n";
            my $flag = system($cmd) ;
            if ($flag != 0)
                {
                $cmd="/share/nas2/genome/bmksoft/tool/fastq_quality_convert/v1.0/fastq_quality_convert -i $file1 -o $newfile -f 2 -t 0\n";
                &run_or_die($cmd) ;
                }  
    }
}
else{
    print "There are three kind of Qphred,please Check\n";
    die $!;
}

print CONF "Qphred\t$Final_Q\n";
foreach my $sample_id(sort keys %{$data_cfg{rawdata}}){
    print CONF "Sample\t$sample_id\n";
    print CONF "fq1\t$data_cfg{rawdata}{$sample_id}{fq1}\n";
    if (exists $data_cfg{rawdata}{$sample_id}{fq2}) {
        print CONF "fq2\t$data_cfg{rawdata}{$sample_id}{fq2}\n";
    }
}
if (exists $data_cfg{Basecall_stat}) {
     print CONF "Basecall_stat\t$data_cfg{Basecall_stat}\n";
}
close CONF;



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub data_cfg_read {
    #&log_current_time("data config check:");
    my ($cfg_file, $data_cfg) = @_;
    my $sample_id;

    open (CFG, $cfg_file) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+$/ or /^#/);
        if (/^Qphred/) {
            $data_cfg->{Qphred} = (split /\s+/,$_)[1];
        }
        if (/^Sample/) {
            $sample_id=(split/\s+/,$_)[1];
        }
        if ($_=~/^fq1\s+/ || $_=~/^fq2\s+/) {
            my $file=(split/\s+/,$_)[1];
            die "$file is not exist!\n" unless (-e $file);

            $data_cfg->{rawdata}{$sample_id}{fq1}=$file if $_=~/^fq1/;
            $data_cfg->{rawdata}{$sample_id}{fq2}=$file if $_=~/^fq2/;
            $data_cfg->{filename}{fq}{$file}="fq1" if $_=~/^fq1/;
            $data_cfg->{filename}{fq}{$file}="fq2" if $_=~/^fq2/;
            if (exists $data_cfg{filename}{sample}{$file}) {
                print "Error,$file is used twice\n";
                die;
            }
            
            $data_cfg->{filename}{sample}{$file}=$sample_id;
        }
		if ($_=~/^Basecall_stat/) {
            my $file=(split/\s+/,$_)[1];
            die "$file is not exist!\n" unless (-e $file);
            $data_cfg->{Basecall_stat}= $file;
        }
    }
    close CFG;

    

    $data_cfg->{sample_num} = scalar keys %{$data_cfg->{rawdata}};
    print "sample_number: $data_cfg->{sample_num}\n";

    for my $s (sort keys %{$data_cfg->{rawdata}}) {
        print "${s}_fq1: $data_cfg->{rawdata}{$s}{fq1}\n";
        print "${s}_fq2: $data_cfg->{rawdata}{$s}{fq2}\n" if (exists $data_cfg{rawdata}{$s}{fq2});
    }
}


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
#############################################################
sub GetTime {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
#############################################################
#&run_or_die($cmd);
sub run_or_die()
{
	my ($cmd) = @_ ;
	&show_log($cmd);
	my $flag = system($cmd) ;
	if ($flag != 0){
		&show_log("Error: command fail: $cmd");
		exit(1);
	}
	&show_log("done.");
	return ;
}
#############################################################
sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}
#############################################################
sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

#############################################################


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:Liuxs <Liuxs\@biomarker.com.cn> 
Description:
Usage:
  Options:
  -cfg <file> 
 
  -h         Help

USAGE
	print $usage;
	exit;
}




