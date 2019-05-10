#!/usr/bin/perl
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use Getopt::Long;
my %opts;

my ($fasq_file1,$fasq_file2,$bam_file,$tempath,$outputpath,$outputname,$datatype,$h);

GetOptions(

           "file1=s"=>\$fasq_file1,
           "file2=s"=>\$fasq_file2,
           "bamfile=s"=>\$bam_file,
           "outputpath=s"=>\$outputpath,
           "outputname=s"=>\$outputname,
           "datatype=s"=>\$datatype,
           "h"=>\&USAGE,

           )or &USAGE;

$bam_file = &ABSOLUTE_DIR($bam_file);
$outputpath = &ABSOLUTE_DIR($outputpath);
if ($datatype eq "PE"){
	configPE();
}elsif ($datatype eq "SE"){
	configSE();
}else {
	&USAGE
}

# ----------------------------------------------------------
# sub function
# ----------------------------------------------------------
sub configPE {
	open (IN,"$Bin/FusionPE.Linux.config")||die "open infile failed!\n";
	open (OUT,">$outputpath/FusionPE.config")||die "open or creat outfile failed!\n";
	while (<IN>) {
	$/ = "\n";
	if (/^\/s/){
	    $_ = "$bam_file";
	    print OUT $_,"\n";
	    }
	elsif (/^TempPath/){
	    $_ = "$outputpath/tmp";
	    print OUT "TempPath=".$_,"\n";
	}
	elsif (/^OutputPath/){
	    $_ = "$outputpath";
	    print OUT "OutputPath=".$_,"\n";
	}
	elsif (/^OutputName/){
		$_ ="$outputname";
	    print OUT "OutputName=".$_,"\n";
	}
	else {
		print OUT $_,"\n";
	}
}
close IN;
close OUT;
}
sub configSE {
	open (IN,"$Bin/SEconfig.config")||die "open infile failed!\n";
	open (OUT,">output.config")||die "open or creat outfile failed!\n";
	$num =0;
	while (<IN>) {
	$/ = "\n";
	if (/^\/I/){
		$num +=1;
		if ($num ==1){
	    $_= "$Bin/$fasq_file1";
		}
		else{
	    $_ = "$Bin/$fasq_file2";
		}
		print OUT $_,"\n";
	}
	elsif (/^Gzip/){
		$_ = "Gzip=False";
		print OUT $_,"\n";
	}
	elsif (/^TempPath/){
	    $_ = "$Bin/$tempath";
	    print OUT "TempPath=".$_,"\n";
	}
	elsif (/^OutputPath/){
	    $_ = "$Bin/$outputpath";
	    print OUT "OutputPath=".$_,"\n";
	}
	elsif (/^OutputName/){
		$_ ="$Bin/$outputname";
	    print OUT "OutputName=".$_,"\n";
	}
	else {
		print OUT $_,"\n";
	}
}
close IN;
close OUT;
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;
	$cur_dir =~ s/\n$//;
	my ($in)=@_;
	my $return="";

	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;
		$dir =~ s/\n$// ;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;
		$return =~ s/\n$// ;
	}
	else
	{
		warn "Warning just for file and dir [$in]\n";
		exit;
	}

	chdir $cur_dir;
	return $return;
}
##############################################################################################################
sub USAGE {
    my $usage=<<"__USAGE__";
#------------------------------------------------------------------------------------
Program: $Script
Version: $version
Contact: <gaom\@biomarker.com.cn>
   Data: 2016-04-07
Fuction: the script is used to get the FusionGene list by FusionMap.
  Usage:
        --file1    <STR>   input file, FASTA format, only requried in SE modle
        --file2    <STR>   input file, FASTA format, only requried in SE modle
        --bamfile          input file,BAM format,only requried in PE modle
        --tempath 		   the soft need a tempfile path during runing
        --outputpath       outfile, path, requried
        --outputname	   outputname, txt format, requried
        --datatype		   PE or SE, out variant, requried
        --h   			   h, help, choose
Example:
    perl $Script --file1 example1.fa --file2 example2.fa --bamfile example.bam --tempath tempath --outputpath outputpath --outputname outputname --datatype PE


#------------------------------------------------------------------------------------
__USAGE__
    print $usage;
    exit;
}
