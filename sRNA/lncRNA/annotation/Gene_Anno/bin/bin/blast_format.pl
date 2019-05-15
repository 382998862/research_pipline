#!/usr/bin/perl -w
use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0";
my @Times = localtime();
use newPerlBase;
my %config=%{readconf("$Bin/../../db_file.cfg")}; 

#######################################################################################
my $Time_Start = sub_format_datetime(localtime(time()));
print STDOUT "$Script start at:[$Time_Start]\n";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($findir,$foutdir,$fshelldir,$fmiddledir,$Nr,$Nt,$Swissprot,$TrEMBL,$Cog,$Kog,$Kegg,$Eval,$queues,$eggNOG);
GetOptions(
				"help|?" =>\&USAGE,
				"id:s"=>\$findir,
				"middir:s"=>\$fmiddledir,
				"od:s"=>\$foutdir,
				"shdir:s"=>\$fshelldir,
				"nr"=>\$Nr,
				"nt"=>\$Nt,
				"swissprot"=>\$Swissprot,
				"eggNOG"=>\$eggNOG,
				"trembl"=>\$TrEMBL,
				"cog"=>\$Cog,
				"Kog"=>\$Kog,
				"kegg"=>\$Kegg,
				"eval"=>\$Eval,
				"queue:s"=>\$queues,
				) or &USAGE;
&USAGE unless ($findir and $foutdir and $fmiddledir);
$findir = Cwd::realpath($findir);
$foutdir = Cwd::realpath($foutdir);
$fmiddledir = Cwd::realpath($fmiddledir);

if (defined $fshelldir) {
	$fshelldir = Cwd::realpath($fshelldir);
	mkdir($fshelldir) if (!-d $fshelldir);
}else{
	$fshelldir = getcwd;
}
$queues||="middle.q";
$Eval=$Eval || "1e-5";
$findir =~ s/\/$//;
my $Cpu = 50;
my $notename = ` hostname `;
chomp $notename;

my $dirname = dirname($findir);

mkdir("$dirname/work_sh") if (!-d "$dirname/work_sh");
mkdir($foutdir) if (!-d $foutdir);
mkdir($fmiddledir) if (!-d $fmiddledir);
my $blast_parser = "$Bin/blast_parser.pl";
my $seq_file_name;

####################
my @blastm0files=glob("$findir/*.blast");
##get the the best hit
if (@blastm0files) {
    $seq_file_name = basename("$blastm0files[0]");
    open TAB1,">$fshelldir/Covertm0_2Tabbest.sh" || die "$!";
    open TAB2,">$fshelldir/Covertm0_2Tab.sh" || die $!;

    foreach my $blastfile (sort @blastm0files) {
        my $basename = basename($blastfile);
        print TAB1 "perl $blast_parser -nohead -eval $Eval -tophit 1 -m 0 -topmatch 1 $blastfile > $fmiddledir/$basename.tab.best \n" ;
        print TAB2 "perl $blast_parser -nohead -eval $Eval -tophit 50 -m 0 -topmatch 1 $blastfile > $fmiddledir/$basename.tab \n" ;
    }
    close TAB1;
    close TAB2;
	`sh $fshelldir/Covertm0_2Tabbest.sh`;
	`sh $fshelldir/Covertm0_2Tab.sh`;
    #&qsubOrDie("$fshelldir/Covertm0_2Tabbest.sh",$queues,$Cpu,"1G");
    #&qsubOrDie("$fshelldir/Covertm0_2Tab.sh",$queues,$Cpu,"1G");
}

####################
my @blastm7files=glob("$findir/*.blast.xml");
##get the the best hit
if (@blastm7files) {
    $seq_file_name = basename("$blastm7files[0]");
    open TAB3,">$fshelldir/Covertm7_2Tabbest.sh" || die "$!";
    open TAB4,">$fshelldir/Covertm7_2Tab.sh" || die $!;

    foreach my $blastfile (sort @blastm7files) {
        my $outprefix=basename($blastfile);
        $outprefix=~s/\.xml$//;
        print TAB3 "perl $blast_parser -nohead -eval $Eval -tophit 1 -topmatch 1 -m 7 $blastfile > $fmiddledir/$outprefix.tab.best \n" ;
        print TAB4 "perl $blast_parser -nohead -eval $Eval -tophit 50 -topmatch 1 -m 7 $blastfile > $fmiddledir/$outprefix.tab \n" ;
    }
    close TAB3;
    close TAB4;

	`sh $fshelldir/Covertm7_2Tabbest.sh`;
	`sh $fshelldir/Covertm7_2Tab.sh`;
    #&qsubOrDie("$fshelldir/Covertm7_2Tabbest.sh",$queues,$Cpu,"1G");
    #&qsubOrDie("$fshelldir/Covertm7_2Tab.sh",$queues,$Cpu,"1G");
}

##################### Cat Tab OUT
$seq_file_name =~ s/\.(\d+?)\.fa\.\S+//;

if (defined $Nr) {
	system("cat $fmiddledir/*.nr.blast.tab.best >$foutdir/$seq_file_name.nr.blast.tab.best");
	system("cat $fmiddledir/*.nr.blast.tab >$foutdir/$seq_file_name.nr.blast.tab");
}
if(defined $Nt){
	system("cat $fmiddledir/*.nt.blast.tab.best >$foutdir/$seq_file_name.nt.blast.tab.best");
	system("cat $fmiddledir/*.nt.blast.tab >$foutdir/$seq_file_name.nt.blast.tab");
}
if(defined $Swissprot){
	system("cat $fmiddledir/*.Swissprot.blast.tab.best >$foutdir/$seq_file_name.Swissprot.blast.tab.best");
	system("cat $fmiddledir/*.Swissprot.blast.tab >$foutdir/$seq_file_name.Swissprot.blast.tab");
}
if(defined $TrEMBL ){
	system("cat $fmiddledir/*.TrEMBL.blast.tab.best >$foutdir/$seq_file_name.TrEMBL.blast.tab.best");
	system("cat $fmiddledir/*.TrEMBL.blast.tab >$foutdir/$seq_file_name.TrEMBL.blast.tab");
}
if(defined $Cog ){
	system("cat $fmiddledir/*.Cog.blast.tab.best >$foutdir/$seq_file_name.Cog.blast.tab.best");
	system("cat $fmiddledir/*.Cog.blast.tab >$foutdir/$seq_file_name.Cog.blast.tab");
}
if(defined $Kog ){
	system("cat $fmiddledir/*.Kog.blast.tab.best >$foutdir/$seq_file_name.Kog.blast.tab.best");
	system("cat $fmiddledir/*.Kog.blast.tab >$foutdir/$seq_file_name.Kog.blast.tab");
}
if(defined $Kegg ){
	system("cat $fmiddledir/*.Kegg.blast.tab.best >$foutdir/$seq_file_name.Kegg.blast.tab.best");
	system("cat $fmiddledir/*.Kegg.blast.tab >$foutdir/$seq_file_name.Kegg.blast.tab");
}
if(defined $eggNOG){
	system("cat $fmiddledir/*.eggNOG.blast.tab.best >$foutdir/$seq_file_name.eggNOG.blast.tab.best");
	system("cat $fmiddledir/*.eggNOG.blast.tab >$foutdir/$seq_file_name.eggNOG.blast.tab");
}




#######################################################################################
my $Time_End = sub_format_datetime(localtime(time()));
print STDOUT "\n$Script Done at: [$Time_End]\t\tTotal elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub sub_format_datetime {   #Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub USAGE {
	my $usage=<<"USAGE";
	Program:$Script
	Version:$version	[2013/4/1]
	Contact:Sun Huaiyu <sunhy\@biomarker.com.cn>
	Descriptions: Get fixed format file according to blast result file.
	Options:
		-id	<dir>	input blast result dir                              required
		-od	<dir>	output dir                                          required
		-middir	<dir>	dir where put the middle result file            required
		-shdir	<dir>	dir where put shell script  default[./]         optional
		-nr              search against Nr annotation                   optional
		-nt              search against Nt annotation                   optional
		-swissprot       search against SwissProt annotation            optional
		-eggNOG          search against eggNOG annotation               optional
		-trembl          search against TrEMBL annotation               optional
		-cog             search against Cog annotation                  optional
		-kog             search against KOG annotation                  optional
		-kegg            search against Kegg annotation                 optional
		-eval       float or exponent,to filter the alignments which worse than the E-value cutoff,default 1e-5,optional
		-h	Help

USAGE
	print $usage;
	exit;
}
