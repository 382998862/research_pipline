# Writer:      huangls

my $ver="1.0.0";

use strict;
use Cwd qw(abs_path);
use Getopt::Long;
my $BEGIN=time();
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $programe_dir=basename($0);
my $path=dirname($0);

############ 流程名称，版本，及固定配置文件
use lib "/share/nas2/genome/bmksoft/tool/newPerlBase/v1.0";
use newPerlBase;


######################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作


our %opts;
GetOptions(\%opts,"anno=s","key=s","od=s","peakann=s","h" );

sub help{
print <<"Usage End.";

contact:huangls

Usage:
-peakann          peak gene list  files dir                           must be given;
-anno             Annotation dir (All genes Result dir for enrichment  optional;
                  must be given if defined KEGG & GO enrichment analysis;
                  should contain *.ko; *.path; *.GO_tree.xls;
-key
-od               Out dir                                         must be given;
-h                help document
Usage End.
exit;
}

if(!defined($opts{od}) || defined($opts{h})||!defined($opts{key})||!defined($opts{peakann})){
	&help();
	exit;
}


###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\n[$Time_Start] $Script start ... \n\n";

################
my $idir=abs_path($opts{idir});
mkdir $opts{od} unless (-d  $opts{od});
my $outdir=abs_path($opts{od});

my $enrichment_dir;
my $all_go_file;
my $gene=abs_path($opts{peakann});
$enrichment_dir=abs_path($opts{anno});
$all_go_file=(glob "$enrichment_dir/*.GO.list.txt")[0];

mkdir "$outdir/$opts{key}" unless (-d "$outdir/$opts{key}");

# kegg go enrich
my$cmd = "perl $Bin/ann/KeggGo_enrich_map_web.pl -d $gene ";
$cmd .= "-k $opts{key} -i $enrichment_dir -o $outdir/$opts{key} ";
runOrDie("$cmd");

# DEG anno
$cmd = "perl $Bin/DEG_Annotation.pl -DEG $gene";
$cmd .= " -id $enrichment_dir -od $outdir/$opts{key} -key $opts{key} ";
runOrDie("$cmd");
mkdir("$outdir/$opts{key}/Graph") unless(-d "$outdir/$opts{key}/Graph");

runOrDie("perl $Bin/ann/kegg_enrichment_plot.pl -enrich_file $outdir/$opts{key}/pathway/kegg_enrichment/$opts{key}.KEGG.stat -od $outdir/$opts{key}/Graph -key $opts{key}");
runOrDie("perl $Bin/ann/draw_GO_DAG_map.pl -All_GO $all_go_file -DEG_list $gene -od $outdir/$opts{key}/Graph -key $opts{key}");
runOrDie("perl $Bin/ann/draw_KEGG_histogram.pl --ipf $outdir/$opts{key}/pathway/kegg_enrichment/$opts{key}.KEGG.xls --opd $outdir/$opts{key}/pathway/kegg_enrichment/ --prf $opts{key}.KEGG");



###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $programe_dir Time :[$Time_End]\n\n";
&Runtime($BEGIN);

sub Runtime
{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total elapsed time: ${t}s\n";
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub ABSOLUTE_DIR
{ #$pavfile=abs_path($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;chomp $return;
	}
	else
	{
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}
