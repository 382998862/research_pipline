#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my (@list, @label, $prefix, $odir);

GetOptions(
                "help|?" =>\&USAGE,
                "lst:s"=>\@list,
                "lab:s"=>\@label,
                "pre:s"=>\$prefix,
                "od:s" =>\$odir,
                ) or &USAGE;
&USAGE unless (@list and @label);

chomp(my $wd=`pwd`);
$prefix ||= 'circRNASets';
$odir   ||= $wd;
system "mkdir -p $odir" unless (-d $odir);
$odir = &ABSOLUTE_DIR($odir);

if (@list!=@label or @list < 2 or @list > 5) {
    print "ERROR: Too more or too less DEG sets.\n";
    &USAGE;
}

print STDOUT "\n[".&GetTime($BEGIN_TIME)."] ".basename($0)." start ...\n";
# ------------------------------------------------------------------
# plot 
# ------------------------------------------------------------------
my @color=("'dodgerblue'","'goldenrod1'","'darkorange1'","'seagreen3'","'orchid3'");
my %circRNA;
my $list_content;
my $label_content;
my $color_content;

for my $i (1..@list) {
    open (LIST, "<$list[$i-1]") or die $!;

    while (<LIST>) {
        chomp;
        next if (/^\s*$/ or /^#/);
        my $id = (split /\s+/)[0];
        $id = "'".$id."'";
        push @{$circRNA{$label[$i-1]}}, $id;
    }

    close LIST;
}

for my $i (0..@label-1) {
    $list_content.= "$label[$i] <- c(".(join ", ",@{$circRNA{$label[$i]}}).")\n";
    $label_content.= "$label[$i] = $label[$i], ";
    $color_content.= "$color[$i], ";
}

$list_content =~ s/\n$//;
$label_content =~ s/, $//;
$color_content =~ s/, $//;

my $more_opts = "";

if (@label == 5) {
    $more_opts.= "    cex = 0.5,\n";
    $more_opts.= "    cat.cex = 0.6,\n";
    $more_opts.= "    margin = 0.1,\n";
    $more_opts.= "    cat.dist = c(0.20, 0.25, 0.20, 0.20, 0.25),\n";
    $more_opts.= "    scaled = FALSE,\n";
}
elsif (@label == 4) {
    $more_opts.= "    cex = 0.6,\n";
    $more_opts.= "    cat.cex = 0.7,\n";
    $more_opts.= "    margin = 0.08,\n";
    $more_opts.= "    scaled = FALSE,\n";
}
elsif (@label == 3) {
    $more_opts.= "    cex = 0.7,\n";
    $more_opts.= "    cat.cex = 0.8,\n";
    $more_opts.= "    margin = 0.06,\n";
    $more_opts.= "    euler.d = FALSE,\n";
    $more_opts.= "    scaled = FALSE,\n";
}
elsif (@label == 2) {
    $more_opts.= "    cex = 0.8,\n";
    $more_opts.= "    cat.cex = 0.9,\n";
    $more_opts.= "    margin = 0.05,\n";
    $more_opts.= "    cat.pos = 180,\n";
    $more_opts.= "    euler.d = FALSE,\n";
    $more_opts.= "    scaled = FALSE,\n";
}

my $R_script = <<"SCRIPT";
#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript
# import library
library(grid)
library(VennDiagram)

#init deg lists
$list_content

#plot venn
venn.diagram(
    x = list($label_content),
    filename = "${prefix}_venn.png",
    imagetype="png",
    fill = c($color_content),
    height=800,
    width=800,
    resolution=200,
    units='px',
    lwd=0,
$more_opts
);

SCRIPT

open (RS, ">$odir/$prefix.venn.r") or die;
print RS $R_script;
close RS;

my $Rscript = "/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript";
system "cd $odir/ && $Rscript $prefix.venn.r";
#######################################################################################
print STDOUT "\n[".&GetTime(time())."] ".basename($0)." done. Total elapsed time: ".(time()-$BEGIN_TIME)."s\n\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
sub ABSOLUTE_DIR { #$pavfile=&ABSOLUTE_DIR($pavfile);
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

################################################################################################################
sub GetTime {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub USAGE {
    my $script = &basename($0);
    my $usage=<<"USAGE";
 ProgramName:   $script
     Version:   $version
     Contact:   Simon Young <yangxh\@biomarker.com.cn>/<simonyoung8824\@gmail.com>
Program Date:   2014-08-07
 Description:   This program is used to draw venn diagram of 2-5 differentially expressed gene sets
       Usage:
        --lst       DEG ids list, the first column is DEG ids
        --lab       label of DEG set

        Options:
        --pre       output file prefix, ['DEGsets']
        --od        output directory, [./]
        --help      help

     Example:
        perl $script --lst T1_vs_T2.DEG.final.xls --lab T1_vs_T2 --lst T1_vs_T3.DEG.final.xls --lab T1_vs_T3 --pre DEG2sets_venn
        perl $script --lst T1_vs_T2.DEG.final.xls --lab T1_vs_T2 --lst T1_vs_T3.DEG.final.xls --lab T1_vs_T3 --lst T1_vs_T4.DEG.final.xls --lab T1_vs_T4

USAGE
    print $usage;
    exit;
}
