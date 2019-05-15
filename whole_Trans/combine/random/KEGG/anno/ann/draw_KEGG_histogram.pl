#! /usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME = time();
my $version = "1.0.0";
##########################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($input, $odir, $prefix, $kegg_map_file,$png,$top_show);

GetOptions(
    "help|?" =>\&USAGE,
    "ipf:s" => \$input,
    "opd:s" => \$odir,
    "prf:s" => \$prefix,
    "map:s" => \$kegg_map_file,
    "png:s" => \$png,
    "top:i" => \$top_show,
) or &USAGE;
&USAGE unless ($input and $prefix);

my $wd = `pwd`; chomp($wd);
$odir ||= $wd;
$png ||= 'Y';
$top_show ||= 50;
$kegg_map_file ||= "$Bin/br08901.keg";
system "mkdir -p $odir" unless (-d $odir);

&USAGE unless ($png eq 'Y' or $png eq 'N');
print STDOUT "\n[".&GetTime($BEGIN_TIME)."] $Script start ...\n";
# --------------------------------------------------------------------
# parse kegg map file 
# --------------------------------------------------------------------
my %kegg_map_tree;

open (MAP, $kegg_map_file) or die $!;
$/ = "\n!\n"; <MAP>;

while (<MAP>) {
    chomp;
    next if (/^\s*$/ or /^#/);
    my @B_hierarchy = (split /\nB\s+/);
    my ($A_hierarchy_name) = $B_hierarchy[0] =~/A<b>(.+)<\/b>/;

    for (my $b=1; $b<@B_hierarchy; $b++) {
        my @C_hierarchy = (split /\nC\s+/,$B_hierarchy[$b]);
        my $B_hierarchy_name = $C_hierarchy[0];

        for (my $c=1; $c<@C_hierarchy; $c++) {
            my $C_hierarchy_name = (split /\s+/,$C_hierarchy[$c],2)[1];
            $kegg_map_tree{$A_hierarchy_name}{$B_hierarchy_name}{$C_hierarchy_name} = 1;
        }
    }

}

close MAP;

$/ = "\n";
#print Dumper(\%kegg_map_tree);
# --------------------------------------------------------------------
# parse DEG KEGG pathway enrichment result
# --------------------------------------------------------------------
my (%annotated_gene, %kegg_ortholog_groups);
my %annotation;
my (%annotated_gene_saved, %kegg_ortholog_groups_saved);
my ($pathway_number, $gene_number);

open (ANN, "<$input") or die $!;

while (<ANN>) {
    chomp;
    next if (/^\s*$/);

    if (/^#/) {
        $gene_number = (split /\t/)[1] if ($. == 1);
        die "ERROR: input file is illegal!\n" if ($gene_number =~/[^\d]/);
        next;
    }

    my @col = (split /\t+/);
    die "ERROR: input file is illegal, please check!\n" if (@col != 8);
    my ($pathway_id, $gene_num, $Gs, $Ks)=(split /\t+/)[0,2,-2,-1];

    my %tmp_Ks;

    foreach my $k (split /\+/,$Ks) {
        $tmp_Ks{$k} = 1;
    }

    $annotated_gene{$pathway_id} = $gene_num;
    $kegg_ortholog_groups{$pathway_id} = scalar keys %tmp_Ks;
    $annotation{$pathway_id} = "$Gs\t$Ks";
}

close ANN;

## reserve top 20% or top 50
my @annotated_pathway = sort { $annotated_gene{$b} <=> $annotated_gene{$a} } keys %annotated_gene;
$top_show = (@annotated_pathway >= $top_show) ? $top_show : @annotated_pathway;

for (my $i=0; $i<$top_show; $i++) {
#for (my $i=0; $i<int(@annotated_pathway*0.20); $i++) {
    $annotated_gene_saved{$annotated_pathway[$i]} = $annotated_gene{$annotated_pathway[$i]};
    $kegg_ortholog_groups_saved{$annotated_pathway[$i]} = $kegg_ortholog_groups{$annotated_pathway[$i]};
}

$pathway_number = keys %annotated_gene_saved;

#print Dumper(\%annotated_gene_saved);
# --------------------------------------------------------------------
# DEG KEGG enrichmented pathway tree stat
# --------------------------------------------------------------------
my %hierarchy_g;
my %hierarchy_k;

open (STAT,">$odir/$prefix.tree.stat") or die $!;
print STAT "#KEGG Category\tUnigene Count(K Number Count)\n";

foreach my $typeA (keys %kegg_map_tree) {
    my $A_print = '';

	foreach my $typeB (keys %{$kegg_map_tree{$typeA}}) {
        my $B_print = '';

		foreach my $typeC (keys %{$kegg_map_tree{$typeA}{$typeB}}) {
            next unless (exists $annotation{$typeC});
            my ($Gs,$Ks) = (split /\t/,$annotation{$typeC})[0,1];
            my @Gs = (split /;/,$Gs);
            my @Ks = (split /\+/,$Ks);

            for my $g (@Gs) {
                $hierarchy_g{$typeA}{$g} = 1;
                $hierarchy_g{$typeB}{$g} = 1;
                $hierarchy_g{$typeC}{$g} = 1;
            }

            for my $k (@Ks) {
                $hierarchy_k{$typeA}{$k} = 1;
                $hierarchy_k{$typeB}{$k} = 1;
                $hierarchy_k{$typeC}{$k} = 1;
            }

            my $C_gene_count = keys %{$hierarchy_g{$typeC}};
            my $C_knum_count = keys %{$hierarchy_k{$typeC}};
            $B_print .= "\t\t$typeC\t$C_gene_count \($C_knum_count\)\n" if ($C_gene_count != 0);
		}

        my $B_gene_count = keys %{$hierarchy_g{$typeB}};
        my $B_knum_count = keys %{$hierarchy_k{$typeB}};
        $A_print .= "\t$typeB\t$B_gene_count \($B_knum_count\)\n$B_print"  if ($B_gene_count != 0);
	}

    my $A_gene_count = keys %{$hierarchy_g{$typeA}};
    my $A_knum_count = keys %{$hierarchy_k{$typeA}};
    print STAT "$typeA\t$A_gene_count \($A_knum_count\)\n$A_print" if ($A_gene_count != 0);
}

close STAT;
# --------------------------------------------------------------------
# plot svg histogram
# --------------------------------------------------------------------

my $width = 800;
my $hight = 600;

my $top_margin = 40;
my $bottom_margin = 60;

my $path_hight = ($hight-$top_margin-$bottom_margin)/$pathway_number;
my $path_box_hight = $path_hight*0.80;
my $X_start = $width*3/8;
my $Y_start = 40;

my $X_start1 = $width*7/8;

my $lab_x=$width*5/8;
my $lab_y=$hight-30;

my %colors=(
    'Metabolism'=>'green',
    'Genetic Information Processing'=>'pink',
    'Environmental Information Processing'=>'purple',
    'Cellular Processes'=>'yellow',
    'Organismal Systems'=>'blue',
    'Human Diseases'=>'red',
    'Drug Development'=>'gray'
);

open (SVG,">$odir/$prefix.svg") or die $!;

print SVG '<?xml version="1.0" standalone="no"?>'."\n";
print SVG '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"'."\n";
print SVG '"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">'."\n";
print SVG "<svg width=\"$width\" height=\"$hight\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">"."\n";

foreach my $typeA (keys %kegg_map_tree) {
	my $Y_st=$Y_start;
	my $color=$colors{$typeA};
    my $plot_A = 0;

	foreach my $typeB (sort keys %{$kegg_map_tree{$typeA}}) {
		foreach my $typeC(sort keys %{$kegg_map_tree{$typeA}{$typeB}}) {
			my $path_id=$typeC;
			next unless exists $annotated_gene_saved{$path_id};
			my $path_gene_num = $annotated_gene_saved{$path_id};
            my $path_gene_per = sprintf("%.2f",$path_gene_num*100/$gene_number)."%";
			my $path_k_num=$kegg_ortholog_groups{$path_id};
			my $path_width=$path_gene_num/$gene_number*$width/2*5;

            my $path_id_x = $X_start+$path_width+10;
            my $path_id_y = $Y_start+$path_hight*2/3;
            my $path_num_x = $X_start-20;
            my $path_num_y = $Y_start+$path_hight*2/3;

            if ($path_gene_num/$gene_number >= 0.2) {
                $path_width = 0.2*$width/2*5;
                $path_id_x = $X_start+$path_width-60;
                $path_gene_num = "$path_gene_num \($path_gene_per\)";
            }

            print SVG "<rect x=\"$X_start\" y=\"$Y_start\" width=\"$path_width\" height=\"$path_box_hight\" style=\"fill:$color;stroke-width:1;stroke:black\"/>\n";
            print SVG "<text x=\"$path_id_x\" y=\"$path_id_y\" font-size=\"9\" font-family=\"Arial\" style=\"text-anchor:start\">$path_gene_num</text>\n";
            print SVG "<text x=\"$path_num_x\" y=\"$path_num_y\" font-size=\"9\" font-family=\"Arial\" style=\"text-anchor:end\">$path_id</text>\n";

			$Y_start += $path_hight;
            $plot_A = 1;
		}
	}

    my $X_st=$X_start1+10;
	print SVG "<line x1=\"$X_st\" y1=\"$Y_st\" x2=\"$X_st\" y2=\"$Y_start\" style=\"stroke:$color;stroke-width:2\"/>\n";

	my $category_bar_x=$X_st-10;
	my $category_bar_y=($Y_st+$Y_start)/2+5;
	print SVG "<text x=\"$category_bar_x\" y=\"$category_bar_y\" font-size=\"10\" font-family=\"Arial\" style=\"text-anchor:end\">$typeA</text>\n" if ($plot_A == 1);
}

print SVG "<line x1=\"$X_start\" y1=\"$Y_start\" x2=\"$X_start1\" y2=\"$Y_start\" style=\"stroke:black;stroke-width:2\"/>\n";

foreach my $i (0,1,2,3,4) {
	my $x_start=$X_start+$width/2*0.75*$i/3;
	my $y_end=$Y_start+5;
	print SVG "<line x1=\"$x_start\" y1=\"$Y_start\" x2=\"$x_start\" y2=\"$y_end\" style=\"stroke:black;stroke-width:1\"/>\n";
	my $test=$i*5;
	my $text_y=$y_end+10;
	print SVG "<text x=\"$x_start\" y=\"$text_y\" font-size=\"10\" font-family=\"Arial\" style=\"text-anchor:middle\" >$test%</text>\n";
}

print SVG "<text x=\"$lab_x\" y=\"$lab_y\" font-size=\"10\" font-family=\"Arial\" style=\"text-anchor:middle\" >Annotated Genes</text>\n";
print SVG '</svg>';

close SVG;

# --------------------------------------------------------------------
# convert svg to png
# --------------------------------------------------------------------
#my $svg_converter="/share/nas1/litc/tools/svg2xxx_release/svg2xxx";
#my $svg_converter="/share/nas1/yangxh/bin/tools/distributing_svg_4.74/svg2xxx_release/svg2xxx"; # 2014-12-11 ~ 
my $svg_converter="/share/nas2/genome/bmksoft/tool/svg2xxx/v1.0/svg2xxx";
if ($png eq 'Y') {
    chdir $odir;
    system "$svg_converter $prefix.svg ./";
}

##########################################################################################
print STDOUT "\n[".&GetTime(time())."] $Script done. Total elapsed time: ".(time()-$BEGIN_TIME)."s\n";
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
####################################################################################################
sub GetTime {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

####################################################################################################
sub USAGE {
    my $usage =<<"USAGE";
#-------------------------------------------------
Program: $Script
Version: $version
Contact: Simon Young <simonyoung8824\@gmail.com>
   Data: 2014-09-10
Fuction: draw KEGG kegg pathway enrichment histogram.
  Usage: 
    --ipf  <STR>    input file, DEG KEGG pathway enrichment result.
    --opd  <STR>    output directory.                               [./]
    --prf  <STR>    prefix of output files.
    --map  <STR>    KEGG pathway maps file.                         [\$Bin/br08901.keg]
    --top  <INT>    number of pathway to show.                      [50]
    --png  <STR>    convert svg to png. 'Y'|'N'                     ['Y']

Example:
    perl $Script --ipf T4_vs_T1.KEGG.xls --prf T4_vs_T1.KEGG --top 40
    perl $Script --ipf Anno_enrichment/pathway/kegg_enrichment/T4_vs_T1.KEGG.xls --opd Anno_enrichment/pathway/kegg_enrichment/ --prf T4_vs_T1.KEGG

#-------------------------------------------------
USAGE
    print $usage;
    exit;
}
