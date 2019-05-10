#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd 'abs_path';
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;

my (@fIn, $gff, $fKey, $dOut, $log);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s{,}" =>\@fIn,
				"gff:s"  =>\$gff,
				"d:s"    =>\$dOut,
				) or &USAGE;
&USAGE unless (@fIn and $gff);

$dOut ||= "./";     
mkdir $dOut unless (-d $dOut);
$dOut   = abs_path($dOut);

my %snp = ();
my %snp_print = ();
my %sampleLabel = ();
my %density_stat; #	$density_stat{sample_id}{gene_id}=snp_num_per_kb;

#===============================================================
# Get Data
#===============================================================

foreach my $fIn (@fIn) {
	my ($sample) = basename($fIn) =~/(\w+)\.snp$/;
	$sampleLabel{$sample} = 1;
	
	open (IN,"$fIn") or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/) ;
		my ($chro, $pos, $refBase, @snpInfo) = split;

		my $pos_int = pos_unit($pos);
		push @{$snp{$chro}{$pos_int}{'pos_unit'}}, $pos;
		
		#@snpInfo: REF     Sample_Base     Qual    Depth   Genotype 
		$snp{$chro}{$pos}{'sam'}{$sample} = \@snpInfo;
		$snp{$chro}{$pos}{'processed'}    = 0;
		$snp{$chro}{$pos}{'refBase'}      = $refBase;
	}
	close (IN);
}

my @sample_list = sort {$a cmp $b} keys %sampleLabel;

#===============================================================
# Process
#===============================================================
open (OUT,">$dOut/sam.merge.snp") or die $!;
$|=1;
print OUT join("\t",
	'#Chro',
	'Pos',
	'Genic/Intergenic',
	'RefBase',
	(map {
		($_."_Base", $_."_Qual", $_."_Depth", $_."_Genotype");
	} @sample_list),
),"\n";

open (IN, $gff) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/) ;

	my ($chro, undef, $feture, $start, $end, undef, $strand, undef, $attribute) = split;
	next if ($feture ne "gene") ;

	my ($geneID) = $attribute =~/ID=([^;]+)/;
	($start, $end) = sort {$a <=> $b} ($start, $end);
    my $gene_len = $end - $start;
	my @pos = find_snps_in_gene($chro, $start, $end, \%snp);
    my %snp_support_sam;

#	print join("\t",@pos),"\n";
#	<STDIN>;

	if (@pos) {
        foreach my $snpPos (@pos) {
#            print OUT join("\t",
#                $chro,
#                $snpPos,
#                $geneID,
#                $snp{$chro}{$snpPos}{'refBase'},
#                (map {
#                    if (defined $snp{$chro}{$snpPos}{'sam'}{$_}){
#                        @{$snp{$chro}{$snpPos}{'sam'}{$_}};
#                    }else{
#                        (".", ".", ".", ".");
#                    }
#                } @sample_list),
#            ),"\n";
#--------------------------------------------------------------
            my $snp_inf = join("\t",$chro,$snpPos,);
            my $sample_inf = join("\t",
                $snp{$chro}{$snpPos}{'refBase'},
                (map {
                    if (defined $snp{$chro}{$snpPos}{'sam'}{$_}){
                        @{$snp{$chro}{$snpPos}{'sam'}{$_}};
                    }else{
                        (".", ".", ".", ".");
                    }
                } @sample_list),
            );
            $snp_print{$snp_inf}{$geneID} = $sample_inf;

            for my $s (@sample_list) {
                if (defined $snp{$chro}{$snpPos}{'sam'}{$s}){
                    $snp_support_sam{$s}++;
                }
            }

            $snp{$chro}{$snpPos}{'processed'} = 1;
        }

        for my $s (@sample_list) {
            if (exists $snp_support_sam{$s}) {
                $density_stat{$s}{$geneID} = int($snp_support_sam{$s}*1000/$gene_len+1);
            } else {
#                $density_stat{$s}{$geneID} = 0;
            }
        }

        $density_stat{'All'}{$geneID} = int(@pos*1000/$gene_len+1);
    } else { ## no snp found in this gene 
#        for my $s (@sample_list) {
#            $density_stat{$s}{$geneID} = 0;
#        }
#        $density_stat{'All'}{$geneID} = 0;

        next;
    }

}

# deal with lines like the following
#chr1    2742603 GRMZM2G085885   T       G       99      96      K       G       99      84      K       .       .       .       .       .       .       .       .
#chr1    2742603 GRMZM2G086242   T       G       99      96      K       G       99      84      K       .       .       .       .       .       .       .       .

for my $snp (sort keys %snp_print) {
    my $gene_ids = join(";",(sort keys %{$snp_print{$snp}}));
    my $sam_inf = $snp_print{$snp}{(sort keys %{$snp_print{$snp}})[0]};
    print OUT "$snp\t$gene_ids\t$sam_inf\n";
}

#
# output intergenic snps 
#
foreach my $chro (keys %snp) {
	foreach my $snpPos (keys %{$snp{$chro}}) {
		next unless exists $snp{$chro}{$snpPos}{'processed'};
		next if ($snp{$chro}{$snpPos}{'processed'} == 1) ;
		print OUT join("\t",
			$chro,
			$snpPos,
			'intergenic',
			$snp{$chro}{$snpPos}{'refBase'},
			(map {
				if (defined $snp{$chro}{$snpPos}{'sam'}{$_}){
					@{$snp{$chro}{$snpPos}{'sam'}{$_}};
				}else{
					(".", ".", ".", ".");
				}
			} @sample_list),
		),"\n";
	}
}

close (IN) ;
close (OUT) ;

#-------------------------------------------------
# snp density stat & plot 
#-------------------------------------------------
open (STAT,">$dOut/AllSample.SNP_density.stat") or die $!;
print STAT "#Sample\tInterval\tGeneNum\n";

for my $s (sort keys %density_stat) {
    my $print;
    my %cut;

    for my $g (keys %{$density_stat{$s}}) {
        if ($density_stat{$s}{$g}==0) {
            $cut{'0'}++;
        } elsif ($density_stat{$s}{$g}>0 and $density_stat{$s}{$g}<=8) {
            my $key = ($density_stat{$s}{$g}-1)."-"."$density_stat{$s}{$g}";
            $cut{$key}++;
        } else {
            $cut{'8~'}++;
        }
    }

    for my $c ('0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8~') {
        if (exists $cut{$c}) {
            $print .= "$s\t$c\t$cut{$c}\n";
        } else {
            $print .= "$s\t$c\t0\n";
        }
    }

    print STAT "$print";
}

close STAT;

# plot 
my $Rscript = "/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript";
#print "$Rscript $Bin/dodgedBar.r --infile $dOut/AllSample.SNP_density.stat --outfile $dOut/AllSample.SNP_density.png --x.col 2 --group.col 1 --y.col 3 --x.lab \"SNP Number per Kb\" --group.lab \"Sample\" --y.lab \"Number of Gene\" --title.lab \"SNP Density\" --legend.col 1 \n";
system "$Rscript $Bin/dodgedBar.r --infile $dOut/AllSample.SNP_density.stat --outfile $dOut/AllSample.SNP_density.png --x.col 2 --group.col 1 --y.col 3 --x.lab \"SNP Number per Kb\" --group.lab \"Sample\" --y.lab \"Number of Gene\" --title.lab \"SNP Density\" --legend.col 1 >/dev/null 2>&1 ";


# ==============================================================
# sub function
# ==============================================================
sub find_snps_in_gene {#
	my ($chro, $start, $end, $ref_snp) = @_;

	my $start_int = pos_unit($start);
	my $end_int   = pos_unit($end);
	my %tmp;
	for(my $i=$start_int; $i<=$end_int; $i++)
	{
		next unless exists $snp{$chro}{$i}{'pos_unit'};
		my @Pos = @{$snp{$chro}{$i}{'pos_unit'}};
		foreach (@Pos)
		{
			$tmp{$_}++;
		}
	}
	my @sortPos = sort {$a<=>$b} keys %tmp;
	
	my @res = ();

	foreach my $snpPos (@sortPos) {
		last if ($snpPos > $end) ;
		next if ($snpPos < $start) ;
		push @res, $snpPos;
	}

	return @res;
}

sub pos_unit {
	my ($pos) = @_;
	my $unit  = 1e4;
	my $pos_int = int($pos/$unit);
	return $pos_int;
}

sub USAGE {#
	my $usage=<<"USAGE";
Usage:
  -i	<files>	SNP files, sep by whitespace, snp format, forced
  -gff	<file>	GFF file, forced
  -d	<str>	Directory where output file produced,optional,default [./]
  -h		Help

USAGE
	print $usage;
	exit;
}
