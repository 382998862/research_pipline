#!/usr/bin/env perl
use	strict;
use	Getopt::Long;
use	Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw/sum/;
my ($gtf,$out);
GetOptions(
			"help|?"=>\&USAGE,
			"gtf:s"=>\$gtf,
			"out:s"=>\$out,
			)or &USAGE;
&USAGE unless ($gtf);
$gtf=&ABSOLUTE_DIR($gtf);
####################
my %code=(
	"e" => "sense_lncRNA",
	"o" => "sense_lncRNA",
	"u" => "lincRNA",
	"i" => "intronic_lncRNA",
	"x" => "anti-sense_lncRNA",
);
open(IN,$gtf) or die $!;
my $i;
my %trans;
while (<IN>) {
    chomp;
	next if (/^#/ || /^$/);
	my ($trans_id,$type);
	my $line = $_;
	my @lines = split ("\t",$line);
	if ($lines[2]=~/transcript/){
		$i++;
		$line =~/transcript_id\s\"([^\"]+)\";/;
		$trans_id = $1;
		$line =~/class_code\s\"([^\"]+)\";/;
		$type =$1;
		my $class=$code{$type};
		$trans{$class}{$trans_id} = 1;
	}

}
close IN;
open(OUT,">$out") or die $!;
foreach my $lnc (sort keys %trans){
	my @id=keys %{$trans{$lnc}};
	my $num=@id;
	my $per=sprintf "%.2f",($num/$i);
	print OUT "#$lnc\t$num($per%)\n";
}
print OUT "#lncRNA_ID\tlncRNA_type\n";
foreach my $ln (sort keys %trans){
	my @ID=keys %{$trans{$ln}};
	foreach my $id (@ID){
		print OUT "$id\t$ln\n";
	}
}
close OUT;
###############################################################
sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	use Cwd 'abs_path';
    my ($in) = @_;
    my $ret = abs_path($in);
    if ( -e $ret ){
        return $ret;
    }else{
        warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
        exit;
    }
}

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:  $Script
 

Usage:  perl $Script -gtf lncRNA.gtf -out lncRNA_class.txt
Options:
		

USAGE
	print $usage;
	exit;
}
