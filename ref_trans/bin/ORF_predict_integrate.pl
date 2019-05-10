#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use threads;
use FindBin qw($Bin $Script);
use Config::General;
use File::Basename qw(basename dirname);
use Cwd qw(abs_path getcwd);
use newPerlBase;

#######################################################################################
my $BEGIN_TIME=time();
my $Title="Ref_Trans_v3.0";
my $version="v1.0.0";
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fasta_file,$odir,$pfam,$known_peptide_file,$top_strand_only,$prefix,$verbose,$new_cds_predict);
GetOptions(
				"fa:s"=>\$fasta_file,
				"od:s"=>\$odir,
				"pfam:s"=>\$pfam,
				"known_pep:s"=>\$known_peptide_file,
				"ss:s"=>\$top_strand_only,
				"prefix:s"=>\$prefix,
				"verb"=>\$verbose,
				"help|h" =>\&USAGE,
				"new_cds_predict:s"=>\$new_cds_predict,
				) or &USAGE;
&USAGE unless ($fasta_file and $odir and $known_peptide_file and $pfam);

########
#my $TransDecoder = ${&readconf("$Bin/../config/db_file.cfg")}{TransDecoder};
my $TransDecoder ="$Bin/annotation/Gene_Anno/bin/TransDecoder/TransDecoder.pl";
###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";  

$known_peptide_file=&ABSOLUTE_DIR($known_peptide_file);
exit unless (-s $known_peptide_file) ;
die "ERROR: Cannot locate pfam database at: $pfam.\n" unless (-s $pfam);
system "mkidr -p $odir" unless(-d $odir);
$odir=&ABSOLUTE_DIR($odir);
my $cds_predict_dir="$odir/CDS_Predict";
system "mkdir $cds_predict_dir" unless -d $cds_predict_dir;
$prefix ||="Unknown";

my @fasta=split /;/,$fasta_file;
my $nt_seq_file;	
if (@fasta>1) {
	my $cat_fa;
	for (@fasta) {
		$_=&ABSOLUTE_DIR($_);
		$cat_fa.="$_ ";
	}
	$nt_seq_file="$cds_predict_dir/longest_transcript.fa";
	system "cat $cat_fa >$nt_seq_file";
} else {
	$nt_seq_file=&ABSOLUTE_DIR($fasta_file);
}
my $cmd;
if (defined $new_cds_predict and -f $new_cds_predict) {
    $cmd= "cp $new_cds_predict  $cds_predict_dir/$prefix.best_candidates.final.pep.fa";
	print "$cmd\n";
	system "$cmd";
	
}
else{
	$cmd= "perl $TransDecoder -trans $nt_seq_file -pref $prefix -od $cds_predict_dir/ -pfam $pfam ";
	$cmd.="-ss "if (defined $top_strand_only) ;
	print "$cmd\n";
	system "$cmd";
}
##################modification of fasta file header format and extract the longest peptide(protein)
##################fasta header:gene_id\t$transcript_id; only this format is acceptted;
#########get transcript id
my %gene_trans;
open IN,"<",$nt_seq_file or die $!;
$/='>';
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my $head=(split/\n+/,$_,2)[0];
	my ($gene,$transcript)=(split/\s+/,$head);
	$gene_trans{$gene}=$transcript;
}
close IN;

#########record peptide(protein) seq
my %pep_seq;
open IN,"<","$cds_predict_dir/$prefix.best_candidates.final.pep.fa" or die $!;
while(<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my ($head,$seq)=split/\n+/,$_,2;
	my $id=(split/\s+/,$head)[0];
	my $gene_id=(split /\|/,$id)[0];
	$seq=~s/\s+//g;
	my $length=length($seq);
	if (!defined $pep_seq{$gene_id}) {
		$pep_seq{$gene_id}=$seq;
	} elsif (defined $pep_seq{$gene_id} and $length> length($pep_seq{$gene_id})) {
		$pep_seq{$gene_id}=$seq;
	}
}
$/='\n';
close IN;

#########output peptide(protein) seq
my $new_peptide_file="$cds_predict_dir/$prefix.new.longest_transcript.pep";
open OUT,">",$new_peptide_file or die $!;
for my $id(sort keys %pep_seq) {
	print OUT ">$id\t$gene_trans{$id}\n$pep_seq{$id}\n";
}
close OUT;

######integrate
system "cat $known_peptide_file $new_peptide_file > $cds_predict_dir/${prefix}.longest_transcript.pep.fa";
system "mv $cds_predict_dir/${prefix}.longest_transcript.pep.fa $odir";

#system "rm -rf $cds_predict_dir" unless(defined $verbose);

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";
&Runtime($BEGIN_TIME);

###########subs
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

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
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

sub USAGE
{
	print <<"	Usage End.";
	Description:
		Function : use Getorf predict transcript cDNA and the Pep sequence;
		Version  : $version.
		Writer   : liut <liut\@biomarker.com.cn>
		Usage    :
		-fa
		   unigene seq fa, seperated by comma when two or more files are provided, required;
		-od
		   output directory, required;
		-known_pep
		   known peptide(protein) file, which will be integrated with predicted unigene peptide(protein) file, required; 
		-pfam
			path to Pfam database, *.hmm, optional;
		-ss
		   strand-specific, only analyzes top strand, optional;
		-prefix
		   prefix for finally integrated peptide(protein) file (default: "longest_protein"), optional;
		-verb
			maintain intermediate files and directories,optional;
		-help
		    Help document
	Usage End.

	exit;
}
