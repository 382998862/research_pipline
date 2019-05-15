#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use threads;
use FindBin qw($Bin $Script);
use Config::General;
use File::Basename qw(basename dirname);
#use Cwd qw(abs_path getcwd);
use newPerlBase;

#######################################################################################
my $BEGIN_TIME=time();
my $Title="Ref_Trans_v3.0";
my $version="v1.0.0";
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fasta_file,$odir,$pfam,$verbose,$cfg,$prefix,$newpep);
$odir="./";
GetOptions(
				"fa:s"=>\$fasta_file,
				"newpep:s"=>\$newpep,
				"od:s"=>\$odir,
				"pfam:s"=>\$pfam,
				"cfg:s"=>\$cfg,
				"verb"=>\$verbose,
				"help|h" =>\&USAGE,
				) or &USAGE;
&USAGE unless ( $odir and $cfg);

###############Time
	my $TransDecoder="$Bin/TransDecoder/TransDecoder.pl";
$pfam ||= "/share/biocloud-compute/anno_database/Pfam/201703/Pfam-A.hmm";
$prefix ||="Unknown";
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";  
system "mkdir -p $odir" unless(-d $odir);
$odir=&ABSOLUTE_DIR($odir);
$cfg=&ABSOLUTE_DIR($cfg);
exit unless (-s $cfg) ;
my %CFG;
&LOAD_PARA($cfg,\%CFG);
exit unless (exists $CFG{'Known_anno'});
my $known_peptide_file;
$known_peptide_file="$CFG{Known_anno}/Known.longest_transcript.pep.fa";
if (!-f $known_peptide_file) {
    $prefix="Known";
	system "mkdir -p $odir/Known_CDS" unless(-d "$odir/Known_CDS");
	&Trans_Decoder($CFG{'Known_unigene'},"$odir/Known_CDS",$CFG{'Lib_type'},$prefix);
	$known_peptide_file="$odir/Known_CDS/$prefix.longest_transcript.pep";
}
my $new_peptide_file;
if (defined $newpep && defined $fasta_file){
	$prefix="New";
	system "mkdir -p $odir/CDS_Predict" unless(-d "$odir/CDS_Predict");
	&deal_pep_fa($newpep,"$odir/CDS_Predict/$prefix.longest_transcript.pep",$fasta_file);
	$new_peptide_file="$odir/CDS_Predict/$prefix.longest_transcript.pep";
}elsif (!defined $newpep && defined $fasta_file ) {
    $prefix="New";
	system "mkdir -p $odir/CDS_Predict" unless(-d "$odir/CDS_Predict");
	&Trans_Decoder($fasta_file,"$odir/CDS_Predict",$CFG{'Lib_type'},$prefix);
	$new_peptide_file="$odir/CDS_Predict/$prefix.longest_transcript.pep";
}elsif(!defined $fasta_file ){
	$new_peptide_file="";
}

######integrate
system "cat $known_peptide_file $new_peptide_file > $odir/$CFG{Project_key}.longest_transcript.pep.fa";

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";
&Runtime($BEGIN_TIME);

###########subs
sub deal_pep_fa{
    my ($pep,$out,$fa)=@_;
	my %gene_trans;
	open IN,"<",$fa or die $!;
	$/='>';
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		my $head=(split/\n+/,$_,2)[0];
		my ($gene,$transcript)=(split/\s+/,$head);
		$gene_trans{$gene}=$transcript;
	}
	close IN;
    my %pep_seq;
	open IN,"<",$pep or die $!;
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

	open OUT,">",$out or die $!;
	for my $id(sort keys %pep_seq) {
		print OUT ">$id $gene_trans{$id}\n$pep_seq{$id}\n" if (exists $pep_seq{$id});
		print OUT ">$id $gene_trans{$id}\n--\n" unless (exists $pep_seq{$id});
	}
	close OUT;
}
sub  Trans_Decoder{
	my ($fa,$od,$lib_type,$pre)=@_;
	my $cmd;
	$cmd= "$TransDecoder -trans $fa -pref $pre -od $od/ -pfam $pfam -ss";
	#$cmd.="-ss "if ( $lib_type ne "fr-unstranded") ;
	print "$cmd\n";
	system "$cmd";
	
	##################modification of fasta file header format and extract the longest peptide(protein)
	##################fasta header:gene_id\t$transcript_id; only this format is acceptted;
	#########get transcript id
	my %gene_trans;
	open IN,"<",$fa or die $!;
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
	open IN,"<","$od/$pre.best_candidates.final.pep.fa" or die $!;
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
	my $peptide_file="$od/$pre.longest_transcript.pep";
	open OUT,">",$peptide_file or die $!;
	for my $id(sort keys %pep_seq) {
		print OUT ">$id $gene_trans{$id}\n$pep_seq{$id}\n";
	}
	close OUT;
}
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
sub LOAD_PARA {
	my $para_file= shift;
	my $para= shift;
    my %db_node;
    my $db_idx = 0;
	my $error_status = 0;
	open IN,$para_file || die "fail open: $para_file";
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r$//;
		next if (/^$/ or /^\#/) ;
		my ($para_key,$para_value) = split(/\s+/,$_);
        #$para_value =~s/\/lustre|\/share/$db_node{$random}/ unless ($para_key eq 'mRNA');
        $para->{$para_key} = $para_value;
	}
	close IN;
	die "\nExit due to error Undefined parameter\n" if ($error_status) ;
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
		   nexgene seq fa, seperated by comma when two or more files are provided,optional ;
		-newpep
		    nexgene  pep fa ,optional ;
		-od
		   output directory, required;
		-cfg
		   detail.cfg file, required; 
		-pfam
			path to Pfam database, *.hmm, optional;
		-verb
			maintain intermediate files and directories,optional;
		-help
		    Help document
	Usage End.

	exit;
}
