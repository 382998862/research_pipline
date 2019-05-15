#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use threads;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
#my %config=%{readconf("$Bin/../../project.cfg")};

#my $version=$config{version};
my ($blast_out,$circ2disease,$circBank,$database,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"blast|b:s"=>\$blast_out,
				"circ2disease|cd:s"=>\$circ2disease,
				"circBank|cb:s"=>\$circBank,
				"dkey|k:s"=>\$database,
				"od|o:s"=>\$od
				) or &USAGE;
&USAGE unless ($blast_out and ($circ2disease || $circBank) and $od);

$blast_out = abs_path($blast_out);
$od = abs_path($od);

my(%disease,%bank,$head,$null);
if(defined $circ2disease){
	$circ2disease = abs_path($circ2disease);
	$database="Circ2disease";
#CircRNA_Name	Synonyms	Chr	Start-End	Strand	Genome_Version	Gene_Symbol	Sequence	Type	Disease_Name	Methods	Sample	Location	Expression_Pattern	miRNA_Validated	miRNA_Target_Validated	RBP_Validated	RBP_Target_Validated	Functional_Description	PubMed_ID	Year	Journal	Title	Note
	
	open(IN,"$circ2disease") or die $!;
	while(<IN>){
		chomp;
		my @tmp = split /\t/,$_;
		if(/^#/){
			my @ti = ("circRNA_name");
			splice(@tmp,0,1,@ti);
			splice(@tmp,2,3);
			splice(@tmp,4,1);
			$head=join("\t",@tmp);
			next;
		}
		splice(@tmp,2,3);
		splice(@tmp,4,1);
		$disease{$tmp[0]}=join("\t",@tmp);
		my @alias = split /; /,$tmp[1];
		foreach my $a (@alias){
			$disease{$a}=join("\t",@tmp);
		}
	}
	close(IN);
	
	$null = "\t--" x 20 if($database =~ /Circ2disease/i);
	my $out = "$od/circRNA_$database\_annotation.xls";
	open(IN,"$blast_out") or die $!;
	open(OUT,">$out") or die $!;
	while(<IN>){
		chomp;
		my @tmp = split /\t/,$_;
		if(/^#/){
			splice(@tmp,4);
			print OUT join("\t",@tmp),"\t$head\n";next;
		}
		my $id = (split /\|/,$tmp[1])[0];
		splice(@tmp,4);
		my $circ = join("\t",@tmp);
		if(exists $disease{$id}){
			print OUT "$circ\t$disease{$id}\n";
		}else{
			print OUT "$circ"."$null\n";
		}
	}
	close(IN);
	close(OUT);
}

#id	circRNA_id	chr	start	end	strand	spliced_seq_length	annotation	best_transcript	gene_symbol	mm9_circRNA_id	mm9_circbase_info	mm9_circbase_squence	mRNA_size	orf_size	fickett_score	hexamer_score	coding_prob	m6A_name
if(defined $circBank){
	$circBank = abs_path($circBank);
	$database = "circBank";
	open(BANK,"$circBank") or die $!;
	while(<BANK>){
		chomp;
		my @tmp = split /\t/,$_;
		if(/^#/){
			my @ti = ("circBank_id");
			splice(@tmp,0,1,@ti);
			splice(@tmp,1,6);
			splice(@tmp,7,1);
			$head = join("\t",@tmp);
			
			next;
		}
		my $id = $tmp[1];
		splice(@tmp,1,6);
		splice(@tmp,7,1);
		my $info = join("\t",@tmp);
		$bank{$id}=$info;
	}
	close(BANK);
	
	$null = "\t--" x 13 if($database =~ /circBank/i);
	my $out = "$od/circRNA_$database\_annotation.xls";
	open(IN,"$blast_out") or die $!;
	open(OUT,">$out") or die $!;
	while(<IN>){
		chomp;
		my @tmp = split /\t/,$_;
		if(/^#/){
			splice(@tmp,4);
			print OUT join("\t",@tmp),"\t$head\n";next;
		}
		my $id = (split /\|/,$tmp[1])[0];
		splice(@tmp,4);
		my $circ = join("\t",@tmp);
		if(exists $bank{$id}){
			print OUT "$circ\t$bank{$id}\n";
		}else{
			print OUT "$circ"."$null\n";
		}
	}
	close(IN);
	close(OUT);
}

########################################################
sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Contact:	Nie<niepy\@biomarker.com.cn>
Description:	This program is used to abstract informations from database(Circ2disease and circBank)...

Usage:
  Options:
	-h		 Help
	-blast|b	<file>	known circRNA blast out
	-circ2disease|cd	<file>	Circ2Disease/Circ2Disease_Association.txt,circBank
	-circBank|cb	<file>	circBank/circBank_annotation.txt
	-dkey|k		<str>	database name
	-od|o		<file>	outdir

USAGE
	print $usage;
	exit;
}
