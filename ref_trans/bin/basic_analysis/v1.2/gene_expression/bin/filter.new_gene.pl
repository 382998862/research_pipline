#!usr/bin/perl -w
use strict;
use FindBin qw($Bin $Script);
use Cwd 'abs_path';

die "Usage:<species> <genome> <newgene dir><one exon>\n" unless (@ARGV==3 or @ARGV==4) ;
my ($species,$genome,$new_dir,$one_exon);
	($species,$genome,$new_dir)=@ARGV if @ARGV==3;
	($species,$genome,$new_dir,$one_exon)=@ARGV if @ARGV==4;
$genome  = abs_path($genome);
$new_dir = abs_path($new_dir);

my (%filtered_1, %filtered);
&filter_orf;
&filter_exon;

open IN, "$new_dir/$species.newGene_final_tracking.list";
open OUT, ">", "$new_dir/$species.newGene_filtered_final_tracking.list";
while(<IN>)
{
	chomp;
	my @data = split /\t+/;
	next unless exists $filtered{$data[7]};
	print OUT "$_\n";
}
close IN;
close OUT;

#system "rm -r $new_dir/tmp_orf";

sub filter_orf
{
	chdir($new_dir);
	system "perl $Bin/gff_gene2transcript_v1.0.pl $genome $new_dir/$species.newGene_final.gff $species.newGene";
	my $fa_file = "$new_dir/$species.newGene.transcript.fa";
	mkdir("tmp_orf") unless -d 'tmp_orf';
	system "perl $Bin/Getorf_Extract.pl -fa $fa_file -od $new_dir/tmp_orf";
	open (IN,"$new_dir/tmp_orf/$species.newGene.transcript.orf") || die $!;
	my %len;
	while(<IN>)
	{
		chomp;
		next unless /^>/;
		my $trans_id = (split /\s+/)[0];
		$trans_id =~ />(\S+)_(\d+)$/;
		$trans_id = $1;
		my $seq = <IN>;             chomp $seq;
		my $len = length($seq);
		$len{$trans_id} = $len unless exists $len{$trans_id};
		$len{$trans_id} = ($len > $len{$trans_id}) ? $len : $len{$trans_id};
	}
	close IN;
	foreach my $trans_id(keys %len)
	{
		$filtered_1{$trans_id}++ if $len{$trans_id} >= 50;
	}
}
sub filter_exon
{
	my $gff_file = "$new_dir/$species.newGene_final.gff";
	open (IN, $gff_file) || die $!;
	my (%info, %exon_num, %exon_parent);
	while(<IN>)
	{
		chomp;
		my @data = split /\t+/;
		if($data[2] eq 'gene')
		{
			$data[-1] =~ /ID=(\S+)$/;
			my $gene_id = $1;
			$info{$gene_id}{'gene'}{$gene_id} = "$_\n";
		}elsif($data[2] eq 'mRNA')
		{
			$data[-1] =~ /ID=(\S+);Parent=(\S+)$/;
			my ($trans_id, $gene_id) = ($1, $2);
			$exon_parent{$trans_id} = $gene_id;
			$info{$gene_id}{'mRNA'}{$trans_id} = "$_\n";
		}elsif($data[2] eq 'CDS')
		{
			$data[-1] =~ /Parent=(\S+)$/;
			my $trans_id = $1;
			my $gene_id = $exon_parent{$trans_id};
			$info{$gene_id}{'mRNA'}{$trans_id} .= "$_\n";
			$exon_num{$trans_id}++;
		}
	}
	close IN;
	my $filter_gff = "$new_dir/$species.newGene_final.filtered.gff";
	open OUT, ">", $filter_gff;
	foreach my $gene_id(sort keys %info)
	{
		my $lines = $info{$gene_id}{'gene'}{$gene_id};
		my $flag;
		foreach my $trans_id(sort keys %{$info{$gene_id}{'mRNA'}})
		{
			next unless exists $filtered_1{$trans_id};
			my $exon_num = $exon_num{$trans_id};
			next if ((!defined $one_exon) and ($exon_num == 1)) ;
			$filtered{$trans_id}++;
			$lines .= $info{$gene_id}{'mRNA'}{$trans_id};
			$flag++;
		}
		print OUT $lines if $flag;
	}
	close OUT;
}

