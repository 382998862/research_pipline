#!usr/bin/perl -w
use strict;

use Getopt::Long;

my $usage="
##*************************************************************************************##
Description: This scirpt is used to check the gff file and genome fa file whether matched.

Options:
      -i_gff     input gff file
      -i_fa      input genome file of .fa format
      -o_gff     output gff file after inspecting
      -o_fa      output genome file of .fa format after inspecting

Users can only check gff or fa with only one type of file parameters supplied.
##*************************************************************************************##
\n";


my ($gff,     $genome);
my ($gff_out, $genome_out);


GetOptions
(
	"i_gff:s" =>\$gff,               "i_fa:s"  =>\$genome,
    "o_gff:s" =>\$gff_out,           "o_fa:s"  =>\$genome_out,
)|| die $usage;

die $usage unless $gff && $gff_out && $genome && $genome_out;     ## some different mode can be set to do different process, this function could be added

my %gff_info;
read_gff();
my %chr_info;
read_genome();

foreach my $chr (keys %gff_info)
{
	print "$chr\tThere is chromosome from gff dosen't appear in genome!\n" unless exists $chr_info{$chr};
}


open (OUT, ">", $gff_out) || die $!;
foreach my $chr (keys %gff_info)
{
	foreach my $gene_id (sort {$gff_info{$chr}{'gene'}{$a}->[0]<=>$gff_info{$chr}{'gene'}{$b}->[0]} keys %{$gff_info{$chr}{'gene'}})
	{
		my $block;

		my $flag;
	
		my $gene_line  = $gff_info{$chr}{'gene'}{$gene_id}->[1];       
		my @trans_ids  = keys %{$gff_info{$chr}{$gene_id}};

		unless (@trans_ids)
		{
			print "gene_id: $gene_id have no trans, please check!\n";
			next;
		}
		$block .= $gene_line."\n";

		foreach (@trans_ids)
		{
			my $trans_line = $gff_info{$chr}{'mRNA'}{$_};         
			my @cds_lines = sort {$gff_info{$chr}{$_}{'CDS'}{$a}<=>$gff_info{$chr}{$_}{'CDS'}{$b}} keys %{$gff_info{$chr}{$_}{'CDS'}};
			unless (@cds_lines)
			{
				print "trans_id: $_ have no cds, please check!\n";
				next;
			}
			$block .= $trans_line."\n";
			my $cds_lines = join "\n", @cds_lines;         
			$block .= $cds_lines."\n";
	
			$flag++;
		}
		unless ($flag)
		{
			print "gene_id: $gene_id gene-mRNA-cds information not intact, please check!\n";
			next;
		}
		print OUT $block;
	}
}
close OUT;


sub read_genome
{
	open (IN, $genome) || die $!;            open OUT, ">", $genome_out;
	while (<IN>)
	{
		chomp;
		if(/>(\S+)/)  
		{
			$chr_info{$1}++;
			print OUT ">$1\n";
			next;   
		}
		next if /^\s*$/;
		print OUT "$_\n"; 
	}
	close IN;                                   close OUT;
}

sub read_gff
{
	open IN, $gff;
	while(<IN>)
	{
		chomp;        next if (/^\s*$/ || /^#/);
		my @data = split /\t+/;               next unless scalar(@data) == 9;
		if($data[2] eq 'gene')
		{
			my $id = gene_check($data[8]);
			my $gene_line = join "\t", (@data[0..7], "ID=$id");
			$gff_info{$data[0]}{'gene'}{$id} = [$data[3],$gene_line];
		}elsif($data[2] eq 'mRNA')
		{
			my ($id, $parent) = mRNA_check($data[8]);
			my $mRNA_line = join "\t", (@data[0..7], "ID=$id;Parent=$parent");
			$gff_info{$data[0]}{'mRNA'}{$id} = $mRNA_line;
			$gff_info{$data[0]}{$parent}{$id}++;
		}elsif($data[2] eq 'CDS')
		{
			my $parent = CDS_check($data[8]);
			my $cds_line = join "\t", (@data[0..7], "Parent=$parent");
			$gff_info{$data[0]}{$parent}{'CDS'}{$cds_line}= $data[3];
		}
	}
	close IN;
}

sub strand_check
{
	my ($strand_para, $line) = @_;
	unless ($strand_para eq '+' || $strand_para eq '-')
	{
		die "There is unexpected character in strand column\n$line\n";
	}
}

sub position_check
{
	my ($position, $line) = @_;
	unless ($position =~ /^\d+$/)
	{
		die "Position must be integer!\n$line\n";
	}
}

sub gene_check
{
	my ($info, $line) = @_;
	my ($id) = $info=~ /ID=([^;]+)/;
	unless ($id)
	{
		die "Check your file with 'gene' row!\n$line\nIt should contain 'ID=**'\n";
	}
	
	return $id;
}

sub mRNA_check
{
	my ($info, $line) = @_;
	my ($id) = $info =~ /ID=([^;]+)/;
	my ($parent) = $info =~ /Parent=([^;]+)/;
	unless($id && $parent)
	{
		die "Check your file with 'mRNA' row!\n$line\nIt should contain 'ID=**;Parent=**'\n";
	}

	return ($id, $parent);
}

sub CDS_check
{
	my ($info, $line) = @_;
	my ($parent) = $info =~ /Parent=([^;]+)/;
	unless($parent)
	{
		die "Check your file with 'CDS' row!\n$line\nIt should contain 'Parent=**'\n";
	}

	return $parent;
}
