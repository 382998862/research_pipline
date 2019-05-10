#!/usr/bin/perl -w
# 
# Copyright (c) BIO_MARK 2014
# Writer:         Dengdj <dengdj@biomarker.com.cn>
# Program Date:   2014
# Modifier:       Dengdj <dengdj@biomarker.com.cn>
# Last Modified:  2014.
my $ver="1.1";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

######################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作

my %opts;
GetOptions(\%opts,"i=s","o=s","m=s","x=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{i}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver

	Usage:

		-i           vcf file of indel(gatk anno)         <infile>        must be given

		-o           statistic result file prefix         <outfile>       must be given

		-m           minimum indel length for draw, default -10           optional

		-x           maximum indel length for draw, default 10            optional

		-h           Help document

	Usage End.

	exit;
}

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
################
# get parameters
my $infile  = $opts{i} ;
my $outfile = $opts{o} ;
my $min_length = defined $opts{m} ? $opts{m} : -10 ;
my $max_length = defined $opts{x} ? $opts{x} : 10 ;

# reading indel vcf file
my %hindel = ();
my %hstat = ();
my $asamples = &reading_indel_vcf_file($infile, \%hindel, \%hstat);

# output stat result 
&output_stat_result($outfile, $asamples, \%hindel, \%hstat);

# draw indel distribution
&draw_indel_distribution("$outfile.dis");

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

###############Subs
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;
	$cur_dir =~ s/\n$//;
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;
		$dir =~ s/\n$// ;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;
		$return =~ s/\n$// ;
	}
	else
	{
		warn "Warning just for file and dir\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}

# &show_log("txt")
sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}

#&run_or_die($cmd);
sub run_or_die()
{
	my ($cmd) = @_ ;
	&show_log($cmd);
	my $flag = system($cmd) ;
	if ($flag != 0){
		&show_log("Error: command fail: $cmd");
		exit(1);
	}
	&show_log("done.");
	return ;
}


#$asamples = &reading_indel_vcf_file($infile, \%hindel, \%hstat);
sub reading_indel_vcf_file()
{
	my ($infile, $ahindel, $ahstat) = @_ ;
	open (IN, $infile) || die "Can't open $infile, $!\n" ;
	my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format, @samples) = ();
	while (<IN>){
		chomp ;
		next if (m/^\s*$/) ;
		if (m/^\#CHROM/){
			($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@samples)=(split/\s+/,$_);
		}
		next if (m/^\#/);
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format, @types) = split ;
		&stat_each_sample($ref, $alt, $info, \@types, \@samples, $ahindel, $ahstat);
	}
	close(IN);

	return(\@samples) ;
}

#&stat_each_sample($ref, $alt, $info, \@types, \@samples, $ahindel, $ahstat);
sub stat_each_sample()
{
	my ($ref, $alt, $info, $atypes, $asamples, $ahindel, $ahstat) = @_ ;
	my $genome_type = &judge_genome_type($info);
	my ($indel_type, $len) = &judge_indel_type($ref, $alt);
	for (my $i=0; $i<@{$atypes}; $i++){
		my $heterhomo = &judge_heter_type($atypes->[$i]);
		if ($heterhomo ne 'NA'){
			$ahindel->{$asamples->[$i]}->{"genomic"}->{$heterhomo}++ ;
			$ahindel->{$asamples->[$i]}->{"genomic"}->{$indel_type}++ ;
			$ahstat->{$asamples->[$i]}->{"genomic"}->{$len}++ ;
			if ($genome_type eq "cds"){
				$ahindel->{$asamples->[$i]}->{"cds"}->{$heterhomo}++ ;
				$ahindel->{$asamples->[$i]}->{"cds"}->{$indel_type}++ ;
				$ahstat->{$asamples->[$i]}->{"cds"}->{$len}++ ;
			}
		}
	}
	$ahindel->{"all"}->{"genomic"}->{$indel_type}++ ;
	$ahindel->{"all"}->{"cds"}->{$indel_type}++ if ($genome_type eq "cds") ;
	$ahstat->{"all"}->{"genomic"}->{$len}++ ;
	$ahstat->{"all"}->{"cds"}->{$len}++ if ($genome_type eq "cds");

	return ;
}

#my $heter = &judge_heter_type($atypes->[$i]);
sub judge_heter_type()
{
	my ($type) = @_ ;
	my $hethomo_type ;
	my $genotype = (split /\:/, $type)[0] ;
	my ($alle1, $alle2) = split /\//, $genotype ;
	if ($alle1 ne $alle2){
		$hethomo_type = "Heter" ;
	}
	elsif ($alle1 ne '0' && $alle1 ne '.'){
		$hethomo_type = "Homo" ;
	}
	else{
		$hethomo_type = "NA" ;
	}

	return($hethomo_type) ;
}

#my $genome_type = &judge_genome_type($info);
sub judge_genome_type(){
	my ($info) = @_ ;
	my $type = "genomic" ;
	my %hinfo = (
		'NONE'=>'genomic', 'CHROMOSOME'=>'genomic', 'CUSTOM'=>'genomic', 'CDS'=>'cds', 'INTERGENIC'=>'genomic',
		'INTERGENIC_CONSERVED'=>'genomic', 'UPSTREAM'=>'genomic', 'UTR_5_PRIME'=>'genomic', 'UTR_5_DELETED'=>'genomic',
		'START_GAINED'=>'genomic', 'SPLICE_SITE_ACCEPTOR'=>'genomic', 'SPLICE_SITE_DONOR'=>'genomic', 'SPLICE_SITE_REGION'=>'genomic',
		'INTRAGENIC'=>'genomic', 'START_LOST'=>'cds', 'SYNONYMOUS_START'=>'cds', 'NON_SYNONYMOUS_START'=>'cds',
		'GENE'=>'genomic', 'TRANSCRIPT'=>'genomic', 'EXON'=>'cds', 'EXON_DELETED'=>'cds', 'NON_SYNONYMOUS_CODING'=>'cds',
		'SYNONYMOUS_CODING'=>'cds', 'FRAME_SHIFT'=>'cds', 'CODON_CHANGE'=>'cds', 'CODON_INSERTION'=>'cds',
		'CODON_CHANGE_PLUS_CODON_INSERTION'=>'cds', 'CODON_DELETION'=>'cds', 'CODON_CHANGE_PLUS_CODON_DELETION'=>'cds',
		'STOP_GAINED'=>'cds', 'SYNONYMOUS_STOP'=>'cds', 'STOP_LOST'=>'cds', 'RARE_AMINO_ACID'=>'cds', 'INTRON'=>'genomic',
		'INTRON_CONSERVED'=>'genomic', 'UTR_3_PRIME'=>'genomic', 'UTR_3_DELETED'=>'genomic', 'DOWNSTREAM'=>'genomic',
		'REGULATION'=>'genomic'
	);
	if ($info =~ /SNPEFF_EFFECT=(.*?)\;/){
		$type = $hinfo{$1} ;
	}

	return($type);
}

#my ($indel_type, $len) = &judge_indel_type($ref, $alt);
sub judge_indel_type()
{
	my ($ref, $alt) = @_ ;
	my $type ;
	$alt =~ s/\,.*$// ;
	my $len = length($alt) - length($ref) ;
	if ($len < 0){
		$type = "DELETION" ;
	}
	elsif ($len > 0){
		$type = "INSERT" ;
	}
	else{
		print "$ref, $alt is not an indel\n" ;
		exit(1);
	}

	return($type, $len) ;
}

#&output_stat_result($outfile, $asamples, \%hindel, \%hstat);
sub output_stat_result()
{
	my ($outfile, $asamples, $ahindel, $ahstat) = @_ ;
	# stat result
	open (OUT, ">$outfile.stat") || die "Can't creat $outfile, $!\n" ;
	print OUT "#Sample\tCDS-Insertion\tCDS-Deletion\tCDS-Het\tCDS-Homo\tCDS-Total\tGenome-Insertion\tGenome-Deletion\tGenome-Het\tGenome-Homo\tGenome-Total\n" ;
	for (my $i=0; $i<@{$asamples}; $i++){
		my $cds_ins = defined $ahindel->{$asamples->[$i]}->{'cds'}->{'INSERT'} ? $ahindel->{$asamples->[$i]}->{'cds'}->{'INSERT'} : 0 ;
		my $cds_del = defined $ahindel->{$asamples->[$i]}->{'cds'}->{'DELETION'} ? $ahindel->{$asamples->[$i]}->{'cds'}->{'DELETION'} : 0 ;
		my $cds_homo = defined $ahindel->{$asamples->[$i]}->{'cds'}->{'Homo'} ? $ahindel->{$asamples->[$i]}->{'cds'}->{'Homo'} : 0 ;
		my $cds_heter = defined $ahindel->{$asamples->[$i]}->{'cds'}->{'Heter'} ? $ahindel->{$asamples->[$i]}->{'cds'}->{'Heter'} : 0 ;
		my $cds_total = $cds_homo + $cds_heter ;
		my $genome_ins = defined $ahindel->{$asamples->[$i]}->{'genomic'}->{'INSERT'} ? $ahindel->{$asamples->[$i]}->{'genomic'}->{'INSERT'} : 0 ;
		my $genome_del = defined $ahindel->{$asamples->[$i]}->{'genomic'}->{'DELETION'} ? $ahindel->{$asamples->[$i]}->{'genomic'}->{'DELETION'} : 0 ;
		my $genome_homo = defined $ahindel->{$asamples->[$i]}->{'genomic'}->{'Homo'} ? $ahindel->{$asamples->[$i]}->{'genomic'}->{'Homo'} : 0 ;
		my $genome_heter = defined $ahindel->{$asamples->[$i]}->{'genomic'}->{'Heter'} ? $ahindel->{$asamples->[$i]}->{'genomic'}->{'Heter'} : 0 ;
		my $genome_total = $genome_homo + $genome_heter ;
		print OUT $asamples->[$i],"\t", $cds_ins,"\t", $cds_del,"\t", $cds_homo,"\t", $cds_heter,"\t", $cds_total ;
		print OUT "\t", $genome_ins,"\t", $genome_del,"\t", $genome_homo,"\t", $genome_heter,"\t", $genome_total,"\n" ;
	}
	my $total_cds_ins = defined $ahindel->{'all'}->{'cds'}->{'INSERT'} ? $ahindel->{'all'}{'cds'}{'INSERT'} : 0 ;
	my $total_cds_del = defined $ahindel->{'all'}{'cds'}{'DELETION'} ? $ahindel->{'all'}{'cds'}{'DELETION'} : 0 ;
	my $total_cds = $total_cds_ins + $total_cds_del ;
	my $total_genome_ins = defined $ahindel->{'all'}->{'genomic'}->{'INSERT'} ? $ahindel->{'all'}{'genomic'}{'INSERT'} : 0 ;
	my $total_genome_del = defined $ahindel->{'all'}{'genomic'}{'DELETION'} ? $ahindel->{'all'}{'genomic'}{'DELETION'} : 0 ;
	my $total_genome = $total_genome_ins + $total_genome_del ;
	print OUT join("\t", ("#Total", $total_cds_ins, $total_cds_del, "--", "--", $total_cds, $total_genome_ins, $total_genome_del, "--", "--", $total_genome)), "\n" ;
	close(OUT);
	# length distribution
	open (OUT, ">$outfile.dis") || die "Can't creat $outfile.dis, $!\n" ;
	print OUT "#Length\tCDS/Genome\tNumber\n" ;
	for my $len (sort {$a<=>$b} keys %{$ahstat->{'all'}->{'genomic'}}){
		my $genome_num = defined $ahstat->{'all'}{'genomic'}{$len} ? $ahstat->{'all'}{'genomic'}{$len} : 0 ;
		my $cds_num = defined $ahstat->{'all'}{'cds'}{$len} ? $ahstat->{'all'}{'cds'}{$len} : 0 ;
		print OUT "$len\tGenome\t$genome_num\n" ;
		print OUT "$len\tCDS\t$cds_num\n" ;
	}
	close(OUT);

	return ;
}

#&draw_indel_distribution($disfile);
sub draw_indel_distribution()
{
	my ($disfile) = @_ ;
	(my $out_png_file = $disfile) =~ s/$/.png/ ;
	my $cmd = "perl $Bin/indel_distribution.pl -i $disfile -o $out_png_file -m $min_length -x $max_length" ;
	&run_or_die($cmd);

	return ;
}

