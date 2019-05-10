#!/usr/bin/perl
use strict;

($#ARGV==3) or die "Usage: $0 gtf_file as_file -p prefix\n";
my $gtfFile = shift;
my $asFile = shift;
my $prefix = shift; 
($prefix eq "-p") or die "Usage: $0 gtf_file as_file -p prefix\n";
$prefix = shift;

my ($event_id_on,$event_type_on,$gene_id_on,$chr_on,$start_on,$end_on,$pattern_on,$strand_on);
my ($event_id_off,$event_type_off,$gene_id_off,$chr_off,$start_off,$end_off,$pattern_off,$strand_off);

my %SeenEvent;
my %SeenGene;

open(N, ">$prefix.as.nr") or die "cannot open $prefix.summary for writing.\n";
open(S, ">$prefix.as.summary") or die "cannot open $prefix.summary for writing.\n";

my $last_gene_id;
my @TSS;
my @TTS;
my @SKIP_ON;	my @XSKIP_ON;
my @SKIP_OFF;	my @XSKIP_OFF;
my @MSKIP_ON;	my @XMSKIP_ON;
my @MSKIP_OFF;	my @XMSKIP_OFF;
my @IR_ON;	my @XIR_ON;
my @IR_OFF;	my @XIR_OFF;
my @MIR_ON;	my @XMIR_ON;
my @MIR_OFF;	my @XMIR_OFF;
my @AE;		my @XAE;

print N "event_id\tevent_type\tgene_id\tchrom\tevent_start\tevent_end\tevent_pattern\tstrand\n";
print S "gene_id\tchrom\tnumTSS\tTSS_list\tnumTTS\tTTS_list\t";
print S "numSKIP_ON\tSKIP_ON_list\tnumSKIP_OFF\tSKIP_OFF_list\t";
print S "numMSKIP_ON\tMSKIP_ON_list\tnumMSKIP_OFF\tMSKIP_OFF_list\t";
print S "numIR_ON\tIR_ON_list\tnumIR_OFF\tIR_OFF_list\t";
print S "numMIR_ON\tMIR_ON_list\tnumMIR_OFF\tMIR_OFF_list\tnumAE\tAE_list\tevents_total\n";

open(F, "<$asFile") or die "cannot open $asFile.\n";
$_ = <F>;   # skip over header line
while (<F>) {
  chomp;

  /^\d+\t(\S+)\t(\S+)\t/ or die "illegal line. $_";
# next if ($1=~/^X/); # skip over X-versions of the events

  my $this_gene_id = $2;
  if (!defined($last_gene_id) || !($last_gene_id eq $this_gene_id)) {

    if (defined($last_gene_id)) {
      my ($t1,$t2,$t3,$t4,$t5,$t6,$t7,$t8,$t9,$t10,$t11);
      print S "$gene_id_on\t$chr_on\t",
              ($t1=scalar(@TSS)), "\t", join(',',@TSS), "\t",
              ($t2=scalar(@TTS)), "\t", join(',',@TTS), "\t",
              ($t3=scalar(@SKIP_ON)), "\t", join(',',@SKIP_ON), "\t",
              ($t4=scalar(@SKIP_OFF)), "\t", join(',',@SKIP_OFF), "\t",
              ($t5=scalar(@MSKIP_ON)), "\t", join(',',@MSKIP_ON), "\t",
              ($t6=scalar(@MSKIP_OFF)), "\t", join(',',@MSKIP_OFF), "\t",
              ($t7=scalar(@IR_ON)), "\t", join(',',@IR_ON), "\t",
              ($t8=scalar(@IR_OFF)), "\t", join(',',@IR_OFF), "\t",
              ($t9=scalar(@MIR_ON)), "\t", join(',',@MIR_ON), "\t",
              ($t10=scalar(@MIR_OFF)), "\t", join(',',@MIR_OFF), "\t",
              ($t11=scalar(@AE)), "\t", join(',',@AE), "\t",
              (1+ $t1-1+$t2-1+($t3 ? ($t3+$t4):0)+($t5 ? ($t5+$t6):0)+($t7 ? ($t7+$t8):0)+($t9 ? ($t9+$t10):0)+$t11), "\n";
    }
    unlink @TSS; unlink @TTS; unlink @SKIP_ON; unlink @SKIP_OFF;
    unlink @MSKIP_ON; unlink @MSKIP_OFF; unlink @IR_ON; unlink @IR_OFF;
    unlink @MIR_ON; unlink @MIR_OFF; unlink @AE;
    unlink @XSKIP_ON; unlink @XSKIP_OFF;
    unlink @XMSKIP_ON; unlink @XMSKIP_OFF; unlink @XIR_ON; unlink @XIR_OFF;
    unlink @XMIR_ON; unlink @XMIR_OFF; unlink @XAE;

    @TSS = (); @TTS = (); @SKIP_ON = (); @SKIP_OFF = ();
    @MSKIP_ON = (); @MSKIP_OFF = (); @IR_ON = (); @IR_OFF = ();
    @MIR_ON = (); @MIR_OFF = (); @AE = ();
    @XSKIP_ON = (); @XSKIP_OFF = ();
    @XMSKIP_ON = (); @XMSKIP_OFF = (); @XIR_ON = (); @XIR_OFF = ();
    @XMIR_ON = (); @XMIR_OFF = (); @XAE = (); 

    $SeenGene{$this_gene_id} = 1;
  }

  #1000008   SKIP_ON aB.10831        chr1    15823482        15823570        15823482-15823570       +       
  if (/_ON/) {
    /^(\d+)\t(\S+)\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\S)\t/ or die "illegal line. $_"; 
    ($event_id_on,$event_type_on,$gene_id_on,$chr_on,$start_on,$end_on,$pattern_on,$strand_on) = ($1,$2,$3,$4,$5,$6,$7,$8);
    $_ = <F>; chomp;
    /^(\d+)\t(\S+_OFF)\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\S)\t/ or die "illegal line. $_";
    ($event_id_off,$event_type_off,$gene_id_off,$chr_off,$start_off,$end_off,$pattern_off,$strand_off) = ($1,$2,$3,$4,$5,$6,$7,$8);

    next if (defined($SeenEvent{"$event_id_on:$event_id_off"})); 

    print N "$event_id_on\t$event_type_on\t$gene_id_on\t$chr_on\t$start_on\t$end_on\t$pattern_on\t$strand_on\n";
    print N "$event_id_off\t$event_type_off\t$gene_id_off\t$chr_off\t$start_off\t$end_off\t$pattern_off\t$strand_off\n";

    $SeenEvent{"$event_id_on:$event_id_off"} = 1;

    if ($event_type_on eq "SKIP_ON") {
      push (@SKIP_ON,$event_id_on);
      push (@SKIP_OFF,$event_id_off);
    } elsif ($event_type_on eq "MSKIP_ON") {
      push (@MSKIP_ON,$event_id_on);
      push (@MSKIP_OFF,$event_id_off);
    } elsif ($event_type_on eq "IR_ON") {
      push (@IR_ON,$event_id_on);
      push (@IR_OFF,$event_id_off);
    } elsif ($event_type_on eq "MIR_ON") {
      push (@MIR_ON,$event_id_on);
      push (@MIR_OFF,$event_id_off);
    } elsif ($event_type_on eq "XSKIP_ON") {
      push (@XSKIP_ON,$event_id_on);
      push (@XSKIP_OFF,$event_id_off);
    } elsif ($event_type_on eq "XMSKIP_ON") {
      push (@XMSKIP_ON,$event_id_on);
      push (@XMSKIP_OFF,$event_id_off);
    } elsif ($event_type_on eq "XIR_ON") {
      push (@XIR_ON,$event_id_on);
      push (@XIR_OFF,$event_id_off);
    } elsif ($event_type_on eq "XMIR_ON") {
      push (@XMIR_ON,$event_id_on);
      push (@XMIR_OFF,$event_id_off);
    }

  } else {
    /^(\d+)\t(\S+)\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\S)\t/ or die "illegal line. $_";
    ($event_id_on,$event_type_on,$gene_id_on,$chr_on,$start_on,$end_on,$pattern_on,$strand_on) = ($1,$2,$3,$4,$5,$6,$7,$8);

    next if (defined($SeenEvent{$event_id_on}));

    print N "$event_id_on\t$event_type_on\t$gene_id_on\t$chr_on\t$start_on\t$end_on\t$pattern_on\t$strand_on\n";

    if ($event_type_on eq "TSS") {
      push (@TSS,$event_id_on);
    } elsif ($event_type_on eq "TTS") {
      push (@TTS,$event_id_on);
    } elsif ($event_type_on eq "AE") {
      push (@AE,$event_id_on);
    } elsif ($event_type_on eq "XAE") {
      push (@XAE,$event_id_on);
    }

    $SeenEvent{$event_id_on} = 1;
  }

  $last_gene_id = $gene_id_on;
}
if (defined($last_gene_id)) {
  my ($t1,$t2,$t3,$t4,$t5,$t6,$t7,$t8,$t9,$t10,$t11);
  print S "$gene_id_on\t$chr_on\t",
          ($t1=scalar(@TSS)), "\t", join(',',@TSS), "\t",
          ($t2=scalar(@TTS)), "\t", join(',',@TTS), "\t",
          ($t3=scalar(@SKIP_ON)), "\t", join(',',@SKIP_ON), "\t",
          ($t4=scalar(@SKIP_OFF)), "\t", join(',',@SKIP_OFF), "\t",
          ($t5=scalar(@MSKIP_ON)), "\t", join(',',@MSKIP_ON), "\t",
          ($t6=scalar(@MSKIP_OFF)), "\t", join(',',@MSKIP_OFF), "\t",
          ($t7=scalar(@IR_ON)), "\t", join(',',@IR_ON), "\t",
          ($t8=scalar(@IR_OFF)), "\t", join(',',@IR_OFF), "\t",
          ($t9=scalar(@MIR_ON)), "\t", join(',',@MIR_ON), "\t",
          ($t10=scalar(@MIR_OFF)), "\t", join(',',@MIR_OFF), "\t",
          ($t11=scalar(@AE)), "\t", join(',',@AE), "\t",
          (1+ $t1-1+$t2-1+($t3 ? ($t3+$t4):0)+($t5 ? ($t5+$t6):0)+($t7 ? ($t7+$t8):0)+($t9 ? ($t9+$t10):0)+$t11), "\n";

}

# Now add genes that are expressed, but don't have splice variation (i.e., single exon; others?)
open(G, "<$gtfFile") or die "cannot open GTF file $gtfFile.\n";
while (<G>) {
  chomp;

  next if !/\ttranscript\t/;

  m/gene_id \"(\S+)\"; transcript_id/ or die "illegal GTF transcript line. $_";
  my $gene_id = $1;
  /^(\S+)\t/; 
  my $chr = $1; 

  next if defined($SeenGene{$gene_id});

  print S "$gene_id\t$chr\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t1\n";

  $SeenGene{$gene_id} = 1; 
}

close(F);
close(N);
close(S);
close(G);
