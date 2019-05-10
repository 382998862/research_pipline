#!/usr/bin/perl
use strict;

($#ARGV==2) or die "Usage: $0 set1,set2,... -s set_id\n";

my $tmp = shift;
my @SETS  = split ',', $tmp;
my $SUFFIX = shift;
($SUFFIX eq "-s") or die "Usage: $0 set1,set2,... -s set_id\n";
$SUFFIX = shift;
my $MODEL = "$SETS[0].$SUFFIX";

my %FPKM;

foreach my $set (@SETS) {
  open(F, "<$set.$SUFFIX") or die "could not open file $set.\n";
  $_ = <F>; # skip first line, but first check
  /^event_id/ or print STDERR "header line expected. $_"; 

  while (<F>) {
    chomp;

    /^(\d+)\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t(\S+)$/ or die "illegal line in set $set.$_";
    $FPKM{"$set:$1"} = $2;
  } close(F);
  print STDERR "done $set.\n";
}
 
my @gene = ();
my $gene_id;
my %TSS_gene;
my %TTS_gene;
my %SKIP_event;
my %MSKIP_event;
my %IR_event;
my %MIR_event;
my %AE_gene;


open(F, "<$MODEL") or die "could not open file $MODEL.\n";

# print header line
print "event_id\tevent_type\tgene_id\tchrom\tevent_start\tevent_end\tevent_pattern\tstrand";
foreach my $set (@SETS) { print "\t$set.fpkm"; }
foreach my $set (@SETS) { print "\t$set.ratio"; }
for (my $i=0; $i<scalar(@SETS); $i++) {
  for (my $j=$i+1; $j<scalar(@SETS); $j++) {
    print "\t" . $SETS[$i] . "-" . $SETS[$j];
  }
}
print "\tnum_stages\n";


$_ = <F>; # skip first line, but check first
/^event_id/ or print STDERR "unexpected header line. $_";
while (<F>) {
  chomp;

  /^(\d+)\t\S+\t(\S+)\t/ or die "illegal line in model. $_";
  my $this_event_id = $1;
  my $this_gene_id = $2;

  if (defined($gene_id) && !($gene_id eq $this_gene_id)) {
    # process events from the previous gene

    for (my $i=0; $i<scalar(@gene); $i++) {
      my $l = $gene[$i];
      next if ($l=~/_OFF/);

      my %Ratio;

      $l =~ /^(\d+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)$/ or die "illegal line in model.$_";
      my ($event_id,$event_type,$gene_id,$chrom,$start,$end,$pattern,$ori,$fpkm) = ($1,$2,$3,$4,$5,$6,$7,$8,$9);

      print "$event_id\t$event_type\t$gene_id\t$chrom\t$start\t$end\t$pattern\t$ori";
      foreach my $set (@SETS) {
        print "\t$set:", $FPKM{"$set:$event_id"};
      }

      my $numStages = scalar(@SETS);
      ## calculate ratios, separately in each set, and print
      foreach my $set (@SETS) {
        if ($l=~/\tTSS/) {
          $Ratio{$set} = ($TSS_gene{$set}>=0.0000001) ? $FPKM{"$set:$event_id"}/$TSS_gene{$set} : -1;
        } elsif ($l=~/\tTTS/) {
          $Ratio{$set} = ($TTS_gene{$set}>=0.0000001) ? $FPKM{"$set:$event_id"}/$TTS_gene{$set} : -1;
        } elsif ($l=~/\tAE/) {
          $Ratio{$set} = ($AE_gene{$set}>=0.0000001) ? $FPKM{"$set:$event_id"}/$AE_gene{$set} : -1;
        } elsif ($l=~/\tSKIP_ON/) {
          $gene[$i+1] =~ /^(\d+)\t\S+_OFF\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+$/ or die "illegal line (OFF) in model.$_";
          my $event_id_2 = $1;

          $Ratio{$set} = ($SKIP_event{"$set:$event_id:$event_id_2"}>=0.0000001) ? 
                         ($FPKM{"$set:$event_id"}/(1.0*$SKIP_event{"$set:$event_id:$event_id_2"})) : -1;

        } elsif ($l=~/\tMSKIP_ON/) {
          $gene[$i+1] =~ /^(\d+)\t\S+_OFF\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+$/ or die "illegal line (OFF) in model.$_";
          my $event_id_2 = $1;

          $Ratio{$set} = ($MSKIP_event{"$set:$event_id:$event_id_2"}>=0.0000001) ?
                         $FPKM{"$set:$event_id"}/$MSKIP_event{"$set:$event_id:$event_id_2"} : -1;
        } elsif ($l=~/\tIR_ON/) {
          $gene[$i+1] =~ /^(\d+)\t\S+_OFF\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+$/ or die "illegal line (OFF) in model.$_";
          my $event_id_2 = $1;

          $Ratio{$set} = ($IR_event{"$set:$event_id:$event_id_2"}>=0.0000001) ?
                         $FPKM{"$set:$event_id"}/$IR_event{"$set:$event_id:$event_id_2"} : -1;
        } elsif ($l=~/\tMIR_ON/) {
          $gene[$i+1] =~ /^(\d+)\t\S+_OFF\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+$/ or die "illegal line (OFF) in model.$_";
          my $event_id_2 = $1;

          $Ratio{$set} = ($MIR_event{"$set:$event_id:$event_id_2"}>=0.0000001) ?
                         $FPKM{"$set:$event_id"}/$MIR_event{"$set:$event_id:$event_id_2"} : -1;
        }

        if ($Ratio{$set}>=0) { printf "\t$set:%1.6f", $Ratio{$set}; } else { print "\t$set:NA"; $numStages--; }
      }        

      # now print out differences between pairs of tissues; this will need to be augmented with
      # some statistical test for multiple samples
      for (my $i=0; $i<scalar(@SETS); $i++) {
        for (my $j=$i+1; $j<scalar(@SETS); $j++) {
          print "\t" . $SETS[$i] . "-" . $SETS[$j] . ":";
          if ($Ratio{$SETS[$i]}>=0 && $Ratio{$SETS[$j]}>=0) {
             print $Ratio{$SETS[$i]}-$Ratio{$SETS[$j]};
          } else {
             print "NA";
          }
        }
      }

      unlink %Ratio;

      print "\t$numStages\n";
    }  # end process previous gene

    # now reset counters and sums for the new gene 
    foreach my $set (@SETS) {
      $TSS_gene{$set} = 0; $TTS_gene{$set} = 0; $AE_gene{$set} = 0;
      $SKIP_event{$set} = 0; $MSKIP_event{$set} = 0; $IR_event{$set} = 0; $MIR_event{$set} = 0;
    }

    @gene = ();
  }

  push @gene, $_;

  if (/\tTSS/) {
    foreach my $set (@SETS) { $TSS_gene{$set} += $FPKM{"$set:$this_event_id"}; }
  } elsif (/\tTTS/) {
    foreach my $set (@SETS) { $TTS_gene{$set} += $FPKM{"$set:$this_event_id"}; }
  } elsif (/\tAE/) {
    foreach my $set (@SETS) { $AE_gene{$set} += $FPKM{"$set:$this_event_id"}; }
  } elsif (/\tSKIP_ON/) {
    $_ = <F>; chomp;
    /^(\d+)\t.*SKIP_OFF\t.*/ or die "illegal line to follow SKIP_ON. $_"; 
    my $this_event_id2 = $1;

    push @gene, $_;

    foreach my $set (@SETS) {
      $SKIP_event{"$set:$this_event_id:$this_event_id2"} = $FPKM{"$set:$this_event_id"} + $FPKM{"$set:$this_event_id2"};
    }
  } elsif (/\tMSKIP_ON/) {
    $_ = <F>; chomp;
    /^(\d+)\t.*MSKIP_OFF\t.*/ or die "illegal line to follow MSKIP_ON. $_";
    my $this_event_id2 = $1;

    push @gene, $_;

    foreach my $set (@SETS) {
      $MSKIP_event{"$set:$this_event_id:$this_event_id2"} = $FPKM{"$set:$this_event_id"} + $FPKM{"$set:$this_event_id2"};
    }
  } elsif (/\tIR_ON/) {
    $_ = <F>; chomp;
    /^(\d+)\t.*IR_OFF\t.*/ or die "illegal line to follow IR_ON. $_";
    my $this_event_id2 = $1;

    push @gene, $_;

    foreach my $set (@SETS) {
      $IR_event{"$set:$this_event_id:$this_event_id2"} = $FPKM{"$set:$this_event_id"} + $FPKM{"$set:$this_event_id2"};
    }
  } elsif (/\tMIR_ON/) {
    $_ = <F>; chomp;
    /^(\d+)\t.*MIR_OFF\t.*/ or die "illegal line to follow MIR_ON. $_";
    my $this_event_id2 = $1;

    push @gene, $_;

    foreach my $set (@SETS) {
      $MIR_event{"$set:$this_event_id:$this_event_id2"} = $FPKM{"$set:$this_event_id"} + $FPKM{"$set:$this_event_id2"};
    }
  } elsif (/_OFF/) {
    die "Error: Unexpected _OFF without _ON. $_";
  } else {
    die "Unrecognized event type code. $_";
  }
 
  $gene_id = $this_gene_id;
}
close(F);

# process last gene
if (defined($gene_id)) {
  # process events from the previous gene

  for (my $i; $i<scalar(@gene); $i++) {
    my $l = $gene[$i];
    next if ($l=~/_OFF/);

    my %Ratio;

    $l =~ /^(\d+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)$/ or die "illegal line in model.$_";
    my ($event_id,$event_type,$gene_id,$chrom,$start,$end,$pattern,$ori,$fpkm) = ($1,$2,$3,$4,$5,$6,$7,$8,$9);

    print "$event_id\t$event_type\t$gene_id\t$chrom\t$start\t$end\t$pattern\t$ori";
    foreach my $set (@SETS) {
      print "\t$set:", $FPKM{"$set:$event_id"};
    }

    ## calculate ratios, separately in each set, and print
    my $numStages = scalar(@SETS);
    foreach my $set (@SETS) {
      if ($l=~/\tTSS/) {
        $Ratio{$set} = ($TSS_gene{$set}>=0.0000001) ? $FPKM{"$set:$event_id"}/$TSS_gene{$set} : -1;
      } elsif ($l=~/\tTTS/) {
        $Ratio{$set} = ($TTS_gene{$set}>=0.0000001) ? $FPKM{"$set:$event_id"}/$TTS_gene{$set} : -1;
      } elsif ($l=~/\tAE/) {
        $Ratio{$set} = ($AE_gene{$set}>=0.0000001) ? $FPKM{"$set:$event_id"}/$AE_gene{$set} : -1;
      } elsif ($l=~/\tSKIP_ON/) {
        $gene[$i+1] =~ /^(\d+)\tSKIP_OFF\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+$/ or die "illegal line in model.$set:$i: ", $gene[$i];
        my $event_id_2 = $1;
        $Ratio{$set} = ($SKIP_event{"$set:$event_id:$event_id_2"}>=0.0000001) ? 
                       ($FPKM{"$set:$event_id"}/(1.0*$SKIP_event{"$set:$event_id:$event_id_2"})) : -1;
      } elsif ($l=~/\tMSKIP_ON/) {
        $gene[$i+1] =~ /^(\d+)\tMSKIP_OFF\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+$/ or die "illegal line in model.$_";
        my $event_id_2 = $1;
        $Ratio{$set} = ($MSKIP_event{"$set:$event_id:$event_id_2"}>=0.0000001) ?
                       $FPKM{"$set:$event_id"}/$MSKIP_event{"$set:$event_id:$event_id_2"} : -1;
      } elsif ($l=~/\tIR_ON/) {
        $gene[$i+1] =~ /^(\d+)\tIR_OFF\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+$/ or die "illegal line in model.$_";
        my $event_id_2 = $1;
        $Ratio{$set} = ($IR_event{"$set:$event_id:$event_id_2"}>=0.0000001) ?
                       $FPKM{"$set:$event_id"}/$IR_event{"$set:$event_id:$event_id_2"} : -1;
      } elsif ($l=~/\tMIR_ON/) {
        $gene[$i+1] =~ /^(\d+)\tMIR_OFF\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+$/ or die "illegal line in model.$_";
        my $event_id_2 = $1;
        $Ratio{$set} = ($MIR_event{"$set:$event_id:$event_id_2"}>=0.0000001) ?
                       $FPKM{"$set:$event_id"}/$MIR_event{"$set:$event_id:$event_id_2"} : -1;
      }

      if ($Ratio{$set}>=0) { printf "\t$set:%1.6f", $Ratio{$set}; } else { print "\t$set:NA"; $numStages--; }
    }

    # now print out differences between pairs of tissues; this will need to be augmented with
    # some statistical test for multiple samples
    for (my $i=0; $i<scalar(@SETS); $i++) {
      for (my $j=$i+1; $j<scalar(@SETS); $j++) {
        print "\t" . $SETS[$i] . "-" . $SETS[$j] . ":";
        if ($Ratio{$SETS[$i]}>=0 && $Ratio{$SETS[$j]}>=0) { 
           print $Ratio{$SETS[$i]}-$Ratio{$SETS[$j]};
        } else {
           print "NA";
        }
      }
    }
  
    unlink %Ratio;

    print "\t$numStages\n"; 
  }  # end process last gene

}

