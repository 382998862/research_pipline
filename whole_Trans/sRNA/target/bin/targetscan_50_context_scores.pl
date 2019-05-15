#!/usr/local/bin/perl -w

#######################################################################
#
# Copyright(c) 2007,2008 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Author: Joe Rodriguez, Robin Ge, Kim Walker, and George Bell
#         Bioinformatics and Research Computing
#         wibr-bioinformatics@wi.mit.edu
#
# Version: 5.0
#
# Comment: This program calculates TargetScan context scores.
#          It produces output as displayed in TargetScan Release 5.0
#
# This code is available from http://www.targetscan.org/cgi-bin/targetscan/data_download.cgi?db=vert_50
#
#
#######################################################################

# Basic ideas:
#
# 1 - Get names of all input files: 
#		- predicted targets (output from targetscan_50.pl)
#		- UTRs (used by targetscan_50.pl to predicted targets)
#		- mature miRNAs with miRNA name, famly name, seed region, species ID, and mature sequence
# 2 - Get mature miRNAs and UTRs into memory
# 3 - Read through predicted targets file, and for each site 
#		- get all miRNAs in that miRNA family and species [if there are none, we can't calculate context score]
#		- extract short subsequence from UTR to use for predicted consequential pairing
#		- do and score alignment 
#		- extract long subsequence from UTR and get local AU contribution
#		- get position contribution
#		- sum 4 contributions to get total context score
#		- and print out these data
# 4 - Once all context scores have been calculated, go back and calculate context score percentiles.
# 5 - Print out all data, including percentile ranks, into final file

########################  Constants  ########################

# Needed for percentile rank
use POSIX qw(floor ceil);

# link human-readable site types to site type code
$siteTypeDescription2num{'8mer'} = 3;
$siteTypeDescription2num{'m8'} = 2;
$siteTypeDescription2num{'1a'} = 1;
$siteTypeNum2description{3} = "8mer";
$siteTypeNum2description{2} = "m8";
$siteTypeNum2description{1} = "1a";

# Minimum distance to end of CDS.  If site is closer, context scores are not calculated
$MIN_DIST_TO_CDS = 15;
$TOO_CLOSE_TO_CDS = "too_close";
$DESIRED_UTR_ALIGNMENT_LENGTH = 23;

# Number of digits after decimal point for context scores
$DIGITS_AFTER_DECIMAL = 3;

# Print a dot every this number of input lines (so we can see progress)
$DOT_EVERY_THIS_NUM_LINES = 50;

$siteType2siteTypeContribution{3} = -0.310;
$siteType2siteTypeContribution{2} = -0.161;
$siteType2siteTypeContribution{1} = -0.099;

###
###  Regression info from Grimson et al., 2007 supp data
###

$siteType2localAUcontributionSlope{3} = -0.64;
$siteType2localAUcontributionSlope{2} = -0.5;
$siteType2localAUcontributionSlope{1} = -0.42;

$siteType2localAUcontributionIntercept{3} = 0.365;
$siteType2localAUcontributionIntercept{2} = 0.269;
$siteType2localAUcontributionIntercept{1} = 0.236;

$siteType2positionContributionSlope{3} = 0.000172;
$siteType2positionContributionSlope{2} = 0.000091;
$siteType2positionContributionSlope{1} = 0.000072;

$siteType2positionContributionIntercept{3} = -0.07;
$siteType2positionContributionIntercept{2} = -0.037;
$siteType2positionContributionIntercept{1} = -0.032;

$siteType2threePrimePairingContributionSlope{3} = -0.0041;
$siteType2threePrimePairingContributionSlope{2} = -0.031;
$siteType2threePrimePairingContributionSlope{1} = -0.0211;

$siteType2threePrimePairingContributionIntercept{3} = 0.011;
$siteType2threePrimePairingContributionIntercept{2} = 0.067;
$siteType2threePrimePairingContributionIntercept{1} = 0.046;

# Make a ceiling for the distance between a site and the end of the UTR (for position contribution)
$maxDistToNearestEndOfUTREnd = 1500;

########################  Beginning of real code  ########################

getUsage();
getFileFormats();
checkArguments();

# Get mature miRNA data
readMiRNAs();

# Read file of aligned UTRs
readUTRs();

# Get a set of empty strings needed for correction of pairing
getReplaceLength(20);

# Open file for output
open (CONTEXT_SCORES_OUTPUT_TEMP, ">$contextScoresOutputTemp") || die "Cannot open $contextScoresOutputTemp for writing: $!";

# Read file of targets predcited by targetscan_50.pl
# and get context score(s) for sites line by line
readTargets();

# Read all context scores we calculated and get percentile ranks
readContextScoresFromFiles("$contextScoresOutputTemp");
getPercentileRanks();
# Read all data again, add percentile ranks, and print data to final file
addPercentileRanksToOtherData($contextScoresOutputTemp, $contextScoreFileOutput);

unlink($contextScoresOutputTemp);
print STDERR "Deleted temporary file ($contextScoresOutputTemp)\n";
print STDERR "\nAll done!  -  See $contextScoreFileOutput\n";

########################  Subroutines  ########################

sub readTargets
{
	# Read file with info about mature miRNAs, and get this in memory
	
	my $lineNum = 0;
	
	print STDERR "Reading targets file ($predictedTargetsFile)...";
	
	open (PREDICTED_TARGETS, $predictedTargetsFile) || die "Cannot open $predictedTargetsFile for reading: $!";
	while (<PREDICTED_TARGETS>)
	{
		# Gene_ID	miRNA_family_ID	species_ID	MSA_start	MSA_end	UTR_start	UTR_end	Group_num	Site_type	miRNA in this species	Group_type	Species_in_this_group	Species_in_this_group_with_this_site_type
		# CDC2L6	miR-1/206	9606	2827	2833	2588	2594	1	1a	x	1a	9606
		# FNDC3A	let-7/98	10090	2329	2335	2024	2030	31	1a	x	1a	10090 10116 9031 9606 9615
	
		chomp;
		my @f = split (/\t/, $_);

		my $transcriptID = $f[0];
		my $miRNA_familyID = $f[1];
		my $speciesID = $f[2];

		my $utrStart = $f[5];	# UTR start
		my $utrEnd = $f[6];	# UTR end
		my $groupNum = $f[7];

		# print "$transcriptID :: $miRNA_familyID :: $f[8]\n";

		# Type of individual site (not always the same as the group type)
		my $siteType = $siteTypeDescription2num{$f[8]};

		my $familySpecies = "$miRNA_familyID\t$speciesID";		

		if ( $mature_seq{$familySpecies} )
		{
			##  If there are annotated miRNAs in this miRNA family in this species,
			##  extract piece of UTR that we'll need to align with the mature miRNA

			my $subseqForAlignment = extractSubseqForAlignment($transcriptID, $miRNA_familyID, $speciesID, $utrStart, $utrEnd, $siteType);

			###
			###  Extract long subsequence from UTR and get local AU contribution
			###

			my $localAUcontribution = getLocalAUcontribution($transcriptID, $speciesID, $siteType, $utrStart, $utrEnd);

			###
			###  Get position contribution
			###

			my $positionContribution = getPositionContribution($transcriptID, $speciesID, $siteType, $utrStart, $utrEnd);

			##
			##  Get list of miRNAs for this miRNA family + species
			##

			##
			##  Get all miRNAs (if any) for this family and this species			
			##

			for (my $i = 0; $i < $#{$mature_seq{$familySpecies}} + 1; $i++)
			{
				my $thisMiRNA = @{$mature_seq{$familySpecies}}[$i];

				# print "$familySpecies ==> $thisMiRNA + UTR ($subseqForAlignment)\n";

				my ($matureMiRNAid, $matureMiRNA) = split /\t/, $thisMiRNA;

				# Modify $subseqForAlignment based on length of mature miRNA
				my ($finalSubseqForAlignment, $matureMiRNAForAlignment) = 
					modifySubseqForAlignment($matureMiRNA, $matureMiRNAid, $subseqForAlignment, $siteType);

				# print "$miRNA_familyID ($matureMiRNAid) :: $transcriptID :: $speciesID\n";
				# print "$subseqForAlignment\n$matureMiRNAForAlignment\n\n";

				###
				###  Predict consequential pairing (alignment) and get 3' pairing contribution
				###

				my ($threePrimePairingContribution, $slignedUTR, $alignmentBars, $alignedMatureMiRNA, $threePrimePairingScore) = get3primePairingContribution($siteType, $finalSubseqForAlignment, $matureMiRNAForAlignment);
				
				# This site is too close to the CDS: don't take context scores for real
				if ( $utrStart < $MIN_DIST_TO_CDS )
				{
					$threePrimePairingContribution = 
					$localAUcontribution = 
					$positionContribution = 
					$totalcontextScore = $TOO_CLOSE_TO_CDS;
				}
				else
				{
					###
					###  Sum 4 contributions to get total context score
					###

					$totalcontextScore = $siteType2siteTypeContribution{$siteType} + $localAUcontribution + $positionContribution + $threePrimePairingContribution;
				}

				###  
				###  Print out results for this miRNA  
				###  

				print CONTEXT_SCORES_OUTPUT_TEMP "$transcriptID";
				print CONTEXT_SCORES_OUTPUT_TEMP "\t$speciesID";
				print CONTEXT_SCORES_OUTPUT_TEMP "\t$matureMiRNAid";
				print CONTEXT_SCORES_OUTPUT_TEMP "\t$siteTypeNum2description{$siteType}";
				print CONTEXT_SCORES_OUTPUT_TEMP "\t$utrStart";
				print CONTEXT_SCORES_OUTPUT_TEMP "\t$utrEnd";
				print CONTEXT_SCORES_OUTPUT_TEMP "\t$threePrimePairingContribution";
				print CONTEXT_SCORES_OUTPUT_TEMP "\t$localAUcontribution";
				print CONTEXT_SCORES_OUTPUT_TEMP "\t$positionContribution";
				print CONTEXT_SCORES_OUTPUT_TEMP "\t$totalcontextScore";
				print CONTEXT_SCORES_OUTPUT_TEMP "\tNA";
				print CONTEXT_SCORES_OUTPUT_TEMP "\t$slignedUTR";
				print CONTEXT_SCORES_OUTPUT_TEMP "\t$alignmentBars";
				print CONTEXT_SCORES_OUTPUT_TEMP "\t$alignedMatureMiRNA";
				print CONTEXT_SCORES_OUTPUT_TEMP "\t$miRNA_familyID";
				print CONTEXT_SCORES_OUTPUT_TEMP "\t$groupNum";
				print CONTEXT_SCORES_OUTPUT_TEMP "\n";
			}
		}
		else
		{
			# This is a site predicted by comparative genomics and there is no annotated miRNA in this famly in this species
			# Ignore this site, as context score is irrelevant 
		}			
	}
	# Print a dot every $DOT_EVERY_THIS_NUM_LINES input lines
	if ($lineNum % $DOT_EVERY_THIS_NUM_LINES == 0)
	{
		print STDERR ".";
	}
	
	close (PREDICTED_TARGETS);
	close (CONTEXT_SCORES_OUTPUT_TEMP);
	
	print STDERR " done\n";
}


sub readUTRs
{
	# Read file with info about mature miRNAs, and get this in memory
	
	print STDERR "Reading UTRs file ($UTRfile)...";
	
	open (ALIGNED_UTRS, $UTRfile) || die "Cannot open $UTRfile for reading: $!";
	while (<ALIGNED_UTRS>)
	{
		# BMP8B	9606	GUCCACCCGCCCGGC
		# BMP8B	9615	-GUG--CUGCCCACC
	
		chomp;
		my @f = split (/\t/, $_);
		
		my $transcriptID = $f[0];
		my $speciesID = $f[1];
		if($speciesID eq "0") { next; } 	# sequences with speciesID = 0 are consensus; ignore them
		my $seq = $f[2];
		
		# Remove gaps from alignment
		$seq =~ s/-//g;
		
		if ($seq)	# Skip any sequences that are all gaps
		{
			# Convert T's to U's
			$seq =~ s/T/U/g;
			$seq =~ s/t/u/g;

			my $geneIDspecies = "$transcriptID\t$speciesID";
			
			$utr_seq{$geneIDspecies} = $seq;
		}
	}
	close (ALIGNED_UTRS);
	
	print STDERR " done\n";
}

sub readMiRNAs
{
	# Read file with info about mature miRNAs, and get this in memory

	print STDERR "Reading miRNA file ($miRNAfile)...";
	
	open (MATURE_MIRNAS, $miRNAfile) || die "Cannot open $miRNAfile for reading: $!";
	while (<MATURE_MIRNAS>)
	{
		# miR-1/206	GGAAUGU	9606	hsa-miR-1	UGGAAUGUAAAGAAGUAUGUAU
		# let-7/98	GAGGUAG	10090	mmu-let-7a	UGAGGUAGUAGGUUGUAUAGUU
	
		chomp;
		my @f = split (/\t/, $_);
		
		# Make key of miRNA family seed region and species ID
		my $miRNA_familyID = $f[0];
		my $speciesID = $f[1];
		my $familySpecies = "$miRNA_familyID\t$speciesID";
		
		# Make hash linking $familySpecies to mature miRNA sequence + mature miRNA name 
		my $matureMiRNAname = $f[2];
		my $matureMiRNAseq = $f[3];
		my $matureMiRNAnameSeq = "$matureMiRNAname\t$matureMiRNAseq";

		if ( $mature_seq{$familySpecies} )
		{
			push (@{$mature_seq{$familySpecies}}, $matureMiRNAnameSeq);
		}
		else
		{
			@{$mature_seq{$familySpecies}} = $matureMiRNAnameSeq;
		}
	}
	close (MATURE_MIRNAS);
	
	print STDERR " done\n";
}

sub getUsage
{
	$usage = <<EODOCS;

	Description: Calculate context scores predicted miRNA targets
		     using TargetScan methods. 

	USAGE:
		$0 miRNA_file UTR_file PredictedTargets_file ContextScoresOutput_file

	EXAMPLE:
		$0 miR_for_context_scores.txt UTR_Sequences_sample.txt targetscan_50_output.txt context_scores_50_output.txt

	Required input files:
		miRNA_file       => mature miRNA data [not used by targetscan_50.pl]
		UTR_file         => aligned UTRs (same as for targetscan_50.pl)
		PredictedTargets => output from targetscan_50.pl

	Output file:
		ContextScoresOutput => Lists context scores and contributions

	For a description of input file formats, type
		$0 -h

	Authors: Joe Rodriguez, Robin Ge, Kim Walker, and George Bell,
	         Bioinformatics and Research Computing
	Version: 5.0 
	Copyright (c) The Whitehead Institute of Biomedical Research 

EODOCS
}

sub checkArguments
{
	# Check for input and output file arguments
	if ($ARGV[0] && $ARGV[0] eq "-h")
	{
		print STDERR "$usage";
		print STDERR "$fileFormats";
		exit (0);
	}
	elsif (! $ARGV[2])
	{
		print STDERR "$usage";
		exit (0);
	}
	elsif (! -e $ARGV[0])	# miRNA file not present
	{
		print STDERR "\nI can't find the file $ARGV[0]\n";
		print STDERR "which should contain the miRNA families by species.\n";
		exit;
	}
	elsif (! -e $ARGV[1])	# UTR file not present
	{
		print STDERR "\nI can't find the file $ARGV[1]\n";
		print STDERR "which should contain the aligned UTRs.\n";
		exit;
	}
	elsif (! -e $ARGV[2])	# PredictedTargets file not present
	{
		print STDERR "\nI can't find the file $ARGV[2]\n";
		print STDERR "which should contain the predicted targets file from targetscan_50.pl.\n";
		exit;
	}
	
	# Get the file names
	$miRNAfile = $ARGV[0];
	$UTRfile = $ARGV[1];
	$predictedTargetsFile = $ARGV[2];
	$contextScoreFileOutput = $ARGV[3];
	$contextScoresOutputTemp = "$ARGV[3].tmp";
	
	if (-e $contextScoreFileOutput)
	{
		print STDERR "Should I over-write $contextScoreFileOutput [yes/no]? ";
		$answer = <STDIN>;
		if ($answer !~ /^y/i)	{ exit; }
	}
}

sub getFileFormats
{
	$fileFormats = <<EODOCS;

	** Required input files:
	
	1 - miRNA_file    => mature miRNA information
		
		contains four fields (tab-delimited):
		Family members	Seed+m8	Species ID	miRNA_ID	Mature sequence

			a. miRNA family ID/name
			b. species ID in which this miRNA has been annotated
			c. ID for this mature miRNA
			d. sequence of this mature miRNA
			
		ex:	   
		miR-1/206	9606	hsa-miR-1	UGGAAUGUAAAGAAGUAUGUAU
		let-7/98	10090	mmu-let-7a	UGAGGUAGUAGGUUGUAUAGUU
		
	2 - UTR_file      => Aligned UTRs [same format as for targetscan_50.pl]		

		contains three fields (tab-delimited):
			a. Gene/UTR ID or name
			b. Species ID for this gene/UTR (must match ID in miRNA file)
			c. Aligned UTR or gene (with gaps from alignment)
		ex:
		BMP8B	9606	GUCCACCCGCCCGGC
		BMP8B	9615	-GUG--CUGCCCACC
		
		A gene will typically be represented on multiple adjacent lines.	

	3 - PredictedTargets_file => targets predicted by targetscan_50.pl
	
		contains 13 fields, although some fields are ignored
		
		See targetscan_50_output.txt for sample

EODOCS
}

sub extractSubseqForAlignment
{
	# Extract a subsequence of a UTR to use for prediced consequential pairing

	my ($transcriptID, $miRNA_familyID, $speciesID, $utrStart, $utrEnd, $siteType) = @_;

	my $startEnd = "$utrStart\t$utrEnd";

	#####
	# Step 1 - Calculate coordinates needed to extract utr_seq
	#####

	# depending on the seed_type
	# you want to do different things to the start and end coordinates
	# these start and end coordinates refer to which piece of msa_sequence you want for the RNAfold process

	my $real_start = $utrStart - 16;
	if ($real_start < 0)	{ $real_start = 0; }
	my $real_end = $utrEnd + 1;
	if ($real_start >= $real_end)	{ $real_start = 0; }

	#####
	# Step 2 - Extract utr_seq
	#####

	# Now that we have the right coordinates
	# lets get the length of the seq to be extracted (substr needs length)
	# use $key to pull out the right seq from %utr_seq $key = symbol_tax_id
	
	my $length = $real_end - $real_start;
	my $transcriptIDspecies = "$transcriptID\t$speciesID";

	my $seq = $utr_seq{$transcriptIDspecies};
	my $subseqForAlignment = substr($seq, $real_start, $length);

	# one last thing, we need to double check the length of the sequence
	# If the sequence doesn't equal $length ($length is how long the sequence should be)
	# that means it was either at the very end or the very beginning of the utr_seq
	# so the rest of the script runs correctly, add the correct number of N's till it is $length
	
	my $subseqForAlignmentLength = length($subseqForAlignment);

	$subseqSpacer = "";
	
	if ($DESIRED_UTR_ALIGNMENT_LENGTH > $subseqForAlignmentLength)
	{
		for($i = $subseqForAlignmentLength; $i<= $DESIRED_UTR_ALIGNMENT_LENGTH; $i++)
		{
			$subseqSpacer .= "N";
		}
		$subseqForAlignment = "$subseqSpacer$subseqForAlignment";
	}

	# Site is right at the end of the UTR
	# Remove leading Ns and add trailing Ns

	if (length($seq) < ($real_start + $length))
	{
		my $lengthDiff = $real_start + $length - length($seq);

		for($i = 0; $i < $lengthDiff; $i++)
		{
			$subseqForAlignment .= "N";
			$subseqForAlignment =~ s/^NN/  /;
		}
	}
	
	# print STDERR "*$transcriptID, $miRNA_familyID, $speciesID, $utrStart, $utrEnd, $siteType* ==> $subseqForAlignment\n";
	
	return $subseqForAlignment;
}

sub modifySubseqForAlignment
{
	# Modify $subseqForAlignment based on length of mature miRNA

	my ($matureMiRNA, $matureMiRNAid, $subseqForAlignment, $siteType) = @_;
	my $subseqForAlignmentLength = length($subseqForAlignment);

	my $spacer1Length = length($matureMiRNA) - $DESIRED_UTR_ALIGNMENT_LENGTH;
	my $spacer1 = "";
	my $spacer2 = "";
	my $spacer2Length;

	for (my $i = 0; $i < $spacer1Length; $i++)
	{
		$spacer1 .= " ";
	}

	if ($spacer1Length < 0)
	{
		$spacer2Length = -$spacer1Length;
	}
	else
	{
		$spacer2Length = 0;
	}

	for (my $i = 0; $i < $spacer2Length; $i++)
	{
		$spacer2 .= " ";
	}

	if ($DESIRED_UTR_ALIGNMENT_LENGTH > $subseqForAlignmentLength)
	{
		for(my $i = $subseqForAlignmentLength; $i <= $DESIRED_UTR_ALIGNMENT_LENGTH; $i++)
		{
			$spacer2 .= " ";
		}
	}
		############

	if ($siteType == 1) { chop $spacer2; }

	my $finalSubseqForAlignment = $spacer1;
	$finalSubseqForAlignment .= $subseqForAlignment;

	my $matureMiRNAForAlignment = $spacer2;
	$matureMiRNAForAlignment .= reverse $matureMiRNA;

	my @forAlignment;
	push @forAlignment, $finalSubseqForAlignment, $matureMiRNAForAlignment;
	
	return @forAlignment;
}

sub getLocalAUcontribution
{
	my ($transcriptID, $speciesID, $siteType, $utrStart, $utrEnd) = @_;

	my $totalUpScore = 0;
	my $totalDownScore = 0;
	my $utrUpStart;
	my $localAUcontributionMaxRaw = 0;

	my $geneIDspecies = "$transcriptID\t$speciesID";
	$utrSeq = $utr_seq{$geneIDspecies};
	$utrSeqLength = length($utrSeq);

	my $utr5 = $utrStart - 1;
	my $utr3 = $utrSeqLength - $utrEnd;
	
	
	if ($siteType == 1)
	{
		$utrUpStart = $utrStart - 32;
	}
	elsif ($siteType == 2)
	{
		$utrUpStart = $utrStart - 31;
	}
	elsif ($siteType == 3)
	{
		$utrUpStart = $utrStart - 31;
	}
	else
	{ print STDERR "Site type $siteType is not recognized\n"; }

	my $utrDownStart = $utrEnd;
	
	# Get 30 nt upstream of site
	my $utrUp = substr($utrSeq, $utrUpStart, 30);
	# Get 30 nt downstream of site
	my $utrDown = substr($utrSeq, $utrDownStart, 30);
	
	# print "UP:   $utrUp\n";
	# print "DOWN: $utrDown\n";


	if ( length($utrUp) < 30)
	{
		# Too close to 5' (CDS) end of UTR to get all subsequence
		# print STDERR "$transcriptID, $speciesID, $siteType, $utrStart, $utrEnd ==> too close to 5' end; getting score that we can\n";
	}
	elsif ( length($utrDown) < 30)
	{
		# Too close to 3' end of UTR to get all subsequence
		# print STDERR "$transcriptID, $speciesID, $siteType, $utrStart, $utrEnd ==> too close to 5' end; getting score that we can\n";
	}

	my @upScores = ();
	my @downScores = ();

	###
	### 8mer sites
	###
		
	if ($siteType eq '3') 
	{
		# Make an array starting with the position right before the site
		# Make it all uppercase too
		@utrUp3to5 = split(//, reverse (uc ($utrUp)));
		$utrUpLength = @utrUp3to5;
		
		# Get upstream score
		for (my $i = 0; $i <= $#utrUp3to5; $i++)
		{
			$scoreThisPos = 1 / ($i + 1);
			if (($utrUp3to5[$i] eq 'U') || ($utrUp3to5[$i] eq 'A'))
			{
				$totalUpScore += $scoreThisPos;
				push @upScores, $scoreThisPos;
			}
			$localAUcontributionMaxRaw += $scoreThisPos;
		}
	
		@utrDown5to3 = split(//, uc ($utrDown));
		$utrDownLength = @utrDown5to3;
		
		# Get downstream score
		for ($i = 0; $i <= $#utrDown5to3; $i++)
		{
			$scoreThisPos = 1 / ($i + 2);
			if (($utrDown5to3[$i] eq 'U') || ($utrDown5to3[$i] eq 'A'))
			{
				$totalDownScore += $scoreThisPos;
				push @downScores, $scoreThisPos;
			}
			$localAUcontributionMaxRaw += $scoreThisPos;
		}
	}

	if ($siteType eq '2') 
	{
		# Make an array starting with the position right before the site
		# Make it all uppercase too
		@utrUp3to5 = split(//, reverse (uc ($utrUp)));
		$utrUpLength = @utrUp3to5;
		
		# Get upstream score
		for ($i = 0; $i <= $#utrUp3to5; $i++)
		{
			$scoreThisPos = 1 / ($i + 1);
			if (($utrUp3to5[$i] eq 'U') || ($utrUp3to5[$i] eq 'A'))
			{
				$totalUpScore += $scoreThisPos;
				push @upScores, $scoreThisPos;
			}
			$localAUcontributionMaxRaw += $scoreThisPos;
		}
	
		@utrDown5to3 = split(//, uc ($utrDown));
		$utrDownLength = @utrDown5to3;
		
		# Get downstream score
		for ($i = 0; $i <= $#utrDown5to3; $i++)
		{
			$scoreThisPos = 1 / ($i + 1);
			if ($i == 0) { $scoreThisPos = 1 / 2; }

			if (($utrDown5to3[$i] eq 'U') || ($utrDown5to3[$i] eq 'A'))
			{
				$totalDownScore += $scoreThisPos;
				push @downScores, $scoreThisPos;
			}
			$localAUcontributionMaxRaw += $scoreThisPos;
		}
	}	

	if ($siteType eq '1') 
	{
		# Make an array starting with the position right before the site
		# Make it all uppercase too
		@utrUp3to5 = split(//, reverse (uc ($utrUp)));
		$utrUpLength = @utrUp3to5;
		
		# Get upstream score
		for ($i = 0; $i <= $#utrUp3to5; $i++)
		{
			$scoreThisPos = 1 / ($i + 2);
			if (($utrUp3to5[$i] eq 'U') || ($utrUp3to5[$i] eq 'A'))
			{
				$totalUpScore += $scoreThisPos;
				push @upScores, $scoreThisPos;
			}
			$localAUcontributionMaxRaw += $scoreThisPos;
		}
	
		@utrDown5to3 = split(//, uc ($utrDown));
		$utrDownLength = @utrDown5to3;
		
		# Get downstream score
		for ($i = 0; $i <= $#utrDown5to3; $i++)
		{
			$scoreThisPos = 1 / ($i + 2);
			if (($utrDown5to3[$i] eq 'U') || ($utrDown5to3[$i] eq 'A'))
			{
				
				$totalDownScore += $scoreThisPos;
				push @downScores, $scoreThisPos;
			}
			$localAUcontributionMaxRaw += $scoreThisPos;
		}
	}
	
	# Get total score and calculate regression
	$totalLocalAUscoreRaw = $totalUpScore + $totalDownScore ;
	$totalLocalAUscoreFraction = $totalLocalAUscoreRaw / $localAUcontributionMaxRaw;
	
	$totalLocalAUscoreRegression = $totalLocalAUscoreFraction * $siteType2localAUcontributionSlope{$siteType} + $siteType2localAUcontributionIntercept{$siteType};
	
	$totalLocalAUscoreRegression = sprintf("%.${DIGITS_AFTER_DECIMAL}f", $totalLocalAUscoreRegression);
	
	# print "upScores = @upScores\n";
	# print "downScores = @downScores\n";
	
	return $totalLocalAUscoreRegression;
}

sub getPositionContribution 
{
	my ($transcriptID, $speciesID, $siteType, $utrStart, $utrEnd) = @_;

	my $geneIDspecies = "$transcriptID\t$speciesID";
	$utrSeq = $utr_seq{$geneIDspecies};
	$utrSeqLength = length($utrSeq);
	my $distToNearestEndOfUTR;

	my $distTo5primeEndOfUTR = $utrStart - 1;
	my $distTo3primeEndOfUTR = $utrSeqLength - $utrEnd;
	
	if ($distTo5primeEndOfUTR <= $distTo3primeEndOfUTR)
	{
		$distToNearestEndOfUTR = $distTo5primeEndOfUTR;
	}
	else
	{
		$distToNearestEndOfUTR = $distTo3primeEndOfUTR;
	}
	
	# Set max distance to not penalize really long UTRs
	if ($distToNearestEndOfUTR > $maxDistToNearestEndOfUTREnd)
	{
		$distToNearestEndOfUTR = $maxDistToNearestEndOfUTREnd;
	}
		
	my $totalPositionRegression = $distToNearestEndOfUTR * 
		$siteType2positionContributionSlope{$siteType} + 
		$siteType2positionContributionIntercept{$siteType};
		
	$totalPositionRegression = sprintf("%.${DIGITS_AFTER_DECIMAL}f", $totalPositionRegression);
	
	return $totalPositionRegression;
}

sub get3primePairingContribution 
{
	my ($type,$utr,$mirna) = @_;

	$utr =~ s/ //g;
	$mirna =~ s/ //g;
	$utr =~ s/\n//g;
	$mirna =~ s/\n//g;
	$utr = uc $utr;
	$mirna = uc $mirna;
	$utr =~ tr/T/U/;
	$mirna =~ tr/T/U/;
	my %seedinfo = (
		'utrstart'		=>		{ 1 => 8, 2 => 8, 3 => 8 },# 0 based
		'mirnastart'	=>		{ 1 => 7, 2 => 8, 3 => 8 },
		'offset'		=>		{ 1 => 1, 2 => 0, 3 => 1 },
		'overhang'		=>		{ 1 => 1, 2 => 0, 3 => 0 },
		'seedspan'		=>		{ 1 => 6, 2 => 7, 3 => 7 }
		);

	$utr = reverse $utr;
	$mirna = reverse $mirna;

	my($utrNum, $mirnaNum) = (substr($utr, $seedinfo{utrstart}{$type}),substr($mirna, $seedinfo{mirnastart}{$type})); #took part off

	my $maxscore = max(length($utrNum), length($mirnaNum));

	$utrNum =~ tr/AUCGN/12345/;
	$mirnaNum =~ tr/AUCGN/12345/;

	my @UTR = split("", $utrNum);
	my @MIRNA = split("", $mirnaNum);

	my $scorehash;
	my($prevmatch,$tempscore) = (0,0);

	for (my $offset = 0; $offset < $maxscore; $offset++) 
	{
		my $score = 0;
		my $string = "";
		my $tempstring = "";
		my $bestmatch = 0;
		for (my $i = 0; $i <= ($#MIRNA - $offset) && ($i <= $#UTR); $i++)
		{	#impose mirna limit
			if((($UTR[$i] * $MIRNA[$i + $offset]) == 2) || (($UTR[$i] * $MIRNA[$i + $offset]) == 12))
			{
				if(($i + $offset - $seedinfo{overhang}{$type}>= 4) && ($i + $offset - $seedinfo{overhang}{$type} <= 7))
				{
					$tempstring .= "|";
					if($prevmatch == 0)
					{
						$tempscore = 0;
					}
					$tempscore += 1;
				}
				else
				{
					$tempstring .= "|";
					if($prevmatch == 0)
					{
						$tempscore = 0;
					}
					$tempscore += .5;
				}
				$prevmatch++;
			}
			elsif($prevmatch >= 2)
			{
				if($tempscore == $score)
				{
					$string .= $tempstring;
				}
				elsif($tempscore > $score)
				{	# WHICH ONE DO WE TAKE IF EQUAL? don't change score, and leave both...
					$bestmatch = $prevmatch;
					$string =~ s/\|/ /g;
					$string =~ s/X/ /g;
					$string .= $tempstring;
					$score = $tempscore;
				}
				else
				{
					$tempstring =~ s/\|/ /g;
					$tempstring =~ s/X/ /g;
					$string .= $tempstring;
				}
				$string .= " ";
				$tempstring = "";
				$tempscore = 0;
				$prevmatch = 0;
			}
			else
			{
				$tempstring =~ s/\|/ /g;
				$tempstring =~ s/X/ /g;
				$string .= $tempstring;
				$string  .= " ";
				$tempstring = "";
				$tempscore = 0;
				$prevmatch = 0;
			}
		}
		if($prevmatch >= 2)
		{
			if($tempscore == $score)
			{
				$string .= $tempstring;
			}
			elsif($tempscore > $score)
			{
				$bestmatch = $prevmatch;
				$string =~ s/\|/ /g;
				$string =~ s/X/ /g;
				$string .= $tempstring;
				$score = $tempscore;
			}
			$tempscore = 0;
			$prevmatch = 0;
		}
		$score = $score - max(0,(($offset-2)/2));
		$string =~ s/\s([\|X])\s/   /g;
		$string =~ s/^([\|X])\s/  /g;
		$string =~ s/\s([\|X])$/  /g;
		push(@{$scorehash->{$score}}, {offset => $offset,gaploc => 'top',matchstring=>$string});

		$score = 0;
		$tempscore = 0;
		$prevmatch = 0;
		$tempstring = "";
		$bestmatch = 0;
		$string = "";
		for (my $i = 0; ($i <= ($#UTR - $offset)) && ($i <= $#MIRNA); $i++)
		{
			if(($UTR[$i + $offset] * $MIRNA[$i] == 2) || (($UTR[$i + $offset] * $MIRNA[$i]) == 12))
			{	#MATCH
				if(($i - $seedinfo{overhang}{$type}>= 4) && ($i - $seedinfo{overhang}{$type}<= 7))
				{
					$tempstring .= "|";
					if($prevmatch == 0)
					{
						$tempscore = 0;
					}
					$tempscore += 1;
				}
				else
				{
					$tempstring .= "|";
					if($prevmatch == 0)
					{
						$tempscore = 0;
					}
					$tempscore += .5;
				}
				$prevmatch++;
			}
			elsif($prevmatch >= 2)
			{
				if($tempscore == $score)
				{
					$string .= $tempstring;
				}
				elsif($tempscore > $score)
				{	#WHICH ONE DO WE TAKE IF EQUAL? dont change score, and leave both...
					$bestmatch = $prevmatch;
					$string =~ s/\|/ /g;
					$string =~ s/X/ /g;
					$string .= $tempstring;
					$score = $tempscore;
				}
				else
				{
					$tempstring =~ s/\|/ /g;
					$tempstring =~ s/X/ /g;
					$string .= $tempstring;
				}

				$string .= " ";
				$tempstring = "";
				$tempscore = 0;
				$prevmatch = 0;
			}
			else
			{
				$tempstring =~ s/\|/ /g;
				$tempstring =~ s/X/ /g;
				$string .= $tempstring;
				$string  .= " ";
				$tempstring = "";
				$tempscore = 0;
				$prevmatch = 0;
			}
		}
		if($prevmatch >= 2)
		{
			if($tempscore == $score)
			{
				$string .= $tempstring;
			}
			elsif($tempscore > $score)
			{
				$bestmatch = $prevmatch;
				$string =~ s/\|/ /g;
				$string =~ s/X/ /g;
				$string .= $tempstring;
				$score = $tempscore;
			}
			$tempscore = 0;
			$prevmatch = 0;
		}
		$score = $score - max(0,(($offset-2)/2));
		$string =~ s/\s([\|X])\s/   /g;
		$string =~ s/^([\|X])\s/  /g;
		$string =~ s/\s([\|X])$/  /g;
		push(@{$scorehash->{$score}}, {offset => $offset,gaploc => 'bottom',matchstring=>$string});
		$score = 0;
		$tempscore = 0;
		$prevmatch = 0;
		$tempstring = "";
		$bestmatch = 0;
		$string = "";
		}

		#### DONE ALIGNMENT, GET BEST SCORE

		foreach my $score (sort {$b <=> $a} keys %{$scorehash})
		{
			my @outputFields;
			my($i_ret,$offset_ret);
			if(($#{$scorehash->{$score}} == 1) && ($scorehash->{$score}->[0]->{offset} == 0))
			{
				$i_ret = 1;
			}
			elsif(($#{$scorehash->{$score}} > 0))
			{
				#select one with the smallest offset, followed by shortest offset.
				#search each $i, and check offset, if previous min-i is init and greater than current offset, record min $i

				for(my $i = 0; $i <= $#{$scorehash->{$score}}; $i++)
				{
					if(defined($offset_ret))
					{
						if($scorehash->{$score}->[$i]->{offset} < $offset_ret)
						{
							$i_ret = $i;
							$offset_ret = $scorehash->{$score}->[$i]->{offset};
						}
						elsif($scorehash->{$score}->[$i]->{offset} == $offset_ret)
						{
							if(($scorehash->{$score}->[$i_ret]->{gaploc} eq "bottom") && ($scorehash->{$score}->[$i]->{gaploc} eq "bottom"))
							{
								die "ERROR Two tied scores with same offset, and gaplocation";
							}
							elsif($scorehash->{$score}->[$i]->{gaploc} eq "bottom")
							{
								$i_ret = $i;
								$offset_ret = $scorehash->{$score}->[$i]->{offset};
							}
						}
					}
					else
					{
						$i_ret = $i;
						$offset_ret = $scorehash->{$score}->[$i]->{offset};
					}
				}
			}
			else
			{
				#just print the one
				$i_ret = 0;
			}
			my $i = $i_ret;

			#adding the seedmatch check

			my($utrpre, $mirnapre) = (substr($utr, 0,$seedinfo{utrstart}{$type}),substr($mirna, 0,$seedinfo{mirnastart}{$type})); #took part off
			$utrpre =~ tr/AUCGN/12345/;
			$mirnapre =~ tr/AUCGN/12345/;
			my @UTRpre = split("",$utrpre);
			my @MIRNApre = split("",$mirnapre);
			my $matchpre = (" " x ($seedinfo{overhang}{$type} + 1));
			for (my $i = 1; $i <= $seedinfo{seedspan}{$type}; $i++)
			{
				if(($UTRpre[$i + $seedinfo{overhang}{$type}] * $MIRNApre[$i] == 2) || (($UTRpre[$i + $seedinfo{overhang}{$type}] * $MIRNApre[$i]) == 12))
				{	#MATCH
					$matchpre .= "|";
				}
				else
				{
					# !!! Something wrong: print E for the error	
					# Seed region does not show alignment -- and it always should
					$matchpre .= "E";
				}
			}
			$matchpre .= (" " x $scorehash->{$score}->[$i]->{offset});
			my $string = $matchpre . $scorehash->{$score}->[$i]->{matchstring};

			my $mirnapreoff = " " x $seedinfo{overhang}{$type};

			my $utrout = substr($utr, 0,$seedinfo{utrstart}{$type}) . substr($utr, $seedinfo{utrstart}{$type});
			my $mirnaout = $mirnapreoff . substr($mirna, 0,$seedinfo{mirnastart}{$type}) . substr($mirna, $seedinfo{mirnastart}{$type});

			my $offsetstring = "-" x $scorehash->{$score}->[$i]->{offset};
			if ($scorehash->{$score}->[$i]->{gaploc} eq "top")
			{
				$utrout = substr($utr, 0,$seedinfo{utrstart}{$type}) . "$offsetstring" . substr($utr, $seedinfo{utrstart}{$type});
			}
			elsif($scorehash->{$score}->[$i]->{gaploc} eq "bottom")
			{
				$mirnaout = $mirnapreoff . substr($mirna, 0,$seedinfo{mirnastart}{$type}) . "$offsetstring" . substr($mirna, $seedinfo{mirnastart}{$type});
			}
			my $longest = max(length($utrout),length($mirnaout));
			my $postutr =  " " x ($longest - length($utrout));
			my $postmatch = " " x ($longest - length($string));
			my $postmirna = " " x ($longest - length($mirnaout));
			$utrout .= $postutr;
			$string .= $postmatch;
			$mirnaout .= $postmirna;
			$utrout = reverse($utrout);	# UTR with gaps
			$string = reverse($string); # alignmend bars
			$mirnaout = reverse($mirnaout); # miRNA with gaps

			if ($score < 3)
			{
				# Don't show consequential pairing if this score < 3
				# This is as of 17 Oct 2008

				($utrout, $string, $mirnaout) = getConsequentialPairing($utrout, $string, $mirnaout);
			}

			my $threePrimePairingRegression = $score * 
				$siteType2threePrimePairingContributionSlope{$type} + 
				$siteType2threePrimePairingContributionIntercept{$type};

			# Prepare output including UTR and miRNA with gaps where appropriate
			push @outputFields, $threePrimePairingRegression, $utrout, $string, $mirnaout, $score;

			return @outputFields;
			last;
		}
		return (0,0,0);
}

sub max
{
	my @list = @_;
	my $maximum = shift @list;
	foreach (@list)
	{
		if($maximum < $_)
		{
			$maximum = $_;
		}
	}
	return $maximum;
}

sub getConsequentialPairing
{
	# Don't show consequential pairing if raw 3' pairing score < 3

	my ($utrSeq, $pairing, $mirnaSeq) = @_;
	my $newSubstring;
	my $utrSeqGapless;
	my $mirnaSeqGapless;
	
	# Working with pipes is a pain
	$pairing =~ s/\|/x/g;
	
	# If there's 1 - 2 bits of complementarity in the 3' end of the miRNA-target pair,
	# remove it, replacing it with spaces
	
	if ($pairing =~ /(\S+\s+\S+)\s+(\S+)/ || $pairing =~ /(\S+)\s+(\S+)/)
	{
		my $pairingToRemove = $1;
		my $pairingToRemoveLength = length($pairingToRemove);

		if ($replaceLength{$pairingToRemoveLength})
		{
			my $newSubstring = $replaceLength{$pairingToRemoveLength};
			
			# print "replaceLength pairingToRemoveLength ==> \"$newSubstring\"\n";		
			# print "Need to replace $pairingToRemove ($pairingToRemoveLength chars) from \"$pairing\" and replace with \"$newSubstring\"\n";

			$pairing =~ s/$pairingToRemove/$newSubstring/;
			
			###
			###  Also, remove all gaps from UTR seq and miRNA seq
			###
			
			$utrSeq = removeGapAddLeadingSpace($utrSeq);
			$mirnaSeq = removeGapAddLeadingSpace($mirnaSeq);
		}
	}
	
	# change back to pipes
	$pairing =~ s/x/|/g;

	return ($utrSeq, $pairing, $mirnaSeq);
}

sub removeGapAddLeadingSpace
{
	my $seq = $_[0];
	
	my $seqLengthPre = length($seq);
	$seq =~ s/-//g;
	my $seqLengthPost = length($seq);
	
	# How many gaps were removed?
	my $seqLengthDiff = $seqLengthPre - $seqLengthPost;

	# Add leading spaces to seq if gap was removed
	if ($seqLengthDiff > 0)
	{
		my $seqGapless = "";

		for ($i = 1; $i <= $seqLengthDiff; $i++)
		{
			$seqGapless .= " ";
		}
		$seqGapless .= $seq;
		$seq = $seqGapless;
	}
	return $seq;
}

sub readContextScoresFromFiles
{
	my @contextScoreFiles = @_;

	foreach $contextScoreFile (@contextScoreFiles)
	{
		open (SCORES, $contextScoreFile) || die "Cannot open $contextScoreFile: $!";
		while (<SCORES>)
		{
			chomp;
			my @f = split (/\t/, $_);

			my $miRNA = $f[2];
			my $contextScore = $f[9];

			# check that $contextScore is really a number
			if ($contextScore =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)$/)
			{
				push @{ $miRNAtoContextScores{$miRNA} }, $contextScore;
			}
		}
		close (SCORES);
	}
}

sub getPercentileRanks
{
	foreach $miRNA (sort keys %miRNAtoContextScores)
	{
		# Sort from least negative to more negative
		my @contextScoreList = sort {$b <=> $a} @{ $miRNAtoContextScores{$miRNA} };

		# print "Doing $miRNA ...\n";

		# Initialize
		my $numLower = 0;

		# Go through all the context scores in the list
		for (my $i = 0; $i <= $#contextScoreList; $i++)
		{
			my $scoreListLength = $#contextScoreList + 1;

			if ($contextScoreList[$i - 1] && $contextScoreList[$i] != $contextScoreList[$i - 1])
			{
				# Get the number of sites with a lower score
				$numLower = $i;
			}
			my $percentileRank = 100 * ( $numLower ) / $scoreListLength;

			# Get the floor so we never have a percentile rank of 100
			$percentileRank = floor($percentileRank);

			$miRNAcontextScoreToPctile{$miRNA}{$contextScoreList[$i]} = $percentileRank;

			# print "$miRNA\t$contextScoreList[$i]\t$numLower\t$percentileRank\n";			
		}
	}
}

sub addPercentileRanksToOtherData
{
	my ($contextScoresOutputTemp, $contextScoreFileOutput) = @_;
	
	open (CONTEXT_SCORES_OUTPUT_FINAL, ">$contextScoreFileOutput") || die "Cannot open $contextScoreFileOutput for writing: $!";
	print CONTEXT_SCORES_OUTPUT_FINAL "Gene ID	Species ID	Mirbase ID	Site Type	UTR start	UTR end	3' pairing contribution	local AU contribution	position contribution	context score	context score percentile	UTR region	UTR-miRNA pairing	mature miRNA sequence	miRNA family	Group #\n";

	# Open the file with all the data so far
	open (CONTEXT_SCORES_OUTPUT_TEMP, $contextScoresOutputTemp) || die "Cannot open $contextScoresOutputTemp FOR READING: $!";
	while (<CONTEXT_SCORES_OUTPUT_TEMP>)
	{
		chomp;
		my @f = split (/\t/, $_);

		my $miRNA = $f[2];
		my $contextScore = $f[9];
		my $percentileRank;

		# Check whether this is a miRNA for which we calculated context scores
		# And make sure it's not a "too close to ORF" site
		if ($contextScore =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)$/)
		{
			$percentileRank = $miRNAcontextScoreToPctile{$miRNA}{$contextScore};

			$f[10] = $percentileRank;
		}

		# Add the percentile rank to the rest of the data
		my $lineWithPercentileRank = join "\t", @f;

		print CONTEXT_SCORES_OUTPUT_FINAL "$lineWithPercentileRank\n";
	}
	close(CONTEXT_SCORES_OUTPUT_TEMP);
}

sub getReplaceLength
{
	my $longestString = $_[0];
	
	my $replacementString = "";

	for (my $i = 1; $i <= $longestString; $i++)
	{
		$replacementString .= " ";
	
		$replaceLength{$i} = $replacementString;
	}
}
