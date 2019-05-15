#!/usr/local/bin/perl -w

#######################################################################
# Copyright(c) 2007,2008 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Author: Sanjeev Pillai and George Bell
#         Bioinformatics and Research Computing
#         wibr-bioinformatics@wi.mit.edu
#
# Version: 5.0
#
# Comment: This program predicts miRNA targets using the TargetScanS algorithm.
#          It produces output as displayed in TargetScan Release 5.0
#
# This code is available from http://www.targetscan.org/cgi-bin/targetscan/data_download.cgi?db=vert_50
#
# Version: 1.1 corrects for 8mer 'UTR end' OBOB
# Version: 5.0 (11 Nov 2008)
#          - allows grouping of sites that are 8mer, 7mer-m8 AND 7mer-1a
#          - These are defined as a new type: "8mer+m8+1a"
#          - Group type is now a list of types in the group, rather than what's conserved
#          - An additional column in printed to show species in the group with a species' site type
#          - Handling of memory for large datasets has been improved
#
#######################################################################

# Basic ideas:
#
# 1 - Grab all miRNA info.
# 2 - Read through UTRs, getting those from one gene at a time.
# 3 - Identify miRNA sites for this gene.
# 4 - Group overlapping miRNA sites in different species into a group
#

# Find sites in species even in which the miRNA has not been annotated
# If you don't want to identify these, set the variable to 0.
$FIND_SITES_ALL_SPECIES = 1;
# For site comparison between species, how much of an overlap (number of positions/nt) are required
$REQUIRED_OVERLAP = 2;
# Number of nt at beginning of UTR to mask (since miRNAs don't target right next to CDS)
$BEG_UTR_MASK_LENGTH = 0;
# If $VERBOSE is non-zero, each Gene ID will be printed as it's processed. 
$VERBOSE = 1;

# Show how to summarize groups
getGroupTypeSummaryInfo();

#######################  End of global variables  #######################

# Start program

getUsage();
getFileFormats();
checkArguments();

# Get miRNA family data
readMiRNAs();

# Write conserved group to this file
open (COORDS, ">$coordsFile") || die "Cannot open $coordsFile for writing: $!";
# Print output file header
print COORDS "Gene_ID\tmiRNA_family_ID\tspecies_ID\tMSA_start\tMSA_end\tUTR_start\tUTR_end\tGroup_num\tSite_type\tmiRNA in this species\tGroup_type\tSpecies_in_this_group\tSpecies_in_this_group_with_this_site_type\n";


#######################  Go through UTR file, processsing one gene at a time  #######################


# Initialize
$lastUtrID = "";
# This is the number (ID) of each group of sites; start with 1 and count up.
$groupNum = 0;

open (UTRS, $UTRfile) || die "Cannot open $UTRfile for reading: $!";
while (<UTRS>)
{
	# NM_031304       hg18    -       GGCCCCAC
	
	chomp;
	s/\r//g;	# For Windows and Mac  (When in doubt, convert input files to Unix format)

	if (! /^\s*$/)	# Ignore empty lines
	{
		# Public code format
		($utrID, $thisSpeciesID, $thisUTR) = split (/\t/, $_);

		if ($thisSpeciesID)	# Skip consensus sequence (if present)
		{
			$strand = $strand;	# Ignore this for now

			# Convert from RNA to DNA if needed
			$thisUTR =~ s/T/U/gi;
			
			# Mask beginning of UTR since miRNA can't target UTR right next to CDS
			for (my $i = 0; $i < $BEG_UTR_MASK_LENGTH; $i++)
			{
				$thisUTR =~ s/[ACGTU]/N/i;
			}

			if ($utrID && $utrID ne $lastUtrID)
			{
				if ($VERBOSE)
				{
					print STDERR "Processing $utrID\n";
				}
				$numUTR_IDs++;
				
				# Look for sites in this gene (UTR) and process the results
				processUTRset();

				# Empty out these variables after finishing this set of UTRs
				%groupNumPlusType2speciesList = ();
				%groupNumToGroupType = ();
				%groupNumToSiteTypes = ();
				%groupNumToSiteTypesList = ();
				%groupNumToSpecies = ();
				%groupToInfo = ();
				%pairToDistance = ();
				%siteToGroupNum = ();
				%siteTypeThisSite = ();
				%species2utr = ();
				%speciesStartEnd = ();
				%speciesThisSite = ();
			}

			# Add this UTR to the set
			$species2utr{"$thisSpeciesID"} = $thisUTR;

			$lastUtrID = $utrID; 			
		}
	}
}

# Get the last one
processUTRset();

print STDERR "Finished predicting miRNA targets of $numUTR_IDs genes\n";

############################  Subroutines  ############################


sub readMiRNAs
{
	my ($mirFamID, $mirSeedRegion, $mirSpeciesID);

	open (MIR_FAMILY_DATA, $miRNAfile) || die "Cannot open $miRNAfile for reading: $!";
	while (<MIR_FAMILY_DATA>)
	{
		# let-7/98	GAGGUAG	10090
	
		chomp;
		s/\r//g;	# For Windows and Mac
		
		###################  Public data format  ###################
		
		($mirFamID, $mirSeedRegion, $mirSpeciesID) = split (/\t/, $_);
		
		# Convert from RNA to DNA if needed
		$mirSeedRegion =~ s/T/U/gi;
		
		$mirID2seed{$mirFamID} = $mirSeedRegion;
		
		# Make sure we know which miRNA family is present in each species
		$MIR_ID_SPECIES{"${mirFamID}::$mirSpeciesID"} = 1;
	}
	
	
	# Set patterns for search for (3 types)
	foreach $mirFamID (sort keys %mirID2seed) 
	{
		($miR2match0{$mirFamID} , $miR2match1{$mirFamID}, $miR2match2{$mirFamID}) = get_seeds($mirID2seed{$mirFamID});
	}
}

sub processUTRset
{
	my $speciesIDthisUTR;
	my (@location0, @location1, @location2);
	my (@length0, @length1, @length2);
	my (@number0, @number1, @number2);
	my ($location0, $length0, $number0);
	my ($location1, $length1, $number1);
	my ($location2, $length2, $number2);
	my (@s, @l, @d);

	# Look at each miRNA family
	foreach $mirFamID (sort keys %mirID2seed) 
	{
		# Initialize/reset
		@species2s = @species2l = ();
		@s = @l = @d = ();
	
		# Do one species' UTR at a time
		foreach $speciesIDthisUTR (sort keys %species2utr)
		{
			# Is this miRNA in this species?
			# If so [or if we want to look anyway], look for sites
			
			if ($FIND_SITES_ALL_SPECIES || $MIR_ID_SPECIES{"${mirFamID}::$speciesIDthisUTR"})
			{
				# Look for 7mer-m8 (?) sites for this miRNA
				($location0, $length0, $number0) = getmatches($species2utr{"$speciesIDthisUTR"}, $miR2match0{$mirFamID});
				@location0 = @$location0; @length0 = @$length0; @number0 = @$number0;

				# Look for 7mer-1A (?) sites for this miRNA
				($location1, $length1, $number1) = getmatches($species2utr{"$speciesIDthisUTR"}, $miR2match1{$mirFamID});
				@location1 = @$location1; @length1 = @$length1; @number1 = @$number1;

				# Look for 8mer (?) sites for this miRNA
				($location2, $length2, $number2) = getmatches($species2utr{"$speciesIDthisUTR"}, $miR2match2{$mirFamID});
				@location2 = @$location2; @length2 = @$length2; @number2 = @$number2;

				# Merge these types of sites when possible
				# ($s_ref, $l_ref, $d_ref) = arrcmp(\@location0,\@location1,\@location2,\@length0,\@length1,\@length2,\@number0,\@number1,\@number2);
				($s_ref, $l_ref, $d_ref) = arrcmp(\@location0,\@location1,\@location2,\@length0,\@length1,\@length2,\@number0,\@number1,\@number2);

				# Dereference
				@s = @$s_ref;
				@l = @$l_ref;
				@d = @$d_ref;

				$species2s{$speciesIDthisUTR} = [ @s ]; 
				$species2l{$speciesIDthisUTR} = [ @l ];				
			}
		}

		# If there are any hit(s) for this miRNA in any species, group each orthologous set of sites

		if (%species2s)
		{
			groupSitesThisMiRNA();

			# Finish processing this gene/miR data 
			summarizePrintGroupsThisGeneThisMir();
		}

		# Empty these out
		@outputThisGeneThisMir = ();
		@groupNumThisGeneThisMir = ();		
	}
}

sub groupSitesThisMiRNA
{
	# Initialize
	%speciesStartEnd = ();
	%groupNumToGroupType = ();
	%siteToGroupNum = ();
	my ($spID, @s, @l);
	
	###############  Make a hash (%speciesStartEnd) of all the sites for this miRNA for this gene
	
	foreach $spID (sort keys %species2utr)
	{
		# If there's a hit and this miRNA is annotated in this species
		if ($species2s{$spID} && ($FIND_SITES_ALL_SPECIES || $MIR_ID_SPECIES{"${mirFamID}::$spID"}))
		{
			@s = @{ $species2s{$spID} };
			@l = @{ $species2l{$spID} };
			
			for ($i = 0; $i <= $#s; $i++)
			{
				@posType = split(/\s+/, $s[$i]);

				$end = $posType[0];
				$start = $posType[0] - $l[$i] + 1;
				$type = $posType[1];

				# Make a hash of all sites, linking to type
				$speciesStartEnd{"${spID}::${start}::$end"} = $type;

				# print "groupSitesThisMiRNA =>\t${mirFamID}\t$spID\t${start}-$end\n";
			}
		}
	}
	
	############### Check for position overlap between sites

	# Do an all vs. all comparison to identify overlaps
	foreach $site1 (sort keys %speciesStartEnd)
	{
		foreach $site2 (sort keys %speciesStartEnd)
		{
			@site1Info = split (/::/, $site1);
			@site2Info = split (/::/, $site2);

			$site1Species = $site1Info[0];
			$site1Start = $site1Info[1];
			$site1End = $site1Info[2];

			$site2Species = $site2Info[0];
			$site2Start = $site2Info[1];
			$site2End = $site2Info[2];
			
			# Skip comparison of same-species sites
			if ($site1Species ne $site2Species)
			{
				############  Choose combinations to give overlap  ############

				# Same start and end
				if ($site1Start eq $site2Start && $site1End eq $site2End)
				{
					groupThisPair($site1, $site2);
					
					$pairToDistance{"$site1 $site2"} = $site1End - $site1Start + 1;
				}
				# Same start
				elsif ($site1Start eq $site2Start)
				{
					groupThisPair($site1, $site2);
					
					if ($site1End > $site2End)
					{ $pairToDistance{"$site1 $site2"} = $site1End - $site1Start + 1;}
					else { $pairToDistance{"$site1 $site2"} = $site2End - $site2Start + 1;}
				}			
				# Same end
				elsif ($site1End eq $site2End)
				{
					groupThisPair($site1, $site2);
					
					if ($site1Start < $site2Start)
					{ $pairToDistance{"$site1 $site2"} = $site1End - $site1Start + 1;}
					else { $pairToDistance{"$site1 $site2"} = $site2End - $site2Start + 1;}
					
				}
				# Offset one direction
				#     xxxxxxx
				#    xxxxxxx
				elsif ($site1Start > $site2Start && $site1Start <= $site2End)
				{
					$numOverlapNt = $site2End - $site1Start + 1;
					if ($numOverlapNt >= $REQUIRED_OVERLAP)
					{
						groupThisPair($site1, $site2);
						$pairToDistance{"$site1 $site2"} = $numOverlapNt;
					}				
				}
				# Offset other direction
				#    xxxxxxx
				#     xxxxxxx
				elsif ($site1End >= $site2Start && $site1End < $site2End)
				{
					$numOverlapNt = $site1End - $site2Start + 1;
					if ($numOverlapNt >= $REQUIRED_OVERLAP)
					{
						groupThisPair($site1, $site2);
						$pairToDistance{"$site1 $site2"} = $numOverlapNt;
					}
				}
				# One within the other (with gaps)
				#      xxxxxxx          xxxxxxxxx
				#     xxxxxxxxx          xxxxxxx
				elsif ( ($site1Start > $site2Start && $site1End < $site2End) ||
					($site2Start > $site1Start && $site2End < $site1End) )
				{
					groupThisPair($site1, $site2);
					$pairToDistance{"$site1 $site2"} = $numOverlapNt;
				}				
			}
		}
	}
	
	############  Gather info for all sites for this gene and this miRNA
	
	foreach $thisSite (sort keys %speciesStartEnd)
	{
		@siteAllInfo = split (/::/, $thisSite);
		
		$speciesThisSite = $siteAllInfo[0];
	
		# This site is a group of 1, so no group info yet
		if (! $siteToGroupNum{$thisSite})	
		{
			# If this group hasn't yet been assigned a number, give it one.
			$groupNum++;
			push @groupNumThisGeneThisMir, $groupNum;
			$siteToGroupNum{$thisSite} = $groupNum;
		}
		if (! $groupNumToGroupType{$groupNum})
		{
			# If this group hasn't yet been assigned a type, it's the same as the site
			$groupNumToGroupType{$groupNum} = $speciesStartEnd{$thisSite};
		}

		###########################  5.2.2007

		$groupNumHere = $siteToGroupNum{$thisSite};
		if (! $groupNumToSiteTypes{$groupNumHere})
		{
			# Start a list of site types for this group
			$groupNumToSiteTypes{$groupNumHere} = "$speciesStartEnd{$thisSite}";
		}
		else
		{
			# Add to the list of site types for this group
			$groupNumToSiteTypes{$groupNumHere} .= ";$speciesStartEnd{$thisSite}";
		}
		
		# Make a list of species in which a site type is found (in this group)  
		$groupNumPlusType2speciesList{$groupNumHere}{$speciesStartEnd{$thisSite}} .= "$speciesThisSite ";
		
		if ($speciesStartEnd{$thisSite} eq "8mer")
		{
			$groupNumPlusType2speciesList{$groupNumHere}{"1a"} .= "$speciesThisSite ";
			$groupNumPlusType2speciesList{$groupNumHere}{"m8"} .= "$speciesThisSite ";	
		}
		
		###########################

		# Given the MSA coords, get the corresponding UTR coords
		# Correct 8mer OBOB: Nov 26, 2007
		($utrStart, $utrEnd) = utrcoords($species2utr{$siteAllInfo[0]}, $siteAllInfo[1], $speciesStartEnd{$thisSite});
	
		# Link each group to the species within it
		$groupNumToSpecies{$siteToGroupNum{$thisSite}} .= "$siteAllInfo[0];";
	
		# Is this miRNA annotated in this species? 
		if (! $MIR_ID_SPECIES{"${mirFamID}::$speciesThisSite"})
		{
			$annotated = "";
		}
		else
		{
			$annotated = "x";
		}

		#
		#
		# Format: ABI3BP	miR-103/107	9615	1145	1156	926	932	8404	m8:1a	1a
		
		push @outputThisGeneThisMir, "$lastUtrID\t$mirFamID\t$speciesThisSite\t$siteAllInfo[1]\t$siteAllInfo[2]\t$utrStart\t$utrEnd\t$siteToGroupNum{$thisSite}\t$speciesStartEnd{$thisSite}\t$annotated";

		#
		#
		#
	}
}

sub summarizePrintGroupsThisGeneThisMir
{
	my ($groupNum, $speciesList);
	my (@species,@uniqueSpecies, @speciesThisGroup);

	# Print summary data about species in each group

	foreach $groupNum (sort {$a <=> $b} @groupNumThisGeneThisMir ) 
	{
		$speciesList = $groupNumToSpecies{$groupNum};

		@species = split (/;/, $speciesList);
		@species = sort @species;

		# Make this species list unique (shouldn't be necessary) and sorted
		%uniqueSpeciesThisGroup = ();
		%uniqueSpeciesThisGroup = map { $_ => 1 } @species;
		@uniqueSpecies = keys %uniqueSpeciesThisGroup;
		
		# If species IDs are numbers
		# @speciesThisGroup = sort {$a <=> $b} @uniqueSpecies;
		# If species IDs are not numbers
		@speciesThisGroup = sort @uniqueSpecies;

		###########################  5.2.2007

		if ($groupNumToSiteTypes{$groupNum})
		{
			@siteTypes = split (/;/, $groupNumToSiteTypes{$groupNum});
			
			# Make this site type list unique (shouldn't be necessary) and sorted 5.2.07
			%uniqueSiteTypesThisGroup = ();
			%uniqueSiteTypesThisGroup = map { $_ => 1 } @siteTypes;
			@uniqueSiteTypes = keys %uniqueSiteTypesThisGroup;		
			@uniqueSiteTypes = sort @uniqueSiteTypes;
			
			$uniqueSiteTypesList = join " ", @uniqueSiteTypes;
		}
		else
		{
			$uniqueSiteTypesList = "";
		}
		
		$groupNumToSiteTypesList{$groupNum} = $uniqueSiteTypesList;

		###########################
		
		$groupToInfo{$groupNum} = "@speciesThisGroup";
	}
	
	foreach $dataOneSiteThisGeneThisMir (@outputThisGeneThisMir)
	{
		@f = split (/\t/, $dataOneSiteThisGeneThisMir);
		$groupNumThisSite = $f[7];
		$siteTypeThisSite = $f[8];
		
		###########################  5.2.2007
		
		if ($groupSiteTypesList2GroupType{$groupNumToSiteTypesList{$groupNumThisSite}})
		{
			$groupType = $groupSiteTypesList2GroupType{$groupNumToSiteTypesList{$groupNumThisSite}};
		}
		else	# Something wrong
		{
			print STDERR "Unrecognized group type: $groupNumToSiteTypesList{$groupNumThisSite}\n";
		}
		
		$dataOneSiteThisGeneThisMir .= "\t$groupType";
		$dataOneSiteThisGeneThisMir .= "\t$groupToInfo{$groupNumThisSite}";
		
		# Add the species list for this site type (but only if site type ne group type)
		if ($groupType ne $siteTypeThisSite)
		{
			$dataOneSiteThisGeneThisMir .= "\t$groupNumPlusType2speciesList{$groupNumThisSite}{$siteTypeThisSite}";
		}
		
		print COORDS "$dataOneSiteThisGeneThisMir\n";
	}
}

sub groupThisPair
{
	if (! $siteToGroupNum{$site1} || ! $siteToGroupNum{$site2})
	{
		if ($siteToGroupNum{$site1})	# Site 1 already part of a group
		{
			# Add Site 2 to the group
			$siteToGroupNum{$site2} = $siteToGroupNum{$site1};
		}
		elsif ($siteToGroupNum{$site2})	# Site 2 already part of a group
		{
			# Add Site 1 to the group
			$siteToGroupNum{$site1} = $siteToGroupNum{$site2};
		}
		else
		{
			# Increment the group number
			$groupNum++;
			push @groupNumThisGeneThisMir, $groupNum;
			$siteToGroupNum{$site1} = $groupNum;
			$siteToGroupNum{$site2} = $groupNum;
		}
	}
}

sub getUsage
{
	$usage = <<EODOCS;

	Description: Search for predicted miRNA targets
		     using the modified TargetScanS algorithm. 

	USAGE:
		$0 miRNA_file UTR_file PredictedTargetsOutputFile

	Required input files:
		miRNA_file    => miRNA families by species
		UTR_file      => Aligned UTRs		

	Output file:
		PredictedTargetsOutputFile    => Lists sites using alignment coordinates (MSA and UTR)

	For a description of input file formats, type
		$0 -h

	Authors: George Bell & Sanjeev Pillai, Bioinformatics and Research Computing
	Version: 5.0 
	Copyright (c) The Whitehead Institute of Biomedical Research 

EODOCS
}

sub getFileFormats
{
	$fileFormats = <<EODOCS;

	** Required input files:
	
	1 - miRNA_file    => miRNA families by species
		
		contains three fields (tab-delimited):
			a. miRNA family ID/name
			b. seed region (7mer) for this miRNA
			c. species ID in which this miRNA has been annotated
		ex:	   
		let-7/98	GAGGUAG	10090
		let-7/98	GAGGUAG	9606
		
		A miRNA family that is present in multiple species
		should be represented in multiple lines, one for each species.
		
	2 - UTR_file      => Aligned UTRs		

		contains three fields (tab-delimited):
			a. Gene/UTR ID or name
			b. Species ID for this gene/UTR (must match ID in miRNA file)
			c. Aligned UTR or gene (with gaps from alignment)
		ex:
		BMP8B	9606	GUCCACCCGCCCGGC
		BMP8B	9615	-GUG--CUGCCCACC
		
		A gene will typically be represented on multiple adjacent lines.	

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
		print STDERR "which should contain the Aligned UTRs.\n";
		exit;
	}
	
	$miRNAfile = $ARGV[0];
	$UTRfile = $ARGV[1];
	$coordsFile = $ARGV[2];
	
	if (-e $coordsFile)
	{
		print STDERR "Should I over-write $coordsFile [yes/no]? ";
		$answer = <STDIN>;
		if ($answer !~ /^y/i)	{ exit; }
	}
}

# get_seeds returns the m8 and 1a seed matches for a given seed sequence

sub get_seeds 
{
	my $seedRegion = $_[0];
	
	# get the m8 seed match

	my $seed1 = $seedRegion;
	my $rseed1 = reverse($seed1);
	$rseed1 =~ tr/AUCG/UAGC/; 

	# get the 1a seed match
	
	my $seed2 = $seedRegion;
	$seed2 =~ s/(\w+)\w$/$1/;
	$seed2 =~ s/(\w+)/U$1/;
	my $rseed2 = reverse($seed2);
	$rseed2 =~ tr/AUCG/UAGC/;

	# get m8:1a seed match
	my $rseed0 = $rseed1 . 'A';

	@match0 = split '', $rseed0;
	@match1 = split '', $rseed1;
	@match2 = split '', $rseed2;
			
	# print "$mirFamID\t$mirID2seed{$mirFamID}\t$rseed0\t$rseed1\t$rseed2\n";

	$pattern0 = $pattern1 = $pattern2 = "";
	
	# Add potential gaps between each pair of nts
	
	$count0 = 0;	
	foreach my $mt0 (@match0) {
		$count0++;
		($count0 < 8) ? ($pattern0 .= "$mt0-{0,}") : ($pattern0 .= "$mt0");
	}

	$count1 = 0;	
	foreach my $mt1 (@match1) {
		$count1++;
		($count1 < 7) ? ($pattern1 .= "$mt1-{0,}") : ($pattern1 .= "$mt1");
	}

	$count2 = 0;
	foreach my $mt2 (@match2) {
		$count2++;
		($count2 < 7) ? ($pattern2 .= "$mt2-{0,}") : ($pattern2 .= "$mt2");
	}

	return($pattern0, $pattern1, $pattern2);
}

# get_matches returns the position and lengths of the matches
# in the alignment sequence for a given seed sequence

sub getmatches {
	my ($in,$match) = ($_[0],$_[1]);
	my ($num,$position,@locs,$len,@lens,$bnum,$dnum,$dsh,@dshs);
	$num = 0;

	while ($in =~ m/$match/g) {
		$num++;
		$position =  pos($in);
		push @locs, "$position ";
		$matched = $&;
		$len = length($&);
		($matched =~ m/\-+\w$/)?($dnum=length($&) -1 ):($dnum=0);
		($matched =~ m/^\w(\-+)\w/)?($bnum=length($1)):($bnum=0);
		$dsh = "$bnum-$dnum";
		push @dshs, $dsh;
		push @lens, $len;
		pos($in) -= 5;
	}
	# if (@locs) { print "getmatches =>\t$mirFamID\t@locs\t@lens\t@dshs\n"; }
	
	return (\@locs,\@lens,\@dshs);
}

#####################################################################
# arrcmp compares the 8mer match with the 7mer matches for each alignment-seed 
# match and checks for overlap between the two. The overlapping ones 
# in the 7mer arrays are ignored and all are merged into a single array.
#####################################################################

sub arrcmp 
{
	my ($inp0,$inp1,$inp2,$l0,$l1,$l2,$n0,$n1,$n2) = @_;

	my @input0 =  @$inp0;
	my @input1 =  @$inp1;
	my @input2 =  @$inp2;
	my @len0 = @$l0;
	my @len1 = @$l1;
	my @len2 = @$l2;
	my @num0 = @$n0;
	my @num1 = @$n1;
	my @num2 = @$n2;

	my ($val0,$val1,$val2,$i1,$i2,@merged,@merlen,@dshlen,$bdn,$edn);

	my $c0 = -1;
	foreach $val0 (@input0) 
	{
		$c0++;
		# push @merged, "$val0 m8:1a";
		push @merged, "$val0 8mer";
		push @merlen, $len0[$c0];
		push @dshlen, $num0[$c0];
	}

	my $mer1 = 'F';
	my $c1 = -1;
	foreach $val1 (@input1) {
	$c1++;		
	if ($#input0 > -1) {
		foreach $num0 (@num0) {
			($bdn,$edn) = split /\-/, $num0;
			$mer1 = 'F';
			$i1 = $val1 + $edn +1;
					if (grep(/$i1 /,@merged)) {
				$mer1 = 'T';
				last, if ($mer1 eq 'T');
			}
		}
		}
		if ($mer1 eq 'F') {
			push @merged, "$val1 m8";
                        push @merlen, $len1[$c1];	
			push @dshlen, $num1[$c1];
		}
        }

	my $c2 = -1;
        foreach $val2 (@input2) {
		$c2++;
		my $mer2 = 'F';
		if ($#input0 > -1) {
               		( grep(/$val2 /,@merged))? ($mer2 = 'T'):();
		}
		if ($mer2 eq 'F') {
			push @merged, "$val2 1a";
			push @merlen, $len2[$c2];
			push @dshlen, $num2[$c2];
		}
		}
		
		if (@merged)
		{
			# print "arrcmp: ==> $utrID\t$mirFamID\t$speciesIDthisUTR\t@merged\t@merlen\t@dshlen <==\n";
		}


        return (\@merged,\@merlen,\@dshlen);
}

# utrcoords gets the UTR coordinates for each alignment match

sub utrcoords {
	my ($align,$hend,$type) = ($_[0],$_[1],$_[2]);	
	my $pseq = substr($align,0,$hend);
	my ($nstart,$nend);
	$pseq =~ s/\-//g;
	$nend = length($pseq); 
	
	$nstart = $nend;
	$nend += 6;

	if ($type eq 'm8:1a') 
	{
		# Correct 8mer OBOB: Nov 26, 2007
		$nend += 1;
	}
	
	return ($nstart,$nend); 
}

sub getGroupTypeSummaryInfo
{
	$groupSiteTypesList2GroupType{"1a"} = "1a";
	$groupSiteTypesList2GroupType{"m8"} = "m8";
	$groupSiteTypesList2GroupType{"8mer"} = "8mer";
	$groupSiteTypesList2GroupType{"1a 8mer"} = "8mer+1a";
	$groupSiteTypesList2GroupType{"8mer m8"} = "8mer+m8";
	$groupSiteTypesList2GroupType{"1a 8mer m8"} = "8mer+m8+1a";
	$groupSiteTypesList2GroupType{"1a m8"} = "m8+1a";
}
