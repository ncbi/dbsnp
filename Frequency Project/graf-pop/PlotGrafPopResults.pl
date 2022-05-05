#!/usr/local/bin/perl

my $disclaim = << "EOF";
    ====================================================================================
                               PUBLIC DOMAIN NOTICE
                National Center for Biotechnology Information

        This software/database is a "United States Government Work" under the
        terms of the United States Copyright Act.  It was written as part of
        the author's official duties as a United States Government employee and
        thus cannot be copyrighted.  This software/database is freely available
        to the public for use. The National Library of Medicine and the U.S.
        Government have not placed any restriction on its use or reproduction.
        Although all reasonable efforts have been taken to ensure the accuracy
        and reliability of the software and data, the NLM and the U.S.
        Government do not and cannot warrant the performance or results that
        may be obtained by using this software or data. The NLM and the U.S.
        Government disclaim all warranties, express or implied, including
        warranties of performance, merchantability or fitness for any particular
        purpose.

        Please cite the author in any work or product based on this material.

        Author: Yumi (Jimmy) Jin (jinyu\@ncbi.nlm.nih.gov)
        File Description: script to plot graphs to show subject populations.
        Date: 03/01/2021
    ====================================================================================
EOF

BEGIN {
    use Cwd 'abs_path';
    my ($scriptName, $scriptDir) = ("", "");
    my $scriptFullname = abs_path($0);
    ($scriptDir, $scriptName) = ($1, $2) if ($scriptFullname =~ /(\S+)\/(\S+)/);
    push ( @INC, $scriptDir);
}

use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use strict;
use Time::HiRes qw(gettimeofday);
use Data::Dumper;
use GD;
use GD::Text;
use GD::Graph;
use GD::Graph::colour;
use GD::Graph::lines;

use GraphColors;
use GrafPopFiles;
use GraphParameters;
use PopulationCutoffs;
use GraphTransformation;
use SubjectAncestry;

if (@ARGV < 2) {
    my $usage = GetScriptUsage();
    print "$disclaim\n$usage\n\n";
    exit;
}

my $inFile = $ARGV[0];
my $outFile = $ARGV[1];
my $showSbjs = 0;
$showSbjs = 1 if ($outFile =~ /\.txt$/);

my $param = new GraphParameters();
exit unless ($param);

my $cgi = new CGI;
my $img = new GD::Image($param->{imageWidth}, $param->{imageHeight});
my $colors = new GraphColors($img);

#------------------------- Save some frequently used parameters ----------------------#
my $xMin  = $param->{xMin};
my $xMax  = $param->{xMax};
my $yMin  = $param->{yMin};
my $yMax  = $param->{yMax};
my $xRange  = $param->{xRange};
my $yRange  = $param->{yRange};

my $gxLeft  = $param->{gxLeft};
my $gxRight  = $param->{gxRight};
my $gyTop = $param->{gyTop};
my $gyBottom = $param->{gyBtm};
my $graphWidth = $param->{gWidth};
my $graphHeight = $param->{gHeight};
my $mbFontWidth = $param->{mbFontWidth};
my $mbFontHeight = $param->{mbFontHeight};
my $black = $colors->{black};

#---------------- Global variables, avoid passing to save runtime -------------------#
my @sbjPopScores =  (); # Ancestry scores saved in the input file
my @allSbjPvalues = (); # Subject GD scores (x, y, z, 1) for plotting

#--------------------------- Read subject GrafPop scores  ---------------------------#
my ($popScoreRef, $allPopSbjs, $minSbjSnps, $maxSbjSnps, $meanSbjSnps, $error)
    = GrafPopFiles::ReadGrafPopResults($inFile, $param->{minSnps}, $param->{maxSnps});
if ($error) {
    print "\nERROR: $error\n";
    exit;
}

@sbjPopScores = @$popScoreRef;
my $numSbjs = @sbjPopScores;

#--------------------------- Read subject races from file ---------------------------#
my %sbjRaces = ();
my %allRaces = ();
my $hasRaceInfo = 0;
my $numRaces = 0;
if ($param->{raceFile}) {
    my ($sbjRaceRef, $allRaceRef, $hasRace, $err) = GrafPopFiles::ReadSubjectRaces($param->{raceFile}, $allPopSbjs);
    %sbjRaces = %$sbjRaceRef;
    %allRaces = %$allRaceRef;
    $hasRaceInfo = $hasRace;
    if ($err) {
        print "\nERROR: $err\n";
        exit;
    }
}

#---------------- Get subject counts for different self-reported races --------------#
my %raceSbjCnts = ();
my $unkRace = "NOT REPORTED";
foreach my $sbjNo (0 .. $#sbjPopScores) {
    my $sbj = $sbjPopScores[$sbjNo]->{subject};

    my $race = $unkRace;
    $race = $sbjRaces{$sbj} if ($sbjRaces{$sbj});

    if ($raceSbjCnts{$race}) {
    	$raceSbjCnts{$race}++;
    }
    else {
	    $raceSbjCnts{$race} = 1;
    }
}

# Sort self-reported races by subject counts
my @sortRaces = ();
foreach my $race (sort {$raceSbjCnts{$b} <=> $raceSbjCnts{$a}} keys %raceSbjCnts) {
    my $cnt = $raceSbjCnts{$race};
    push @sortRaces, $race;
}
$numRaces = @sortRaces;

#----------------------- Process input list of races to display -----------------------#
my ($raceNameRef, $raceIdRef, $raceColorRef) = GetDisplayRaces();
my %selRaceNames = %$raceNameRef;
my %raceIdNos = %$raceIdRef;
my %raceColorIds = %$raceColorRef;

#------------------------ Set race and color for each subject -------------------------#
foreach my $sbjNo (0 .. $#sbjPopScores) {
    my %info = %{$sbjPopScores[$sbjNo]};
    my $sbj = $info{subject};
    my $race = $unkRace;
    $race = $sbjRaces{$sbj} if ($sbjRaces{$sbj});
    my $raceId = $selRaceNames{$race};
    my $raceNo = $raceIdNos{$raceId};
    my $colorNo = $raceColorIds{$race};
    my $color = $colors->{raceColors}->[$colorNo];

    my $gd1 = $info{x};
    my $gd2 = $info{y};
    my $gd3 = $info{z};
    my $gd4 = $info{gd4};

    my @pvals = ($gd1, $gd2, $gd3, 1);
    if ($param->{showGd4}) {
	    @pvals = ($gd1, $gd4, $gd3, 1);
    }
    push @allSbjPvalues, \@pvals;

    $sbjPopScores[$sbjNo]->{race} = $race;
    $sbjPopScores[$sbjNo]->{raceNo} = $raceNo;
    $sbjPopScores[$sbjNo]->{color} = $color;
    $sbjPopScores[$sbjNo]->{colorNo} = $colorNo;
}

print "Total $numSbjs samples have ancestry scores\n";

my $cutoffValues = new PopulationCutoffs($img, $param);

if ($param->{rotx}) {
    RotateSubjectsOnx($param->{rotx});
}

#---------------------------------- Plot graph ---------------------------------------#
PlotPopulations($outFile);


#---------------------------------- Subroutines ---------------------------------------#

#
# Plot graph to show subject ancestry
#
sub PlotPopulations
{
    my $outPngFile = shift;

    open(IMG, ">$outPngFile") or die $!;
    binmode IMG;

    # Plot axes, y-axis can be GD2 or GD4
    my $xLabel = "GD1";
    my $yLabel = "GD2";

    if ($param->{rotx}) {
	    $yLabel = "GD2, GD3: rotated by $param->{rotx} degree";
    }
    if ($param->{showGd4}) {
	    $yLabel = "GD4";
    }
    PlotAxes($gyTop, $xLabel, $yLabel);

    # Legends to show what color represents which self-reported ancestry
    my $lgdx = $gxLeft + 10;
    my $lgdy = $gyTop + 10;
    my $lgdSide = 10;
    my $lgdGap = $lgdSide + 10;

    my $snpLgdGap = 105; # Room to show "SNPs/Sbj"
    $snpLgdGap = 85 if ($graphWidth < 401);

    my $hasOther = 0;
    my $numSelSbjs = 0;
    my $numOthSbjs = 0;
    my $raceNo = 1;
    print "\nSelf-reported races/ethnicities\n" if (@sortRaces > 1);
    foreach my $raceId (0 .. $#sortRaces) {
	    my $raceNum = $raceId + 1;
	    my $race = $sortRaces[$raceId];
	    my $cnt = $raceSbjCnts{$race};
        print "$raceNo: $race (n=$cnt)\n" if (@sortRaces > 1);
	    my $dispRace = TruncateDisplayRace($race, $cnt, $lgdSide, $snpLgdGap);
	    my $colorNo = $raceColorIds{$race};
	    my $color = $colors->{raceColors}->[$colorNo];
	    $color = $black unless ($hasRaceInfo);
	    if ($colorNo > 0) {
	        $img->filledRectangle($lgdx, $lgdy, $lgdx + $lgdSide, $lgdy + $lgdSide, $color);
	        $img->string(gdMediumBoldFont, $lgdx + $lgdSide + 10, $lgdy, "$dispRace ($cnt)", $black);
	        $lgdy += $lgdGap;
	        $numSelSbjs += $cnt;
	    }
	    else {
	        $hasOther = 1;
	        $numOthSbjs += $cnt;
	    }

	    $raceNo++;
    }

    if ($hasOther) {
	    my $color = $colors->{raceColors}->[0];
	    $img->filledRectangle($lgdx, $lgdy, $lgdx + $lgdSide, $lgdy + $lgdSide, $color);
	    $img->string(gdMediumBoldFont, $lgdx + $lgdSide + 10, $lgdy, "unselected ($numOthSbjs)", $black);
	    $lgdy += $lgdGap;
    }

    $img->string(gdMediumBoldFont, $lgdx, $lgdy, "Total $numSbjs subjects", $black);

    # Show min and max Anc SNPs with genotypes per subject
    my $snpLgdx = $gxLeft + $graphWidth - $snpLgdGap;
    my $snpLgdy = $gyTop + 10;
    my $snpLgdGap = $mbFontHeight * 1.5;

    $img->string(gdMediumBoldFont, $snpLgdx, $snpLgdy, " SNPs/Sbj:", $black);

    $snpLgdx += $mbFontWidth;
    $snpLgdy += $snpLgdGap;
    my $minSnpText = sprintf("Min  %6d", $minSbjSnps);
    $img->string(gdMediumBoldFont, $snpLgdx, $snpLgdy, $minSnpText, $black);

    $snpLgdy += $snpLgdGap;
    my $maxSnpText = sprintf("Max  %6d", $maxSbjSnps);
    $img->string(gdMediumBoldFont, $snpLgdx, $snpLgdy, $maxSnpText, $black);

    $snpLgdy += $snpLgdGap;
    my $meanSnpText = sprintf("Mean %6.0f", $meanSbjSnps);
    $img->string(gdMediumBoldFont, $snpLgdx, $snpLgdy, $meanSnpText, $black);

    # Plot subjects, one dot for one subject
    PlotSubjects($hasRaceInfo);

    # Plot vertices, cutoff lines
    $cutoffValues->PlotVertices($colors->{black});
    if ($param->{showCutoff}) {
        if ($param->{showGd4}) {
            $cutoffValues->PlotSasCutoffLines($colors->{gold});
        }
        else {
            $cutoffValues->PlotCutoffLines($colors->{gold})
        }
    }

    # Show expected areas for populations specified by the user
    if (!$param->{showGd4} && !$param->{rotx}) {
        my @predAreaPops = split /\,/, $param->{showAreas};
        for my $popId (@predAreaPops) {
            if ($popId > 0 && $popId < 9) {
                $cutoffValues->ShowExpectedArea($popId-1, $colors->{expAreaColor});
            }
        }
    }

    print IMG $img->png;
    print "Graph saved to $outPngFile\n\n";

    exit;
}

#
# Plot all subjects, one dot for each subject
#
sub PlotSubjects
{
    my ($hasRace) = @_;

    # Find the top and bottom positions to show population labels without blocking data points
    my $eurMaxy = 0;
    my $afoMiny = 10;
    my $easMiny = 10;

    my $eurLblLen = length(" European ") * $mbFontWidth;
    my $afoLblLen = length(" Africann ") * $mbFontWidth;
    my $easLblLen = length(" East Asian ") * $mbFontWidth;

    my $eurVx = $param->{vtxCoords}->[0]->[0];
    my $afoVx = $param->{vtxCoords}->[1]->[0];
    my $easVx = $param->{vtxCoords}->[2]->[0];

    # Z-values determine which subjects are in the front
    my %sbjZvals = ();
    foreach my $sbjNo (0 .. $#allSbjPvalues) {
	    my ($xVal, $yVal, $zVal, $junk) = @{$allSbjPvalues[$sbjNo]};
	    $sbjZvals{$sbjNo} = $zVal;

	    $eurMaxy = $yVal if ($xVal > $eurVx - $eurLblLen/2 && $xVal < $eurVx + $eurLblLen/2 && $yVal > $eurMaxy);
	    $afoMiny = $yVal if ($xVal > $afoVx - $afoLblLen/2 && $xVal < $afoVx + $afoLblLen/2 && $yVal < $afoMiny);
	    $easMiny = $yVal if ($xVal > $easVx - $easLblLen/2 && $xVal < $easVx + $afoLblLen/2 && $yVal < $easMiny);
    }

    my @sortSbjNos = ();
    foreach my $sbjNo (sort {$sbjZvals{$b} <=> $sbjZvals{$a}} keys %sbjZvals) {
	    push @sortSbjNos, $sbjNo;
    }

    # Plot the subjects of non-selected ancestries to the back
    my $halfDot = $param->{dotSize}/2;
    foreach my $sbjNo (0 .. $#sortSbjNos) {
	    unless ($sbjPopScores[$sbjNo]->{colorNo}) {
	        my ($xVal, $yVal, $zVal, $junk) = @{$allSbjPvalues[$sbjNo]};

	        if ($xVal > $xMin && $xVal < $xMax && $yVal > $yMin && $yVal < $yMax) {
		        my $color = $sbjPopScores[$sbjNo]->{color};
		        $color = $black unless ($hasRace);
                PlotOneDot($xVal, $yVal, $param->{dotSize}, $color);
	        }
	    }
    }

    # Plot the subjects of selected ancestries to the front
    foreach my $sbjNo (0 .. @sortSbjNos) {
	    if ($sbjPopScores[$sbjNo]->{colorNo}) {
	        my ($xVal, $yVal, $zVal, $junk) = @{$allSbjPvalues[$sbjNo]};

	        if ($xVal > $xMin && $xVal < $xMax && $yVal > $yMin && $yVal < $yMax) {
		        my $color = $sbjPopScores[$sbjNo]->{color};
		        $color = $black unless ($hasRace);
		        PlotOneDot($xVal, $yVal, $param->{dotSize}, $color);
	        }
	    }
    }

    # Show vertex labels
    if (!$param->{showGd4} && !$param->{rotx}) {
        my ($eurLbl, $afoLbl, $easLbl) = ("European", "African", "East Asian");
        my $vetGap = 0.02;
        my ($eurLbly, $afoLbly, $easLbly) = (1.5, 1.05, 1.05);
        $eurLbly = $eurMaxy + $vetGap if ($eurLbly < $eurMaxy + $vetGap);
        $afoLbly = $afoMiny - $vetGap if ($afoLbly > $afoMiny - $vetGap);
        $easLbly = $easMiny - $vetGap if ($easLbly > $easMiny - $vetGap);

        PlotVertexLabel($eurLbl, 0, $eurLbly);
        PlotVertexLabel($afoLbl, 1, $afoLbly);
        PlotVertexLabel($easLbl, 2, $easLbly);
    }
}

#
# Plot the x and y axes, and the tick marks
#
sub PlotAxes
{
    my ($yTop, $xLabelStr, $yLabelStr) = @_;

    my $yBottom = $yTop + $graphHeight;
    my $format = "%4.2f";
    my $tickLblLen = 4;

    # Plot y-axis
    $img->line($gxLeft, $yTop, $gxLeft, $yBottom, $black);
    $img->line($gxRight, $yTop, $gxRight, $yBottom, $black);

    # Plot ticks on y-axis
    my ($yMinorIntv, $yNumMinors, $yTickNoFirstMajor) = GetTickIntervals($param->{yMin}, $param->{yMax});
    my $yTotMinors = int($param->{yRange} / $yMinorIntv + 0.05);
    my $dy = $graphWidth/$yTotMinors;

    my $xMinor = $gxLeft - $param->{minorTickLen};
    my $xMajor = $gxLeft - $param->{majorTickLen};
    my $xString = $gxLeft - $mbFontWidth*$tickLblLen - $param->{majorTickLen} - $param->{tickGap};

    for my $i (0 .. $yTotMinors) {
	    my $val = $param->{yMin} + $yMinorIntv * $i;
	    my $valStr = sprintf($format, $val);
	    my $y = int($yBottom - $dy*$i + 0.5);
        $y = $yBottom if ($y > $yBottom);
        if ($y >= $yTop) {
            $img->line($gxLeft, $y, $xMinor, $y, $black);
            if (($i + 1 - $yTickNoFirstMajor) % $yNumMinors == 0) {
                $img->line($xMinor, $y, $xMajor, $y, $black);
                $img->string(gdMediumBoldFont, $xString, $y-$mbFontHeight/2, $valStr, $black);
             }
        }
    }

    # Draw y-axis label
    my $lblLen = length($yLabelStr);
    my $yLabel = $yBottom - ($graphHeight - $mbFontWidth*$lblLen)/2;
    $img->stringUp(gdMediumBoldFont, $xString-25, $yLabel, $yLabelStr, $black);

    # Plot x-axis
    $img->line($gxLeft, $yBottom, $gxRight, $yBottom, $black);
    $img->line($gxLeft, $yTop, $gxRight, $yTop, $black);

    # Plot x-axis ticks
    my ($xMinorIntv, $xNumMinors, $xTickNoFirstMajor) = GetTickIntervals($param->{xMin}, $param->{xMax});
    my $xTotMinors = int($param->{xRange} / $xMinorIntv + 0.05);
    my $dx = $graphWidth/$xTotMinors;

    my $yMinor = $yBottom + $param->{minorTickLen};
    my $yMajor = $yBottom + $param->{majorTickLen};
    my $yString = $yMajor + 5;

    for my $i (0 .. $xTotMinors) {
	    my $val = $param->{xMin} + $xMinorIntv * $i;
	    my $valStr = sprintf($format, $val);
	    my $x = int($gxLeft + $dx*$i + 0.5);
	    $img->line($x, $yBottom, $x, $yMinor, $black);
        if (($i + 1 - $xTickNoFirstMajor) % $xNumMinors == 0) {
    	    $img->line($x, $yMinor, $x, $yMajor, $black);
    	    $img->string(gdMediumBoldFont, $x-$mbFontWidth*2, $yString, $valStr, $black);
        }
    }

    # Draw x-axis label
    my $lblLen = length($xLabelStr);
    my $xLabel = $gxLeft + ($graphWidth - $mbFontWidth*$lblLen)/2;
    $img->string(gdMediumBoldFont, $xLabel, $yMajor + 25, $xLabelStr, $black);
}

#
# Set tick intervals based on value range
#
sub GetTickIntervals
{
    my ($minVal, $maxVal) = @_;

    my $range = $maxVal - $minVal;
    my $minorIntv = 0.02; # value diff between two adjacent minor ticks
    my $numMinors = 5;    # number of minor ticks in each major tick range

    if ($range < 0.5) {
        $minorIntv = 0.01;
    }
    elsif ($range < 0.2) {
        $minorIntv = 0.005;
        $numMinors = 2;
    }

    my $majorIntv = $minorIntv * $numMinors;
    my $startVal = int(($minVal+0.0001)/$majorIntv) * $majorIntv;

    # Which minor tick is the first major tick
    my $tickNoFirstMajor = int(($startVal + $majorIntv - $minVal + 0.0001) / $minorIntv) + 1;

    return ($minorIntv, $numMinors, $tickNoFirstMajor);
}

#
# Plot population labels for the three vertices of the triangle
#
sub PlotVertexLabel
{
    my ($vLbl, $popNo, $vLbly, $extraGap) = @_;

    my $lblLen = length($vLbl) * $mbFontWidth;
    my ($xVal, $yVal, $zVal, $junk) = @{$param->{vtxCoords}->[$popNo]};
    my $xPos = $gxLeft  + ($xVal - $xMin) * $graphWidth  / ($xMax - $xMin) - $lblLen/2.0;
    my $yPos = $gyBottom - ($vLbly - $yMin) * $graphHeight / ($yMax - $yMin) + $extraGap;
    $yPos -=  $mbFontHeight / 2; # The position of the middle of the label
    if ($xPos > $gxLeft + 10 && $xPos + $lblLen < $gxRight - 10 &&
        $yPos > $gyTop + 10 && $yPos < $gyBottom - 20) {
        $img->string(gdMediumBoldFont, $xPos, $yPos, $vLbl, $black);
    }
}

#
# Check if there is enough room to show the self-reported ancestry.
# If not, truncate the ancestry.
#
sub TruncateDisplayRace
{
    my ($race, $cnt, $lgdSide, $snpLgdGap) = @_;

    my $dispRace = $race;

    my $raceLen = length($race);
    my $cntLen = length($cnt);
    my $lblLen = $raceLen + $cntLen + 4;

    my $extraLen = 10 + $lgdSide + $lblLen * $mbFontWidth + $snpLgdGap - $graphWidth;
    if ($extraLen > 0) {
	    my $cutChars = int($extraLen * 1.0 / $mbFontWidth) + 1;
	    $dispRace = substr($race, $raceLen);
	    $dispRace = substr($race, 0, $raceLen - $cutChars - 3);
	    $dispRace .= "...";
    }

    return $dispRace;
}

#
# Get races selected by user to display on the graph
#
sub GetDisplayRaces
{
    my $selRaces = $param->{selRaces};
    my @origRaceVals = split /\,/, $selRaces; # The race IDs selected to display on the graph
    my @selRaceVals = ();
    for my $raceVal (@origRaceVals) {
        push @selRaceVals, $raceVal if ($raceVal =~ /^\d+$/ && $raceVal > 0 && $raceVal <= $numRaces);
    }

    if (!$selRaces) {                        # show top 10 races by default
        for my $raceNo (0 .. $#sortRaces) {
            my $raceId = $raceNo + 1;
            my $race = $sortRaces[$raceNo];
            push @selRaceVals, $raceId;
            $selRaces .= "$raceId,";
            last if ($raceId >= $colors->{maxShowRaces});
        }
        $selRaces =~ s/\,\s*$//;
    }

    # Create a hash to find race name from the selected race ID.
    my %selRaceNames = ();
    my %selRaceIds = ();
    for my $val (@selRaceVals) {
        my $id = $val - 1;
        my $race = $sortRaces[$id];
        $selRaceIds{$id} = 1;
        $selRaceNames{$race} = $id;
    }
    my $numSelRaces = keys %selRaceIds;

    # Create a hash to find color for a self-reported ancestry
    my %raceColorIds = ();
    my $selRaceNo = 1;
    foreach my $raceId (0 .. $#sortRaces) {
        my $race = $sortRaces[$raceId];

        my $colorNo = 0;
        if (defined $selRaceIds{$raceId}) {
            $colorNo = $selRaceNo;
            $selRaceNo++;
        }
        $colorNo = $colorNo % $colors->{maxShowRaces};
        $raceColorIds{$race} = $colorNo;
    }

    my $raceNo = 0;
    my %raceIdNos = ();
    foreach my $raceId (sort {$a <=> $b} keys %selRaceIds) {
        $raceIdNos{$raceId} = $raceNo;
        $raceNo++;
    }

    return (\%selRaceNames, \%raceIdNos, \%raceColorIds);
}


#
# Plot a filled circle give x, y scores
#
sub PlotOneDot
{
    my ($xVal, $yVal, $size, $color) = @_;

    my $x = $gxLeft   + ($xVal - $xMin) * $graphWidth  * 1.0 / $xRange;
    my $y = $gyBottom - ($yVal - $yMin) * $graphHeight * 1.0 / $yRange;
    $x = int($x + 0.5);
    $y = int($y + 0.5);

    return if ($x < $gxLeft+$size || $x > $gxRight-$size || $y < $gyTop+$size || $y > $gyBottom-$size);

    # Acoid using fill() since:
    # 1. It is hard to find spots to fill when circles overlap
    # 2. It may overspill, especially when the dots are close to the edges
    $img->setPixel($x, $y, $color);

    if ($size > 1) {
        $img->arc($x, $y, $size, $size, 0, 360, $color);
        if ($size > 4) {
            $img->setPixel($x+1, $y, $color);
            $img->setPixel($x-1, $y, $color);
            $img->setPixel($x, $y+1, $color);
            $img->setPixel($x, $y-1, $color);
        }
        if ($size > 5) {
            for my $iSize (5 .. $size-1) {
                $img->arc($x, $y, $iSize, $iSize, 0, 360, $color);
            }
        }
    }
}

#
# Rotate all subjects around the x-axis by an angle
#
sub RotateSubjectsOnx
{
    my $theta = shift;

    my @afoVtxCoords = @{$param->{vtxCoords}->[1]};
    my ($afoVx, $afoVy, $afoVz) = @afoVtxCoords;

    GraphTransformation::MoveShape3D(-$afoVx, -$afoVy, -$afoVz, \@allSbjPvalues);
    GraphTransformation::RotateShape3D(-$theta, "x", \@allSbjPvalues);
    GraphTransformation::MoveShape3D($afoVx, $afoVy, $afoVz, \@allSbjPvalues);
}

sub GetScriptUsage
{
    my $usage = "Usage: PlotGrafPopResults.pl <input file> <output file> [Options]

    Note:
          Input file is the file generated by the C++ grafpop program that includes subject ancestry scores.
          Output file should be a .png file.
          Options should be entered after the two required parameters.

    Options:
        Specify the input file with self-reported subject race information
            -spf     text file with two columns (no header): subject and self-reported population

        Set window size in pixels
            -gw      graph width (500 - 2000, default 800)

        Set graph axis limits, max - min should be between 0.1 and 1.5
            -xmin    min x value
            -xmax    max x value
            -ymin    min y value
            -ymax    max y value

        Set minimum and maximum numbers of genotyped fingerprint SNPs for samples to be processed
            -minsnp  minimum number of SNPs with genotypes
            -maxsnp  maximum number of SNPs with genotypes

        Set population cutoff lines
            -ecut    proportion: cutoff European proportion dividing Europeans from other populations. Default 90%.
            -fcut    proportion: cutoff African proportion dividing Africans from other populations. Default 95%.
            -acut    proportion: cutoff East Asian proportion dividing East Asians from other populations. Default 95%.
            -ohcut   proportion: cutoff African proportion dividing Latin Americans from Other population. Default 14%.
            -fhcut   proportion: cutoff African proportion dividing Latin Americans from African Americans. Default 40%.

        Select some self-reported populations (by IDs) to be highlighted on the graph (for studies with multiple races)
            -pops    comma separated population IDs, e.g., -pops 1,3,4 -> highlight populations #1, #3 and #4

        Select self-reported populations (by IDs) to show areas including 95% dbGaP subjects with genotypes of
        at least 50,000 ancestry SNPs, to help estimate subject populations
            -areas   comma separated dbGaP self-reported population IDs, e.g., -areas 1,3
                         -> show areas that include 95% dbGaP subjects with self-reported populations #1 and #3
                            1: European                 2: African
                            3: East Asian               4: African American
                            5: Latin American 1         6: Latin American 2
                            7: Asian-Pacific Islander   8: South Asian

        Select which score to show on the y-axis
            -gd4:    show GD4 on y-axis (GD4 separates South Asians from Latin Americans and other Asians)

        Set population cutoff lines
            -cutoff: show cutoff lines

        Rotate the plot on x-axis by a certain angle
            -rotx    angle (in degrees) to rotate the GD2/GD3 vs. GD1 plot on the x-axis (0 - 360)

        Set the size (diameter) of each dot that represents each subject
            -dot     dot size in pixels (1 - 10, default 2)";

    return $usage;
}