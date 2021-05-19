package GraphParameters;

use strict;
use Carp;
use GD;
use Getopt::Long;

#----------------------------- Set default values for parameters --------------------------------#

# Width and height of the plot (the box excluding header and footer)
my $gWidth  = 800;
my $gHeight = $gWidth;

my $leftEdge   = 20;
my $rightEdge  = 20;
my $bottomEdge = 50;

my $graphTitleHt = 20;

my $majorTickLen = 6;
my $minorTickLen = 3;
my $xMinorTicks  = 5;
my $yMinorTicks  = 5;

my $mbFontWidth  = gdMediumBoldFont->width;
my $mbFontHeight = gdMediumBoldFont->height;
my $mlFontWidth  = gdLargeFont->width;

my $axisLabelGap = 20;  # Distance between tick labels and axis labels
my $tickGap = 5;        # Gap between the tick and tick label

my $showGd4  = 0;
my $dotSize  = 2;
my $rotx     = 0;
my $selRaces = "";
my $isBeyond = 0;
my $raceFile = "";
my $showCutoff = "";
my $showAreas  = "";

# Left edge of the plot (the rectangular box)
my $gxLeft  = $leftEdge + $mbFontHeight + $axisLabelGap + $mbFontWidth*2 + $tickGap + $majorTickLen;
my $gxRight = $gxLeft + $gWidth;
my $gyTop   = 10;
my $gyBtm   = $gyTop + $gHeight;

# Width and height of the whole graph
my $imageWidth  = $gxLeft + $gWidth + $rightEdge;
my $imageHeight = $gyTop + $gHeight + $graphTitleHt + $bottomEdge;

my $xCutMin = 0;
my $xCutMax = 0;
my $yCutMin = 0;
my $yCutMax = 0;

my $xMin = 1.0;
my $xMax = 1.8;
my $yMin = 1.0;
my $yMax = 1.8;

my $xRange = $xMax - $xMin;
my $yRange = $yMax - $yMin;

my $totNumAncSnps = 100437;
my $minSnps = 0;
my $maxSnps = $totNumAncSnps + 1;

# Set population cutoff values.
my @eurVtxCoord = (1.4785,  1.4460, 0.0000, 1);
my @afoVtxCoord = (1.0500,  1.1000, 0.0000, 1);
my @easVtxCoord = (1.7422,  1.1000, 0.0000, 1);
my @vertexCoords = (\@eurVtxCoord, \@afoVtxCoord, \@easVtxCoord);

my $eurCut      = 90;
my $afoCut      = 95;
my $easCut      = 95;
my $othLatCut   = 14;
my $afaLacCut   = 40;

my $meanSasx    = 1.524446;
my $meanSasy    = 0.079465;
my $sdSasx      = 0.019755;
my $sdSasy      = 0.005485;
my $asnLatCut   = 1.525;

my $numSasSds   = 4;
my $sasCutaVal  = 5;
my $sasCutBasey = $meanSasy - $numSasSds * $sdSasy;

my $asnCutBasex = 1.58;
my $meanAsny    = 0;
my $asnCutaVal  = 30;
my $sasLenCut   = 0.04;  # For separating South Asians from Hispanics using a different score

sub new
{
    my ($class) = @_;

    GetInputParameters();

    $rotx = $rotx % 360;

    # Reset y-axis limits when graph is rotated to leave room for points
    if (!$showGd4 && $rotx) {
        if ($rotx > 105 && $rotx < 255) {
            $yMin = 0.7;
            $yMax = 1.5;
        }
        elsif ($rotx > 75 && $rotx < 285) {
            $yMin = 0.9;
            $yMax = 1.7;
        }
    }

    # Do some (not all) sanity checks
    my $cutErr = "";
    if (defined $xMin && defined $xMax && $xMax - $xMin < 0.0999) {
        $cutErr .=  "\tDifference of x-max ($xMax) and x-min ($xMin) is less than 0.1.\n";
    }
    if (defined $xMin && defined $xMax && $xMax - $xMin > 1.5001) {
        $cutErr .=  "\tDifference of x-max ($xMax) and x-min ($xMin) is greater than 1.5.\n";
    }
    if (defined $yMin && defined $yMax && $yMax - $yMin < 0.0999) {
        $cutErr .=  "\tDifference of y-max ($yMax) and y-min ($yMin) is less than 0.1.\n";
    }
    if (defined $yMin && defined $yMax && $yMax - $yMin > 1.5001) {
        $cutErr .=  "\tDifference of y-max ($yMax) and y-min ($yMin) is greater than 1.5.\n";
    }

    $cutErr .= "\tminsnps ($minSnps) too large\n" if ($minSnps > $totNumAncSnps);
    $cutErr .= "\tminsnps ($minSnps) greater than maxsnps ($maxSnps)\n" if ($minSnps > $maxSnps);

    $cutErr .= "\tecut < 50\n"  if ($eurCut < 50);
    $cutErr .= "\tohcut < 5\n"  if ($othLatCut < 5);
    $cutErr .= "\tfhcut < 5\n"  if ($afaLacCut < 5);

    $cutErr .= "\tecut > 98\n"  if ($eurCut > 98);
    $cutErr .= "\tacut > 98\n"  if ($easCut > 98);
    $cutErr .= "\tfcut > 98\n"  if ($afoCut > 98);
    $cutErr .= "\tohcut > 50\n" if ($othLatCut > 50);
    $cutErr .= "\tfhcut > 80\n" if ($afaLacCut > 80);

    $cutErr .= "\tecut + ohcut < 95\n"  if ($eurCut + $othLatCut < 95);
    $cutErr .= "\tfcut + ohcut < 100\n" if ($afoCut + $othLatCut < 100 && $afoCut > 0);
    $cutErr .= "\tacut + ohcut < 100\n" if ($easCut + $othLatCut < 100 && $easCut > 0);
    $cutErr .= "\tfhcut > fcut\n"       if ($afaLacCut > $afoCut && $afoCut > 0);
    $cutErr .= "\tfhcut + ecut < 100\n" if ($afaLacCut < 100 - $eurCut);

    if ($cutErr) {
        print "\nERROR in population cutoff line settings:\n$cutErr\n";
        return "";
    }

    $showCutoff = 0 if ($rotx);

    # The cutoff line will not be drawn if they are set to be negative
    $afoCut = 200 if ($afoCut < 0);
    $easCut = 200 if ($easCut < 0);

    # When calculating populations, use cutoff values in (0, 1)
    $eurCut    /= 100.0;
    $afoCut    /= 100.0;
    $easCut    /= 100.0;
    $afaLacCut /= 100.0;
    $othLatCut /= 100.0;

    $minSnps = 0 unless ($minSnps);
    $maxSnps = 200000 unless ($maxSnps);

    if ($xCutMin && $xCutMax && $xCutMin > $xCutMax) {
        $xCutMin = 0;
        $xCutMax = 0;
    }

    if ($yCutMin && $yCutMax && $yCutMin > $yCutMax) {
        $yCutMin = 0;
        $yCutMax = 0;
    }

    # When the user doesn't select races to show, nor which regions to get subjects, show all subjects
    my $isSelAll = 1;
    $isSelAll = 0 if ($selRaces || $xCutMin || $xCutMax || $yCutMin || $yCutMax);

    if ($dotSize && $dotSize =~ /^\d+$/) {
        $dotSize = 1 if ($dotSize < 1);
        $dotSize = 10 if ($dotSize > 10);
    }
    else {
        $dotSize = 2;
    }

    my @alfaPops = ("EUR", "AFO", "EAS", "AFA", "LAC", "LEN", "OAS", "SAS");
    my @alfaFullPops = ("European", "African", "East Asian", "African American",
                        "Latin American 1", "Latin American 2", "Asian-Pacific Islander", "South Asian");

    bless {
        leftEdge     => $leftEdge,
        rightEdge    => $rightEdge,
        gWidth       => $gWidth,
        gHeight      => $gHeight,
        gxLeft       => $gxLeft,
        gxRight      => $gxRight,
        gyTop        => $gyTop,
        gyBtm        => $gyBtm,

        tickGap      => $tickGap,
        axisLabelGap => $axisLabelGap,
        majorTickLen => $majorTickLen,
        minorTickLen => $minorTickLen,
        xMinorTicks  => $xMinorTicks,
        yMinorTicks  => $yMinorTicks,
        dotSize      => $dotSize,

        mlFontWidth  => $mlFontWidth,
        mbFontWidth  => $mbFontWidth,
        mbFontHeight => $mbFontHeight,
        imageWidth   => $imageWidth,
        imageHeight  => $imageHeight,

        minSnps      => $minSnps,
        maxSnps      => $maxSnps,

        xMin         => $xMin,
        xMax         => $xMax,
        yMin         => $yMin,
        yMax         => $yMax,
        xRange       => $xRange,
        yRange       => $yRange,

        xCutMin      => $xCutMin,
        xCutMax      => $xCutMax,
        yCutMin      => $yCutMin,
        yCutMax      => $yCutMax,

        eurCut       => $eurCut,
        afoCut       => $afoCut,
        easCut       => $easCut,
        othLatCut    => $othLatCut,
        afaLacCut    => $afaLacCut,

        showGd4      => $showGd4,

        raceFile     => $raceFile,
        rotx         => $rotx,
        selRaces     => $selRaces,
        isBeyond     => $isBeyond,
        showCutoff   => $showCutoff,
        showAreas    => $showAreas,

        meanSasx     => $meanSasx,
        meanSasy     => $meanSasy,
        sdSasx       => $sdSasx,
        sdSasy       => $sdSasy,
        sasCutBasey  => $sasCutBasey,
        sasCutaVal   => $sasCutaVal,
        sasLenCut    => $sasLenCut,
        asnLatCut    => $asnLatCut,

        asnCutBasex  => $asnCutBasex,
        meanAsny     => $meanAsny,
        asnCutaVal   => $asnCutaVal,

        alfaPops     => \@alfaPops,
        alfaFullPops => \@alfaFullPops,

    	vtxCoords    => \@vertexCoords,
    }, $class;
}

#
# Read parameters from command line and overwrite the defaulted ones
#
sub GetInputParameters
{
    my ($inFile, $outFile) = ("", "");
    my %inputParams = ();

    if (@ARGV > 1) {
        $inFile = $ARGV[0];
        $outFile = $ARGV[1];

        for my $i (2 .. $#ARGV) {
            my $arg = $ARGV[$i];
            if ($arg =~ /^-([A-Za-z]\w*)$/) {
                my $param = $1;
                if ($i == $#ARGV) {
                    $inputParams{$param} = 1;
                }
                else {
                    my $val = $ARGV[$i+1];
                    if ($val =~ /^-[A-Za-z]\w*$/) {
                        $inputParams{$param} = 1;
                    }
                    else {
                        $inputParams{$param} = $val;
                        $i++;
                    }
                }
            }
        }
    }

    if ($inFile =~ /^-/ || $outFile =~ /^-/) {
    	($inFile, $outFile) = ("", "");
    }

    foreach my $param (keys %inputParams) {
        SetInputParameter($param, $inputParams{$param});
    }

    if ($showGd4) {
        $yMin = -0.3 unless (defined $inputParams{"ymin"});
        $yMax =  0.5 unless (defined $inputParams{"ymax"});
    }

    if (defined $inputParams{"ecut"} ||
        defined $inputParams{"fcut"} ||
        defined $inputParams{"acut"} ||
        defined $inputParams{"ohcut"} ||
        defined $inputParams{"fhcut"} ) {
        $showCutoff = 1
    }
    return ($inFile, $outFile);
}

sub SetInputParameter
{
    my ($param, $value) = @_;

    $gWidth       =  $value if ($param eq "gw");
    $dotSize      =  $value if ($param eq "dot");

    $minSnps      =  $value if ($param =~ /^minsnp/);
    $maxSnps      =  $value if ($param =~ /^maxsnp/);

    $xMin         =  $value if ($param eq "xmin");
    $xMax         =  $value if ($param eq "xmax");
    $yMin         =  $value if ($param eq "ymin");
    $yMax         =  $value if ($param eq "ymax");

    $xCutMin      =  $value if ($param eq "xcmin");
    $xCutMax      =  $value if ($param eq "xcmax");
    $yCutMin      =  $value if ($param eq "ycmin");
    $yCutMax      =  $value if ($param eq "ycmax");

    $eurCut       =  $value if ($param eq "ecut");
    $afoCut       =  $value if ($param eq "fcut");
    $easCut       =  $value if ($param eq "acut");
    $othLatCut    =  $value if ($param eq "ohcut");
    $afaLacCut    =  $value if ($param eq "fhcut");

    $showGd4      =  $value if ($param eq "gd4");
    $raceFile     =  $value if ($param eq "spf");
    $rotx         =  $value if ($param eq "rotx");
    $selRaces     =  $value if ($param eq "pops");
    $isBeyond     =  $value if ($param eq "isByd");
    $showCutoff   =  $value if ($param eq "cutoff");
    $showAreas    =  $value if ($param eq "areas");

    $gWidth  =  500 if ($gWidth <  500);
    $gWidth  = 2000 if ($gWidth > 2000);
    $gHeight = $gWidth;
    $gxLeft  = $leftEdge + $mbFontHeight + $axisLabelGap + $mbFontWidth*2 + $tickGap + $majorTickLen;
    $gxRight = $gxLeft + $gWidth;
    $gyBtm   = $gyTop + $gHeight;

    $imageWidth  = $gxLeft + $gWidth + $rightEdge;
    $imageHeight = $gyTop + $gHeight + $graphTitleHt + $bottomEdge;

    $xRange  = $xMax - $xMin;
    $yRange  = $yMax - $yMin;

    # Round min, max values to the hundreds
    $xMin = int($xMin * 100 + 0.5) / 100;
    $xMax = int($xMax * 100 + 0.5) / 100;
    $yMax = int($yMax * 100 + 0.5) / 100;
    if ($yMin < 0) {
        $yMin = int($yMin * 100 - 0.5) / 100;
    }
    else {
        $yMin = int($yMin * 100 + 0.5) / 100;
    }
}

1;
