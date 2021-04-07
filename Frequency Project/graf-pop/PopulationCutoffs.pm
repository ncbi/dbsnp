package PopulationCutoffs;

use strict;
use Carp;
use GD;
use GD::Text;
use GD::Graph;
use GD::Graph::colour;
use GD::Graph::lines;
use GraphColors;
use GraphTransformation;

#------------------- Set some global variables for convenience ----------------------------#

my $pi = 3.1415936535;

my ($eurVx, $eurVy) = (0, 0);
my ($afoVx, $afoVy) = (0, 0);
my ($easVx, $easVy) = (0, 0);

my $cutLineLen = 0.065;  # Length of the cutoff lines outside the triangle

sub new
{
    my ($class, $img, $param) = @_;

    my @eurVtxCoord = @{$param->{vtxCoords}->[0]};
    my @afoVtxCoord = @{$param->{vtxCoords}->[1]};
    my @easVtxCoord = @{$param->{vtxCoords}->[2]};

    ($eurVx, $eurVy) = ($eurVtxCoord[0], $eurVtxCoord[1]);
    ($afoVx, $afoVy) = ($afoVtxCoord[0], $afoVtxCoord[1]);
    ($easVx, $easVy) = ($easVtxCoord[0], $easVtxCoord[1]);

    bless {
        img => $img,
        param => $param
    }, $class;
}

#
# Plot the oval expected area for a given population
#
sub ShowExpectedArea
{
    my ($self, $popNo, $color) = @_;

    my @allPopParams = GetPopulationParameters();
    my %popParams = $popNo < 8 ? %{$allPopParams[$popNo]} : ();
    my $sdRatio = $popParams{sdy} / $popParams{sdx};

    my ($prevx, $prevy) = (0, 0);
    for my $i (0 .. 360) {
        my $ang = $i * $pi / 180;
        my $xVal1 = $popParams{meanx} + $popParams{radius} * cos($ang);
        my $yVal1 = $popParams{meany} + $popParams{radius} * sin($ang) * $sdRatio;
        my ($xVal, $yVal) = GraphTransformation::TransformPoint2D($xVal1, $yVal1, -$afoVx, -$afoVy, $popParams{theta});

        my $x = $self->{param}->{gxLeft}+($xVal-$self->{param}->{xMin})*$self->{param}->{gWidth} /$self->{param}->{xRange};
        my $y = $self->{param}->{gyBtm} -($yVal-$self->{param}->{yMin})*$self->{param}->{gHeight}/$self->{param}->{yRange};

        $self->{img}->line($prevx, $prevy, $x, $y, $color) if ($i > 0);

        $prevx = $x;
        $prevy = $y;
    }
}

#
# Plot the vertext triangle
#
sub PlotVertices
{
    my ($self, $color) = @_;

    my @vertexCoords = @{$self->{param}->{vtxCoords}};

    # Check if all vertices are in the drawing area
    my $xMin = $self->{param}->{xMin};
    my $xMax = $self->{param}->{xMax};
    my $yMin = $self->{param}->{yMin};
    my $yMax = $self->{param}->{yMax};
    my $rotAng = $self->{param}->{rotx} ? $self->{param}->{rotx} * $pi / 180 : 0;
    my $eurRoty = $afoVy + ($eurVy - $afoVy) * cos($rotAng);

    for my $popNo (0 .. 2) {
        my ($xVal, $yVal, $zVal, $junk) = @{$vertexCoords[$popNo]};
        $yVal = $eurRoty if ($popNo == 0);
        if ($xVal < $xMin || $xVal > $xMax || $yVal < $yMin || $yVal > $yMax) {
            return;
        }
    }

    my ($x1, $y1) = (0, 0);
    for my $i (0 .. 3) {
        my $popNo = $i % 3;
        my ($xVal, $yVal, $zVal, $junk) = @{$vertexCoords[$popNo]};
        $yVal = $eurRoty if ($popNo == 0);
        my $x2 = $self->{param}->{gxLeft}+($xVal-$xMin)*$self->{param}->{gWidth} /$self->{param}->{xRange};
        my $y2 = $self->{param}->{gyBtm} -($yVal-$yMin)*$self->{param}->{gHeight}/$self->{param}->{yRange};
        $self->{img}->line($x1, $y1, $x2, $y2, $color) if ($i > 0);

        $x1 = $x2;
        $y1 = $y2;
    }
}

#
# Get parameters for plotting the oval expected areas for populations
#
sub GetPopulationParameters
{
    my $self = shift;

    # Columns:
    # 1. beta: the slope of linear regression
    # 2. theta: the angle of the slope
    # 3 and 4: mean and SD of the x values after points are rotated
    # 5 and 6: mean and SD of the y values after points are rotated
    # 7. the length of the horiontal axis of the ellipse
    my @popParamStrs = (
        " 0.165663   0.164172  1.528709  0.008199  1.369508  0.008497  0.025285",
        " 0.481781   0.448967  1.050820  0.005359  1.100100  0.002581  0.011565",
        "-0.533563  -0.490136  1.662749  0.006226  1.424757  0.003582  0.015364",
        " 0.792775   0.670320  1.166605  0.067915  1.097353  0.005798  0.180400",
        " 0.689824   0.603863  1.455597  0.103179  1.100900  0.015606  0.263943",
        "-1.124799  -0.844065  1.219782  0.055360  1.609036  0.014449  0.143183",
        "-0.764712  -0.652851  1.580264  0.026192  1.519036  0.004822  0.062084",
        "-1.429868  -0.960496  1.133701  0.035680  1.620424  0.006446  0.077267");

    my @popParams = ();

    for my $popNo (0 ..7) {
        my $paramStr = $popParamStrs[$popNo];
        $paramStr =~ s/^\s+//;
        my @vals = split /\s+/, $paramStr;

        my %info = ("beta" => $vals[0], "theta" => $vals[1], "meanx" => $vals[2], "sdx" => $vals[3],
                    "meany" => $vals[4], "sdy" => $vals[5], "radius" => $vals[6]);
        push @popParams, \%info;
    }

    return @popParams;
}

#
# Plot cutoff lines to divide different populations
#
sub PlotCutoffLines
{
    my ($self, $color) = @_;

    # East Asian line inside the triangle
    my $easXval1 = $easVx - ($easVx - $afoVx) * (1 - $self->{param}->{easCut});
    my $easYval1 = $easVy;
    my $easXval2 = $easVx - ($easVx - $eurVx) * (1 - $self->{param}->{easCut});
    my $easYval2 = $easVy + ($eurVy - $easVy) * (1 - $self->{param}->{easCut});
    $self->PlotOneCutoffLine($easXval1, $easYval1, $easXval2, $easYval2, $color) if ($self->{param}->{easCut} < 1);

    # African line inside the triangle
    my $afoXval1 = $afoVx + ($easVx - $afoVx) * (1 - $self->{param}->{afoCut});
    my $afoYval1 = $afoVy;
    my $afoXval2 = $afoVx + ($eurVx - $afoVx) * (1 - $self->{param}->{afoCut});
    my $afoYval2 = $afoVy + ($eurVy - $afoVy) * (1 - $self->{param}->{afoCut});
    $self->PlotOneCutoffLine($afoXval1, $afoYval1, $afoXval2, $afoYval2, $color) if ($self->{param}->{afoCut} < 1);

    # European line inside the triangle
    my $eurXval1 = $eurVx - ($eurVx - $afoVx) * (1 - $self->{param}->{eurCut});
    my $eurYval1 = $eurVy - ($eurVy - $afoVy) * (1 - $self->{param}->{eurCut});
    my $eurXval2 = $eurVx + ($easVx - $eurVx) * (1 - $self->{param}->{eurCut});
    my $eurYval2 = $eurYval1;
    $self->PlotOneCutoffLine($eurXval1, $eurYval1, $eurXval2, $eurYval2, $color);

    # Cutoff lines outside the triangle, extended from the above three lines
    $self->PlotOuterCutoff($easXval1, $easYval1, 1, $color);
    $self->PlotOuterCutoff($afoXval1, $afoYval1, 1, $color);
    $self->PlotOuterCutoff($easXval2, $easYval2, 2, $color);
    $self->PlotOuterCutoff($eurXval1, $eurYval1, 3, $color);
    $self->PlotOuterCutoff($eurXval2, $eurYval2, 2, $color);
    $self->PlotOuterCutoff($afoXval2, $afoYval2, 3, $color);

    # Cutoff lines outside the triangle, extended from the other lines inside the triangle
    my $afoXval3 = $easVx - ($easVx - $afoVx) * (1 - $self->{param}->{othLatCut});
    my $afoXval4 = $easVx - ($easVx - $afoVx) * $self->{param}->{othLatCut};
    my $afoYval3 = $afoVy;
    my $afoYval4 = $afoVy;
    my $eurXval3 = $eurVx - ($eurVx - $afoVx) * $self->{param}->{afaLacCut};
    my $eurYval3 = $eurVy - ($eurVy - $afoVy) * $self->{param}->{afaLacCut};

    $self->PlotOuterCutoff($afoXval3, $afoYval3, 1, $color);
    $self->PlotOuterCutoff($afoXval4, $afoYval4, 1, $color);
    $self->PlotOuterCutoff($eurXval3, $eurYval3, 3, $color);

    # Cutoff lines inside the triangle, separating Hispanics from other populations
    my $ofXval1 = $afoXval3;
    my $ofYval1 = $afoVy;
    my $ofXval2 = $eurVx + ($easVx - $eurVx) * $self->{param}->{othLatCut};
    my $ofYval2 = $eurVy - ($eurVy - $easVy) * $self->{param}->{othLatCut};
    my $ofXval3 = $eurVx;
    my $ofYval3 = $afoVy + ($ofYval2 - $afoVy) * ($eurVx - $ofXval1) / ($ofXval2 - $ofXval1);

    my $oaXval1 = $afoXval4;
    my $oaYval1 = $afoVy;
    my $oaXval2 = $eurVx - ($eurVx - $afoVx) * $self->{param}->{othLatCut};
    my $oaYval2 = $eurVy - ($eurVy - $afoVy) * $self->{param}->{othLatCut};

    my $ofhXval = $ofXval2 - ($eurVx - $afoVx) * $self->{param}->{othLatCut};
    my $ofhYval = $ofYval2 - ($eurVy - $afoVy) * $self->{param}->{othLatCut};
    $self->PlotOneCutoffLine($oaXval1, $oaYval1, $ofhXval, $ofhYval, $color);
    $self->PlotOneCutoffLine($ofXval1, $ofYval1, $ofXval3, $ofYval3, $color);

    # Cutoff line separating the two Hispanic sub-populations
    my $afmVx = ($easVx + $afoVx) / 2;
    my $ehXval = $eurVx - ($eurVx - $afmVx) * (1 - $self->{param}->{eurCut});
    my $ehYval = $eurVy - ($eurVy - $afoVy) * (1 - $self->{param}->{eurCut});
    $self->PlotOneCutoffLine($eurVx, $ofYval3, $eurVx, $eurYval1, $color);

    # Cutoff lines separating Hispanics from African Americans
    my $fhXval1 = $eurXval3;
    my $fhYval1 = $eurYval3;
    my $fhXval2 = $eurXval3 + ($easVx - $eurVx) * $self->{param}->{othLatCut};
    my $fhYval2 = $eurYval3 - ($eurVy - $easVy) * $self->{param}->{othLatCut};
    $self->PlotOneCutoffLine($fhXval1, $fhYval1, $fhXval2, $fhYval2, $color);
}

#
# Plot cutoff line outside the triangle, starting from a certain point
#
sub PlotOuterCutoff
{
    my ($self, $x1, $y1, $popId, $color) = @_;

    my ($dx, $dy) = (0, 0);

    if ($popId == 1) {
        $dx = $x1 - $eurVx;
        $dy = $y1 - $eurVy;
    }
    elsif ($popId == 2) {
        $dx = $x1 - $afoVx;
        $dy = $y1 - $afoVy;
    }
    elsif ($popId == 3) {
        $dx = $x1 - $easVx;
        $dy = $y1 - $easVy;
    }

    my $dl = sqrt($dx**2 + $dy**2);
    my $x2 = $x1 + $dx * $cutLineLen / $dl;
    my $y2 = $y1 + $dy * $cutLineLen / $dl;

    $self->PlotOneCutoffLine($x1, $y1, $x2, $y2, $color);
}

#
# Plot one cutoff line given two points
#
sub PlotOneCutoffLine
{
    my ($self, $x1, $y1, $x2, $y2, $color) = @_;

    my $xP1 = $self->{param}->{gxLeft}+($x1-$self->{param}->{xMin})*$self->{param}->{gWidth} /$self->{param}->{xRange};
    my $yP1 = $self->{param}->{gyBtm} -($y1-$self->{param}->{yMin})*$self->{param}->{gHeight}/$self->{param}->{yRange};
    my $xP2 = $self->{param}->{gxLeft}+($x2-$self->{param}->{xMin})*$self->{param}->{gWidth} /$self->{param}->{xRange};
    my $yP2 = $self->{param}->{gyBtm} -($y2-$self->{param}->{yMin})*$self->{param}->{gHeight}/$self->{param}->{yRange};

    $self->{img}->line($xP1, $yP1, $xP2, $yP2, $color);
}

#
# Plot cutoff lines to separate SAS, ASN from other populations, on GD4 vs. GD1 graph
#
sub PlotSasCutoffLines
{
    my ($self, $color) = @_;

    # Plot cutoff line to separate South Asians from Europeans and other populations
    my $halfWid = 0.1; # from center to the rightmost x coord, determining how long to plot
    my $numPts = 1000;

    my ($x1, $y1) = (0, 0);
    for my $i (0 .. $numPts) {
        my $dx = ($i*2 - $numPts)/$numPts * $halfWid;
        my $xVal = $self->{param}->{meanSasx} + $dx;
        my $yVal = $self->{param}->{sasCutBasey} + $self->{param}->{sasCutaVal} * $dx * $dx;

        my $x2 = $self->{param}->{gxLeft}+($xVal-$self->{param}->{xMin})*$self->{param}->{gWidth} /$self->{param}->{xRange};
        my $y2 = $self->{param}->{gyBtm} -($yVal-$self->{param}->{yMin})*$self->{param}->{gHeight}/$self->{param}->{yRange};

        if ($i > 0 && $x2 > $self->{param}->{gxLeft}
                   && $x2 < $self->{param}->{gxRight}
                   && $y2 > $self->{param}->{gyTop}
                   && $y2 < $self->{param}->{gyBtm} ) {
            $self->{img}->line($x1, $y1, $x2, $y2, $color);
        }
        $x1 = $x2;
        $y1 = $y2;
    }

    # Plot cutoff line to separate East Asians from South Asians and other populations
    my $halfHt = 0.08; # determining the length of the segment being plotted
    ($x1, $y1) = (0, 0);
    for my $i (0 .. $numPts) {
        my $dy = ($i*2 - $numPts)/$numPts * $halfHt;
        my $yVal = $self->{param}->{meanAsny} + $dy;
        my $xVal = $self->{param}->{asnCutBasex} + $self->{param}->{asnCutaVal} * $dy * $dy;

        my $x2 = $self->{param}->{gxLeft}+($xVal-$self->{param}->{xMin})*$self->{param}->{gWidth} /$self->{param}->{xRange};
        my $y2 = $self->{param}->{gyBtm} -($yVal-$self->{param}->{yMin})*$self->{param}->{gHeight}/$self->{param}->{yRange};

        if ($i > 0 && $x2 > $self->{param}->{gxLeft}
                   && $x2 < $self->{param}->{gxRight}
                   && $y2 > $self->{param}->{gyTop}
                   && $y2 < $self->{param}->{gyBtm} ) {
            $self->{img}->line($x1, $y1, $x2, $y2, $color);
        }
        $x1 = $x2;
        $y1 = $y2;
    }
}

1;
