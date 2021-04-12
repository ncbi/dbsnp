package GraphColors;

use strict;
use Carp;
use GD;

sub new
{
    my ($class, $img) = @_;

    my $white   = $img->colorAllocate(255, 255, 255);
    my $red     = $img->colorAllocate(255,   0,   0);
    my $maroon  = $img->colorAllocate(255,   0,   0);
    my $green   = $img->colorAllocate(255, 128,   0); # not green anymore
    my $blue    = $img->colorAllocate(  0,   0, 255);
    my $navy    = $img->colorAllocate(  0,   0, 128);
    my $black   = $img->colorAllocate(  0,   0,   0);
    my $gray    = $img->colorAllocate(128, 128, 128);
    my $yellow  = $img->colorAllocate(255, 255,   0);
    my $olive   = $img->colorAllocate(128, 128,   0);
    my $purple  = $img->colorAllocate(128,   0, 128);
    my $magenta = $img->colorAllocate(255,   0, 255);
    my $orange  = $img->colorAllocate(255, 165,   0);
    my $cyan    = $img->colorAllocate(  0, 255, 255);
    my $teal    = $img->colorAllocate(  0, 128, 128);
    my $gold    = $img->colorAllocate(204, 153,  80);
    my $expAreaColor = $img->colorAllocate(  0, 150, 150);

    # Colors for self-reported ancestries
    my $maxShowRaces = 10;
    my @raceColors = ();
    for my $raceNo (1 .. $maxShowRaces) {
        push @raceColors, "";
    }
    $raceColors[0] = $yellow;
    $raceColors[1] = $blue;
    $raceColors[2] = $red;
    $raceColors[3] = $olive;
    $raceColors[4] = $purple;
    $raceColors[5] = $cyan;
    $raceColors[6] = $green;
    $raceColors[7] = $teal;
    $raceColors[8] = $navy;
    $raceColors[9] = $magenta;
    $raceColors[10] = $gold;

    bless {
        white   => $white,
	    red     => $red,
        maroon  => $maroon,
        green   => $green,
        blue    => $blue,
        navy    => $navy,
        black   => $black,
        gray    => $gray,
        olive   => $olive,
        magenta => $magenta,
        orange  => $orange,
        cyan    => $cyan,
        teal    => $teal,
        gold    => $gold,
        expAreaColor => $expAreaColor,
        maxShowRaces => $maxShowRaces,
        raceColors => \@raceColors
    }, $class;
}

1;
