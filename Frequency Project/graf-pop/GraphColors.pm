#    ====================================================================================
#                               PUBLIC DOMAIN NOTICE
#                National Center for Biotechnology Information
#
#        This software/database is a "United States Government Work" under the
#        terms of the United States Copyright Act.  It was written as part of
#        the author's official duties as a United States Government employee and
#        thus cannot be copyrighted.  This software/database is freely available
#        to the public for use. The National Library of Medicine and the U.S.
#        Government have not placed any restriction on its use or reproduction.
#        Although all reasonable efforts have been taken to ensure the accuracy
#        and reliability of the software and data, the NLM and the U.S.
#        Government do not and cannot warrant the performance or results that
#        may be obtained by using this software or data. The NLM and the U.S.
#        Government disclaim all warranties, express or implied, including
#        warranties of performance, merchantability or fitness for any particular
#        purpose.
#
#        Please cite the author in any work or product based on this material.
#    ====================================================================================

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
