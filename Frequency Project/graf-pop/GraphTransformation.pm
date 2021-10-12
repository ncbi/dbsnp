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

package GraphTransformation;

use strict;
use Carp;

my $pi = 3.1415926535;

#
# Move and rotate a 2D point (given x, y coordinates
#
sub TransformPoint2D
{
    my ($x1, $y1, $dx, $dy, $alpha) = @_;

    my $xt = $x1 + $dx;
    my $yt = $y1 + $dy;

    my $x = $xt * cos($alpha) - $yt * sin($alpha);
    my $y = $yt * cos($alpha) + $xt * sin($alpha);

    $x = $x - $dx;
    $y = $y - $dy;

    return ($x, $y);
}

#
# Rotate a 3D point (array reference) around an axis
#
sub RotatePoint3D
{
    my ($pRef, $deg, $axis) = @_; # deg = angle in degree; axis = x, y or z

    my $rad = $deg * $pi / 180;
    my $sind = sin($rad);
    my $cosd = cos($rad);

    my $t1 = [];
    if ($axis =~ /^x/) {
        $t1 = [
                [1, 0,      0,     0],
                [0, $cosd, -$sind, 0],
                [0, $sind,  $cosd, 0],
                [0, 0,      0,     1]  ];
    }
    elsif ($axis =~ /^y/) {
        $t1 = [
                [ $cosd, 0,  $sind, 0],
                [ 0,     1,  0,     0],
                [-$sind, 0,  $cosd, 0],
                [ 0,     0,  0,     1]  ];
    }
    elsif ($axis =~ /^z/) {
        $t1 = [
                [$cosd, -$sind, 0, 0],
                [$sind,  $cosd, 0, 0],
                [0,      0,     1, 0],
                [0,      0,     0, 1]  ];
    }

    my @p2 = TransformPoint3D($pRef, $t1);

    return @p2;
}

#
# Move a 3D point
#
sub MovePoint3D
{
    my ($pRef, $dx, $dy, $dz) = @_;

    my $t1 = [  [1, 0, 0, $dx],
                [0, 1, 0, $dy],
                [0, 0, 1, $dz],
                [0, 0, 0, 1  ]  ];

    my @p2 = TransformPoint3D($pRef, $t1);

    return @p2;
}

#
# Transform a 3D point given a transformation matrix
#
sub TransformPoint3D
{
    my ($pRef, $tRef) = @_;

    my @p = @$pRef;
    my @t = @$tRef;

    my @q = (0, 0, 0, 0);

    for my $i (0 .. 3) {
        for my $j (0 .. 3) {
            $q[$i] += $t[$i]->[$j] * $p[$j];
        }
    }

    return @q;
}

#
# Rotate all points in the global array around an axis
#
sub RotateShape3D
{
    my ($deg, $axis, $dataPts) = @_;

    my $rad = $deg * $pi / 180;
    my $sind = sin($rad);
    my $cosd = cos($rad);

    my $t1 = [];
    if ($axis =~ /x/) {
        $t1 = [
                [1, 0,      0,     0],
                [0, $cosd, -$sind, 0],
                [0, $sind,  $cosd, 0],
                [0, 0,      0,     1]  ];
    }
    elsif ($axis =~ /y/) {
        $t1 = [
                [ $cosd, 0,  $sind, 0],
                [ 0,     1,  0,     0],
                [-$sind, 0,  $cosd, 0],
                [ 0,     0,  0,     1]  ];
    }
    elsif ($axis =~ /z/) {
        $t1 = [
                [$cosd, -$sind, 0, 0],
                [$sind,  $cosd, 0, 0],
                [0,      0,     1, 0],
                [0,      0,     0, 1]  ];
    }

    TransformShape3D($t1, $dataPts);
}

#
# Scale all points in the global array
#
sub ScaleShape3D
{
    my ($sx, $sy, $sz, $dataPts) = @_;
    $sy = $sx unless ($sy);
    $sz = $sx unless ($sz);

    my $t1 = [  [$sx, 0,   0,   0],
                [0,   $sy, 0,   0],
                [0,   0,   $sz, 0],
                [0,   0,   0,   1]  ];

    TransformShape3D($t1, $dataPts);
}

#
# Move all points in the global array
#
sub MoveShape3D
{
    my ($dx, $dy, $dz, $dataPts) = @_;
    my $t1 = [  [1, 0, 0, $dx],
                [0, 1, 0, $dy],
                [0, 0, 1, $dz],
                [0, 0, 0, 1  ]  ];

    TransformShape3D($t1, $dataPts);
}


#
# Transform a list of points
#
sub TransformShape3D
{
    my ($tRef, $dataPts) = @_;

    for my $i (0 .. $#$dataPts) {
        my @p = @{$dataPts->[$i]};
        my @q = TransformPoint3D(\@p, $tRef);
        $dataPts->[$i] = \@q;
    }
}


1;
