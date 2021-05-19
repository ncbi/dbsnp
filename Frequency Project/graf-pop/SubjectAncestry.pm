package SubjectAncestry;

use strict;
use Carp;
use GraphParameters;

sub new
{
    my ($class, $param, $sbjScores, $sbjRaces) = @_;
    my $refType = ref($sbjScores);

    die "\nERROR: failed to create SubjectAncestry: sbjScores is not an array\n" if ($refType ne "ARRAY");
    my $numSbjs = $#$sbjScores + 1;
    print "\nWARNING: sbjScores is an empty\n" if ($numSbjs < 1);
    if ($sbjRaces) {
        die "\nERROR: SubjectAncestry: sbjRaces is not a hash\n" if (ref($sbjRaces) ne "HASH");
        foreach my $sbjNo (0 .. $#$sbjScores) {
            my $sbj = $sbjScores->[$sbjNo]->{subject};
            $sbjScores->[$sbjNo]->{race} = $sbjRaces->{$sbj} if ($sbjRaces->{$sbj});
        }
    }

    bless {
        param => $param,
        sbjGenoPopIds => [],
        numSubjects => $numSbjs,
        sbjScores => $sbjScores
    }, $class;
}

sub ShowPopulationComparison
{
    my $self = shift;

    my %raceGenoIds = ();
    my @sbjPopScores = @{$self->{sbjScores}};
    my %raceSbjs = ();
    my %raceGenoPopCnts = ();

    for my $sbjNo (0 .. $#sbjPopScores) {
        my %info = %{$sbjPopScores[$sbjNo]};
        my $sbj  = $info{subject};
        my $race = $info{race};
        my $genoPopId = $self->{sbjGenoPopIds}->[$sbjNo];
        if ($raceSbjs{$race}) {
            $raceSbjs{$race}++;
            $raceGenoPopCnts{$race}->[$genoPopId]++;
        }
        else {
            $raceSbjs{$race} = 1;
            my @popCnts = ();
            for my $popId (0 .. 9) {
                push @popCnts, 0;
            }
            $popCnts[$genoPopId] = 1;
            $raceGenoPopCnts{$race} = \@popCnts;
        }
    }

    my $maxRaceLen = 0;
    foreach my $race (keys %raceSbjs) {
        my $len = length($race);
        $maxRaceLen = $len if ($len > $maxRaceLen);
    }
    $maxRaceLen = 5 if ($maxRaceLen < 5);

    print "\nThe following table shows the self-reported races/ethnicities (column 'Race')\n" .
          "and population IDs assigned by GrafPop (other columns, Pop9 = Other)\n\n";
    my $cntLen = 6;
    my $raceFmt = "\%$maxRaceLen" . "s";
    my $cntFmt = "\%$cntLen" . "d";

    my $header = sprintf($raceFmt, "Race");
    for my $popId (1 .. 9) {
        $header = sprintf("%s  \%$cntLen" . "s", $header, "Pop$popId");
    }
    $header = sprintf("%s  \%$cntLen" . "s", $header, "Total");
    print "$header\n";
    my $headLen = length($header);
    my $dashLine = "-";
    for my $i (1 .. $headLen) {$dashLine .= "-"};
    print "$dashLine\n";
    foreach my $race (sort {$raceSbjs{$b} <=> $raceSbjs{$a}} keys %raceSbjs) {
        my @popCnts = @{$raceGenoPopCnts{$race}};
        my $line = sprintf($raceFmt, $race);
        for my $popNo (1 .. 9) {
            $line = sprintf("%s  $cntFmt", $line, $popCnts[$popNo]);
        }
        $line = sprintf("%s  $cntFmt", $line, $raceSbjs{$race});
        print "$line\n";
    }
    print "$dashLine\n";
}

#
# Assign populations to subjects based on the ancestry scores
#
sub SetSubjectGenoPopulations
{
    my ($self) = @_;

    my @sbjPopScores = @{$self->{sbjScores}};
    my %subPopSbjs = ();
    $self->{sbjGenoPopIds} = [];

    for my $sbjNo (0 .. $#sbjPopScores) {
        my %info = %{$sbjPopScores[$sbjNo]};
        my $gd1  = $info{x};
        my $gd4  = $info{gd4};
        my $ePct = $info{ePct};
        my $fPct = $info{fPct};
        my $aPct = $info{aPct};
        my $numSnps = $info{snps};

        # Calculate genotype populations based on ancestry proportions
        my $genoPopId = 0;
        if ($numSnps > 1000) {
            $genoPopId = $self->GetGenoPopId($gd1, $gd4, $ePct, $fPct, $aPct);
        }

        push @{$self->{sbjGenoPopIds}}, $genoPopId;
    }
}

#
# Save subjects and ancestry results to the output file
#
sub SaveResults
{
    my ($self, $outFile) = @_;

    my @sbjOutputLines = ();
    my $numShowSbjs = 0;
    my @sbjOutputLines = ();
    my @sbjPopScores = @{$self->{sbjScores}};
    my @sbjGenoPopIds = @{$self->{sbjGenoPopIds}};
    my %subPopSbjs = ();

    for my $sbjNo (0 .. $#sbjPopScores) {
        my %info = %{$sbjPopScores[$sbjNo]};
        my $sbj  = $info{subject};
        my $snps = $info{snps};
        my $race = $info{race};
        my $gd1  = $info{x};
        my $gd2  = $info{y};
        my $gd3  = $info{z};
        my $gd4  = $info{gd4};
        my $ePct = $info{ePct};
        my $fPct = $info{fPct};
        my $aPct = $info{aPct};
        my $numSnps = $info{snps};
        my $genoPopId = $sbjGenoPopIds[$sbjNo];

        # Check if the subject is within the selected area or not
        my $isWithin = 1;
        $isWithin = 0 if ($self->{param}->{xCutMin} && $gd1 < $self->{param}->{xCutMin});
        $isWithin = 0 if ($self->{param}->{yCutMin} && $gd2 < $self->{param}->{yCutMin});
        $isWithin = 0 if ($self->{param}->{xCutMax} && $gd1 > $self->{param}->{xCutMax});
        $isWithin = 0 if ($self->{param}->{yCutMax} && $gd2 > $self->{param}->{yCutMax});

        # Decide whether to include the subject or not based on the selected area
        my $showData = 1;
        $showData = 0 if (($self->{param}->{isBeyond} && $isWithin) || (!$self->{param}->{isBeyond} && !$isWithin));
        if ($self->{param}->{isSelAll} || $showData) {
            my $showGd1  = sprintf("%5.4f", $gd1);
            my $showGd2  = sprintf("%5.4f", $gd2);
            my $showGd3  = sprintf("%6.5f", $gd3);
            my $showGd4  = sprintf("%6.5f", $gd4);
            my $showbPct = sprintf("%5.2f", $fPct);
            my $showwPct = sprintf("%5.2f", $ePct);
            my $showaPct = sprintf("%5.2f", $aPct);
            my $genoPop  = "";
            if ($genoPopId > 0 && $genoPopId < 9) {
                my $genoPopNo = $genoPopId - 1;
                $genoPop = $self->{param}->{alfaFullPops}->[$genoPopNo];
            }
            push @sbjOutputLines,
                "$sbj\t$numSnps\t$race\t$showGd1\t$showGd2\t$showGd3\t$showGd4\t$showbPct\t$showwPct\t$showaPct\t$genoPopId\t$genoPop";

            if ($subPopSbjs{$genoPopId}) {
                $subPopSbjs{$genoPopId}++;
            }
            else {
                $subPopSbjs{$genoPopId} = 1;
            }
            $numShowSbjs++;
        }
    }

    if ($numShowSbjs > 0) {
        open FILE, ">$outFile" or die "\nERROR: Couldn't open $outFile for writing!\n";
        print FILE "Subject\t#SNPs\tSelf-reported ancestry\tGD1\tGD2\tGD3\tGD4\tP_f (%)\tP_e (%)\tP_a (%)\tPopID\tComputed population\n";
        for my $line (@sbjOutputLines) {
            print FILE "$line\n";
        }
        close FILE;
        print "\nTotal $self->{numSubjects} subjects. $numShowSbjs subjects were selected and saved to $outFile\n\n";

        print "\tPopID\t#Subjs\tPopulation\n";
        foreach my $genoPopId (1 .. 8) {
            my $genoPop = $self->{param}->{alfaFullPops}->[$genoPopId-1];
            my $numSbjs = $subPopSbjs{$genoPopId};
            printf("\t%d\t%6d\t%s\n", $genoPopId, $numSbjs, $genoPop);
        }
        print "\n";
    }
    else {
        print "\nNo subjects are found in the specified area.\n\n";
    }
}

#
# Assign population based on the calculated ancestry protortions and GD4
#
sub GetGenoPopId
{
    my ($self, $gd1, $gd4, $ePct, $fPct, $aPct) = @_;
    my $eComp = $ePct / 100;
    my $fComp = $fPct / 100;
    my $aComp = $aPct / 100;
    my $eurVx = $self->{param}->{vtxCoords}->[0]->[0];

    my $genoPopId = 0;

    my $isSas = 0;
    my $meanSasx = $self->{param}->{meanSasx};
    my $meanAsny = $self->{param}->{meanAsny};

    my $y = $self->{param}->{sasCutBasey} + $self->{param}->{sasCutaVal} * ($gd1-$meanSasx)**2;
    $isSas = 1 if ($gd4 > $y);

    my $isAsn = 0;
    my $x = $self->{param}->{asnCutBasex} + $self->{param}->{asnCutaVal} * ($gd4-$meanAsny)**2;
    $isAsn = 1 if ($gd1 > $x);

    if ($eComp > $self->{param}->{eurCut}) {
    	$genoPopId = 1;
    }
    elsif ($fComp > $self->{param}->{afoCut}) {
	    $genoPopId = 2;
    }
    elsif ($aComp > $self->{param}->{easCut}) {
	    $genoPopId = 3;
    }
    elsif ($isSas) {
	    $genoPopId = 8;
    }
    elsif ($fComp < $self->{param}->{othLatCut}) {
        if ($isAsn) {
            $genoPopId = 7;
        }
        else {
            if ($gd1 < $eurVx && $aComp < $self->{param}->{othLatCut}) {
                $genoPopId = 5;
            }
            else {
                if ($gd4 + $gd1 < $self->{param}->{asnLatCut}) {
                    $genoPopId = 6; # Hispanic2 are on the left lower side
                }
                else {
                    # Set pop ID to other. These are not Hispanics,
                    # probably European/Asian admixtures. Set to Other.
                    $genoPopId = 9;
                }
            }
        }
    }
    elsif ($aComp < $self->{param}->{othLatCut}) {
        if ($fComp > $self->{param}->{afaLacCut}) {
            $genoPopId = 4;
        }
        else {
            $genoPopId = 5;
        }
    }
    else {
    	$genoPopId = 9;
    }

    return $genoPopId;
}

1;
