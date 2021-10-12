#!/usr/local/bin/perl

###################################################################################################
#
#  Description: Runs GrafPop medium tests
#  Created by Jimmy Jin on 04/27/2021
#
###################################################################################################
use strict;

use Carp;
use Cwd 'abs_path';
use File::Path qw(make_path remove_tree);
use Getopt::Long::Descriptive;
use Time::HiRes qw(gettimeofday);

my $script = $0;
my $pathToScript = abs_path($script);
if ($pathToScript =~ /(.+)\/(.+)/) {
    $pathToScript = $1;
    $script = $2;
}

if (@ARGV < 1) {
    print "Usage $script <manifest> [is_update_baseline] [is_debug]\n\n";
    exit;
}

my $manifest = $ARGV[0];
my $updateBaseline = @ARGV > 1 ? $ARGV[1] : 0;
my $DEBUG = @ARGV > 2 ? $ARGV[2] : 0;

my $dataDir  = "$pathToScript/data";
my $inputDir  = "$pathToScript/input";
my $outputDir  = "$pathToScript/output";
my $baselineDir  = "$pathToScript/baseline";
my $grafpopDir = $pathToScript;
$grafpopDir =~ s/\/$//;
$grafpopDir =~ s/medium_testing$//;

my @t1 = gettimeofday;

my @allTests = GetAllTestCommands($manifest, $inputDir, $outputDir);

### Run the tests
my $numTests = @allTests;
my $numFailTests = 0;

print "\n";
if ($updateBaseline) {
    my $rmCmd = "rm $baselineDir/*.*";
    print "$rmCmd\n" if ($DEBUG);
    system($rmCmd);
}

my @failTests = ();
for my $testNo (0 .. $#allTests) {
    my %test = %{$allTests[$testNo]};

    my $testName = $test{name};
    my $testCmd = $test{command};
    my $testOutFile = $test{outputFile};
    my $expOutFile = $test{fileExpected};
    my $expMsg = $test{stdoutMsg};

    my $testId = $testNo+ 1;

    print "Running Test No. $testId: $testName\n";

    my ($testPassed, $err, $runTime)  = RunTest($testCmd, $testOutFile, $expOutFile, $expMsg, $outputDir);
    if ($testPassed) {
	print "\t-------------------------------> Completed. Time used: $runTime\n";
    }
    else {
	print "\t-------------------------------> Failed to run test: $err\n";
	push @failTests, "$testName\t$err";
	$numFailTests++;
    }
}

if ($numFailTests > 0) {
    print "\n\n*********************************  ERROR: $numFailTests tests failed! ********************************\n";
    for my $test (@failTests) {
	my ($name, $err) = split /\t/, $test;
	print "Test $name ERROR: $err\n";
    }
    print "\n";
}
else {
    print "\n\n#################################  Done. All $numTests test cases passed! #################################\n";
    print "Output files saved to directory $outputDir/\n";
}

my @t2 = gettimeofday;
my $totTime = GetTimeDifference(\@t1, \@t2);
print "\nTotal time used: $totTime\n\n";

#
#  Build the command and Execute the test case and check the results if necessary
#
sub RunTest
{
    my ($cmd, $testOutFile, $expOutFile, $expMsg, $testOutDir) = @_;

    my @testT1 = gettimeofday;
    my $fullCmd = "$cmd";
    $fullCmd =~ s/\-spf\s+(\S+)/\-spf $inputDir\/$1/;

    my $fullTestOutFile = "$testOutDir/$testOutFile";
    my $fullBaselineFile = "$baselineDir/$testOutFile";

    # Run test
    print "\t$fullCmd\n" if ($DEBUG);
    my $testResults = `$fullCmd`;

    my $passed = 0;
    my $error = "";

    my @outputLines = split /\n/, $testResults;
    for my $line (@outputLines) {
	if ($line =~ /(ERROR: .+)/ && $expMsg !~ /ERROR/) {
	    $error = $line;
	    $passed = 0;
	    last;
	}

	if ($expMsg =~ /^ERROR:*\s+(.+)/) {
	    my $msg = $1;

	    if ($line =~ /ERROR.+$msg/) {
		$passed = 1;
		last;
	    }
	}
	elsif ($line =~ /$expMsg/) {
	    $passed = 1;
	    last;
	}
    }
    $error = "Unexpeted STDOUT message" if (!$passed && !$error);

    if ($passed && $expOutFile) {
	if ($updateBaseline) {
	    my $cpCmd = "cp $fullTestOutFile $baselineDir";
	    print "$cpCmd\n" if ($DEBUG);
	    system($cpCmd);
	}
	else {
	    my $diffCmd = "diff $fullTestOutFile $fullBaselineFile";
	    print "$diffCmd\n" if ($DEBUG);
	    my $diff = `$diffCmd`;
	    if ($diff) {
		$passed = 0;
		$error = "Unexpected results found in output file";
	    }
	}
    }

    my @testT2 = gettimeofday;
    my $runTime = GetTimeDifference(\@testT1, \@testT2);

    return ($passed, $error, $runTime);
}

#
# Get all the test cases from the manifest and do some sanity checks
#
sub GetAllTestCommands
{
    my ($manifest, $inputDir, $outputDir) = @_;

    my @allTests = ();

    my $cmdNo = 1;
    open FILE, $manifest or die "\nERROR: Couldn't open $manifest!\n";
    while(<FILE>) {
        chomp;
        my $line = $_;
        $line =~ s/\s*$//;
        next if ($line =~ /^\s*\#/ || $line !~ /\S/);

        my @vals = split /\s*\|\s*/, $line;
        for my $i (0 .. $#vals) {
            $vals[$i] = s/^\s*//;
            $vals[$i] = s/\s*$//;
        }

        my ($testName, $cmd, $inputFile, $raceFile, $options, $outputFile, $expFile, $expMsg) = split /\s*\|\s*/, $line;

	die "\nERROR: didn't find input file in No. $cmdNo: $testName\n" unless ($inputFile);
	die "\nERROR: didn't find output file in No. $cmdNo: $testName\n" unless ($outputFile);

	$cmd .= " $inputDir/$inputFile";
	$cmd .= " $outputDir/$outputFile";
	$cmd .= " -spf $raceFile" if ($raceFile);
	$cmd .= " $options" if ($options);

	my %test = (name => $testName,
		    command => $cmd,
		    inputFile => $inputFile,
		    raceFile => $raceFile,
		    options => $options,
		    outputFile => $outputFile,
		    fileExpected => $expFile,
		    stdoutMsg => $expMsg);
	push @allTests, \%test;

        $cmdNo++;
    }
    close FILE;

    return @allTests;
}

#
# Return a string showing time difference between two times.
# $t1 and $t2 are references to return thisLine of 'gettimeofday'.
#
sub GetTimeDifference
{
    my ($t1, $t2) = @_;

    my $msg = "";

    my $t1sec = $$t1[0];
    my $t1us  = $$t1[1];
    my $t2sec = $$t2[0];
    my $t2us  = $$t2[1];

    my $ds = $t2sec - $t1sec;
    my $dus = $t2us - $t1us;

    if ($dus < 0) {
        $dus += 1000000;
        $ds -= 1;
    }

    if ($ds < 0 || ($ds == 0 && $dus < 0)) {
        print "Error: T1 > T2\n";
        return $msg;
    }

    my $dds = $ds % 60;
    $dds += $dus/1000000;

    my $dh = int($ds/3600);
    $ds = $ds % 3600;
    my $dm = int($ds/60);

    $msg .= "$dh hour" if ($dh > 0);
    $msg .= "s" if ($dh > 1);
    $msg .= " $dm minute" if ($dh > 0 || $dm > 0);
    $msg .= "s" if ($dm > 1);
    $msg .= " $dds seconds";

    return $msg;
}
