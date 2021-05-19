#!/usr/local/bin/perl -w

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
        File Description: script to extract genotypes of Ancestry SNPs from one or multiple vcf or vcf.gz files.
        Date: 05/06/2021
    ====================================================================================
EOF

use strict;
use warnings;
use Carp;
use Time::HiRes qw(gettimeofday);
use Cwd 'abs_path';
use File::Path qw(make_path remove_tree);

if (@ARGV < 2) {
    print "\n$disclaim\n";
    print "Usage: ExtractAncSnpsFromVcfGz.pl <vcf_or_vcf_gz_file> <output_vcf_file> [keyword]\n\n";
    print "Note:  When 'keyword' is provided, and the name of the input file includes the keyword followed\n";
    print "       by an integer, the script searches the directory for all files that have similar\n";
    print "       name as the input file, with the only difference being the integer after the keyword,\n";
    print "       and extract genotypes from all these files.\n";
    print "\n";
    exit;
}

my $script = $0;
my $pathToScript = abs_path($script);
if ($pathToScript =~ /(.+)\/(.+)/) {
    $pathToScript = $1;
    $script = $2;
}

my $vcfFile = $ARGV[0];
my $outFile = $ARGV[1];
my $keyWord = @ARGV > 2 ? $ARGV[2] : "";

my @t1 = gettimeofday;

unless (-e $vcfFile) {
    print "\nERROR: didn't find file $vcfFile!\n\n";
    exit;
}

my %allVcfFiles = ();

if ($keyWord) {
    my $shortFile = $vcfFile;
    my $fileDir = ".";
    if ($vcfFile =~ /(.+)\/(.+)/) {
        $fileDir = $1;
        $shortFile = $2;
    }

    # Make sure the keyword is included in the file name
    my $fileNo = 0;
    my ($prevStr, $postStr) = ("", "");

    if ($shortFile =~ /(\S+)$keyWord(\d+)(\S+)/) {
        $prevStr = $1;
        $fileNo = $2 + 0;
        $postStr = $3;
    }
    else {
        print "\nERROR: didn't find keyword '$keyWord' before a number in file name!\n\n";
        exit;
    }

    opendir DIR, $fileDir or die "Couldn't open directory $fileDir!\n";
    my @outFiles = readdir DIR;
    for my $file (@outFiles) {
        if ($file =~ /^$prevStr$keyWord(\d+)$postStr$/) {
            my $fileNo = $1 + 0;
            $allVcfFiles{$1} = "$fileDir/$file";
        }
    }
    closedir DIR;
}
else {
    $allVcfFiles{1} = $vcfFile;
}

my $numVcfFiles = keys %allVcfFiles;
print "Found $numVcfFiles vcf files with key word '$keyWord' before an integer.\n\n" if ($keyWord);

my $ancSnpFile = "AncInferSNPs.txt";
$ancSnpFile = "$pathToScript/data/$ancSnpFile";
unless (-e $ancSnpFile) {
    print "\nERROR: didn't find $ancSnpFile\n";
    exit;
}

my ($totAncSnps, $rsSnpRef, $gb37Ref, $gb38Ref) = GetAncestrySnps($ancSnpFile);
my %rsAncSnpIds = %$rsSnpRef;
my %gb37AncSnpIds = %$gb37Ref;
my %gb38AncSnpIds = %$gb38Ref;

print "Found $totAncSnps ancestry SNPs\n";
print "Extracting ancestry SNP genos from $numVcfFiles vcf files ...\n";

my ($totVcfSnps, $totSaveSnps) = (0, 0);
open OUTFILE, ">$outFile" or die "\nERROR: couldn't open $outFile for writing!\n";
my $saveHead = 1;
foreach my $fileNo (sort {$a <=> $b} keys %allVcfFiles) {
    print "File with $keyWord$fileNo: $allVcfFiles{$fileNo}\n" if ($numVcfFiles > 1);
    my ($numVcfSnps, $numVcfAncSnps) = ExtractAncGenoFromFile($allVcfFiles{$fileNo}, $saveHead);
    $saveHead = 0;
    $totVcfSnps += $numVcfSnps;
    $totSaveSnps += $numVcfAncSnps;
}
close OUTFILE;

print "\nExtracted genotypes of total $totSaveSnps SNPs from $totVcfSnps SNPs in $numVcfFiles files.\n";
print "Results saved to $outFile\n\n";

my @t2 = gettimeofday;
my $time = GetTimeDifference(\@t1, \@t2);
print "Time used $time\n";


sub GetAncestrySnps
{
    my $baseAncFile = shift;

    my $ancFile = $baseAncFile;

    my @allPaths = split /:/, $ENV{PATH};
    unshift @allPaths, "data";
    unless (-e $baseAncFile) {
        for my $path (@allPaths) {
            print "Checking $path\n";
            if (-e "$path/$baseAncFile") {
                $ancFile = "$path/$baseAncFile";
                last;
            }
    	}
    }
    die "\nERROR: didn't find ancestry SNP file $baseAncFile!\n" unless (-e $ancFile);
	
    my %rsSnpIds = ();
    my %gb37SnpIds = ();
    my %gb38SnpIds = ();

    my $snpId = 0;
    open FILE, $ancFile or die "Couldn't open $ancFile\n";
    my $head = <FILE>;
    while(<FILE>) {
        chomp;

        my ($chr, $gb37, $gb38, $rs, $ref, $alt, @freqs) = split /\s+/, $_;
        if ($chr && $gb37 && $gb38 && $rs) {
            $rsSnpIds{"rs$rs"} = $snpId;
            $gb37SnpIds{"$chr\t$gb37"} = $snpId;
            $gb38SnpIds{"$chr\t$gb38"} = $snpId;
        }

        $snpId++;
    }
    close FILE;

    return ($snpId, \%rsSnpIds, \%gb37SnpIds, \%gb38SnpIds);
}

sub ExtractAncGenoFromFile
{
    my ($file, $saveHead) = @_;

    my @saveLines = ();

    if ($file =~ /vcf\.gz/ || $file =~ /vcf\.bgz/) {
    	open FILE, "zcat $file |"  or die "gunzip $file: $!";
    }
    elsif ($file =~ /\.vcf$/) {
	    open FILE, $file or die "Couldn't open $file\n";
    }

    my $rowNo = 0;
    my $numAncSnps = 0;
    my $numVcfSnps = 0;
    while(<FILE>) {
        chomp;
        $rowNo++;

        if ($_ =~ /^#/) {
            push @saveLines, $_ if ($saveHead);
            next;
        }

        my ($chr, $pos, $snp) = split /\t/, $_;
        $chr =~ s/chr//;

        if (defined $rsAncSnpIds{$snp} ||
            defined $gb37AncSnpIds{"$chr\t$pos"} ||
            defined $gb38AncSnpIds{"$chr\t$pos"}) {
            push @saveLines, $_;
            $numAncSnps++;
        }

        $numVcfSnps++;
        print "\tChecked $numVcfSnps SNPs. Found $numAncSnps ancestry SNPs\n" if ($numVcfSnps % 5000000 == 0);
    }
    close FILE;
    print "\tChecked $numVcfSnps SNPs. Found $numAncSnps ancestry SNPs\n";

    if ($numAncSnps > 0) {
        for my $line (@saveLines) {
            print OUTFILE "$line\n";
        }
        print "\tSaved $numAncSnps FP SNPs\n";
    }

    return ($numVcfSnps, $numAncSnps);
}

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

    if ($ds == 0) {
        $msg .= "$dus micro seconds";
        return $msg;
    }

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
