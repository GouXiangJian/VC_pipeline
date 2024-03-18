#!/usr/bin/env perl

#date:   2021-11-27
#writer: Xiangjian Gou
#func:   filter SNP according to the minor allele frequency (MAF) and missing rate (MR)
#usage:  perl snpFilter.pl input_file MAF MR output_file_prefix
#e.g.:   perl snpFilter.pl barley.wgs.477.hmp 0.05 0.2 barley.wgs.477

use strict;
use warnings;

$| = 1;

my ($inputFile, $mafThreshold, $mrThreshold, $outputFilePrefix) = @ARGV;

open my $iHMP, '<', $inputFile or die "Error: cannot open file '$inputFile': $!";
open my $oHMP, '>', "$outputFilePrefix-MAF$mafThreshold-MR$mrThreshold.hmp";
open my $oLOG, '>', "$outputFilePrefix-MAF$mafThreshold-MR$mrThreshold.log";
print $oHMP scalar <$iHMP>;

while (<$iHMP>) {
    my @row = split;
    my @snp = @row[11 .. $#row];
    my %baseCount;
    foreach my $base (@snp) {
        if ($base =~ /\A[AGTC]\z/) {
            $baseCount{$base}++;
        }
    }
    my $genotypeCount = keys %baseCount;
    if ($genotypeCount == 2) {
        my ($maxBase, $minBase) = sort { $baseCount{$b} <=> $baseCount{$a} } keys %baseCount;
        my $maxCount = $baseCount{$maxBase};
        my $minCount = $baseCount{$minBase};
        my $MAF = $minCount/($maxCount+$minCount); #MAF: exclude N and heterozygous base
        my $MR = (@snp-$maxCount-$minCount)/@snp; #MR: consider both N and heterozygous base
        if ($MAF >= $mafThreshold and $MR <= $mrThreshold) {
            print $oHMP $_;
        }
    }
    else {
        print $oLOG "$row[0]\t$genotypeCount\n"; #snp id, genotype count
    }
}

close $iHMP;
close $oHMP;
close $oLOG;
