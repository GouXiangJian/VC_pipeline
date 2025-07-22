#!/usr/bin/env perl

#date:   2025-07-22
#writer: Xiangjian Gou (QQ = 862137261)

#load modules
use strict;
use warnings;
use Getopt::Long;

#record version information
my $VERSION = 'variant (SNP/INDEL) calling, step4 v1.2 (2025-07-22)';

#set default options
my $refGenomeFa = "/public/home/xjgou/genome/genome.fa";
my $softwareDir = "/public/home/xjgou/tools";
my $outputDir = "output";
my $jumpGenomicsDBImport = '';
my $jumpGenotypeGVCFs = '';
my $queue = "smp";
my $thread = 10;
my $memory = 100;
my $version;
my $help;

#get options from command line
GetOptions(
    'refGenomeFa=s'            => \$refGenomeFa,
    'softwareDir=s'            => \$softwareDir,
    'outputDir=s'              => \$outputDir,
    'jumpGenomicsDBImport|jm'  => \$jumpGenomicsDBImport,
    'jumpGenotypeGVCFs|jt'     => \$jumpGenotypeGVCFs,
    'queue=s'                  => \$queue,
    'thread=i'                 => \$thread,
    'memory=i'                 => \$memory,
    'version+'                 => \$version,
    'help+'                    => \$help,
);

#describe program information
my $usage = <<__GUIDE__;
####################################################################################################
Function: genotyping for each chromosome, including: GenomicsDBImport and GenotypeGVCFs

Usage: perl step4.pl option1 <value1> option2 <value2> ... optionN <valueN>

Options:
  #Options for path:
  -r  | -refGenomeFa  <STRING> : the reference genomic fasta file. (default: /public/home/xjgou/genome/genome.fa)
  -s  | -softwareDir  <STRING> : the directory of storing all software. (default: /public/home/xjgou/tools)
  -o  | -outputDir    <STRING> : set a directory for storing output information. (default: output)

  #Options for step:
  -jm | -jumpGenomicsDBImport  : no execute GenomicsDBImport
  -jt | -jumpGenotypeGVCFs     : no execute GenotypeGVCFs

  #Options for resources:
  -q  | -queue  <STRING>       : set the queue to use. (default: smp)
  -t  | -thread    <INT>       : set the number of threads to use. (default: 10)
  -m  | -memory    <INT>       : set the size of memory to use. (default: 100 [100GB])

  #Options for other:
  -v  | -version               : show the version information.
  -h  | -help                  : show the help information.
####################################################################################################

__GUIDE__

#output version and help information
die "$VERSION\n" if $version;
die $usage if $help;

#set whether to execute each step
$jumpGenomicsDBImport = '#' if $jumpGenomicsDBImport;
$jumpGenotypeGVCFs = '#' if $jumpGenotypeGVCFs;

#get the all chr name
open my $iFASTA, '<', $refGenomeFa or die "Error: cannot open file '$refGenomeFa': $!";
my @names;
while (<$iFASTA>) {
    if (/\A>/) {
        my ($name) = /\A>(\S+)/;
        push @names, $name;
    }
}
close $iFASTA;

#generate a scheduling script for each chromosome, meanwhile, generate a comprehensive scheduling script
open my $oTOTAL, '>', 'step4.sh';
my $lsfDir = "step4.lsf";
system "mkdir -p $lsfDir";
foreach my $chr (@names) {
    open my $oEACH, '>', "$lsfDir/$chr.lsf";
    lsfInfo($oEACH, $chr, $queue, $thread, $memory, $outputDir, $softwareDir, $refGenomeFa, $jumpGenomicsDBImport, $jumpGenotypeGVCFs);
    close $oEACH;
    print $oTOTAL "bsub < $lsfDir/$chr.lsf\n";
}
close $oTOTAL;

#create a subroutine to write all command into lsf script
sub lsfInfo {
    my ($handle, $chr, $queue, $thread, $memory, $outputDir, $softwareDir, $refGenomeFa, $jumpGenomicsDBImport, $jumpGenotypeGVCFs) = @_;
    my $input = join ' ', map {"-V $_"} grep {/vcf\z/} glob "$outputDir/HaplotypeCallerOut/*";
    my $command = <<__COMMAND__;
#!/bin/bash

#BSUB -q $queue
#BSUB -n $thread
#BSUB -J step4.$chr
#BSUB -o step4.$chr.out
#BSUB -e step4.$chr.err

#run GenomicsDBImport
$jumpGenomicsDBImport mkdir -p $outputDir/GenomicsDBImportOut
$jumpGenomicsDBImport $softwareDir/gatk-4.1.2.0/gatk --java-options "-Xmx${memory}g -Xms${memory}g" GenomicsDBImport --genomicsdb-workspace-path $outputDir/GenomicsDBImportOut/$chr --batch-size 50 -R $refGenomeFa -L $chr $input

#run GenotypeGVCFs
$jumpGenotypeGVCFs mkdir -p $outputDir/GenotypeGVCFsOut
$jumpGenotypeGVCFs $softwareDir/gatk-4.1.2.0/gatk --java-options "-Xmx${memory}g -Xms${memory}g" GenotypeGVCFs -R $refGenomeFa -V gendb://$outputDir/GenomicsDBImportOut/$chr -O $outputDir/GenotypeGVCFsOut/$chr.vcf -new-qual -G StandardAnnotation --use-new-qual-calculator

__COMMAND__
    print $handle $command;
}
