#!/usr/bin/env perl

#date:   2021-11-07
#writer: Xiangjian Gou (QQ = 862137261)

#load modules
use strict;
use warnings;
use Getopt::Long;

#record version information
my $VERSION = 'variant (SNP/INDEL) calling, step1 v1.0 (2021.11.07)';

#set default options
my $germplasmDir = "/public/home/xjgou/base_composition/six_sample/rawReads";
my $refGenomeLib = "/public/home/xjgou/base_composition/six_sample/barleyv2/barleyv2";
my $refGenomeFa = "/public/home/xjgou/base_composition/six_sample/barleyv2/barleyv2.fa";
my $softwareDir = "/public/home/xjgou/tools";
my $outputDir = "output";
my $jumpNGSQC = '';
my $jumpBWA = '';
my $jumpSortSam = '';
my $jumpMarkDuplicates = '';
my $jumpHaplotypeCaller = '';
my $queue = "normal";
my $thread = 4;
my $memory = 40;
my $version;
my $help;

#get options from command line
GetOptions(
    'germplasmDir=s'          => \$germplasmDir,
    'refGenomeLib|rl=s'       => \$refGenomeLib,
    'refGenomeFa|rf=s'        => \$refGenomeFa,
    'softwareDir=s'           => \$softwareDir,
    'outputDir=s'             => \$outputDir,
    'jumpNGSQC|jn+'           => \$jumpNGSQC,
    'jumpBWA|jb+'             => \$jumpBWA,
    'jumpSortSam|js+'         => \$jumpSortSam,
    'jumpMarkDuplicates|jm+'  => \$jumpMarkDuplicates,
    'jumpHaplotypeCaller|jh+' => \$jumpHaplotypeCaller,
    'queue=s'                 => \$queue,
    'thread=i'                => \$thread,
    'memory=i'                => \$memory,
    'version+'                => \$version,
    'help+'                   => \$help,
);

#describe program information
my $usage = <<__GUIDE__;
####################################################################################################
Function: data preprocessing, including: NGSQC, BWA, SortSam, MarkDuplicates, and HaplotypeCaller (non-BQSR)

Usage: perl step1.pl option1 <value1> option2 <value2> ... optionN <valueN>

Options:
  #Options for path:
  -g  | -germplasmDir <STRING> : the directory of storing all germplasms. (default: /public/home/xjgou/base_composition/six_sample/rawReads)
  -rl | -refGenomeLib <STRING> : the reference genomic library. (default: /public/home/xjgou/base_composition/six_sample/barleyv2/barleyv2)
  -rf | -refGenomeFa  <STRING> : the reference genomic fasta file. (default: /public/home/xjgou/base_composition/six_sample/barleyv2/barleyv2.fa)
  -s  | -softwareDir  <STRING> : the directory of storing all software. (default: /public/home/xjgou/tools)
  -o  | -outputDir    <STRING> : set a directory for storing output information. (default: output)

  #Options for step:
  -jn | -jumpNGSQC             : no execute NGSQC
  -jb | -jumpBWA               : no execute BWA
  -js | -jumpSortSam           : no execute SortSam
  -jm | -jumpMarkDuplicates    : no execute MarkDuplicates
  -jh | -jumpHaplotypeCaller   : no execute HaplotypeCaller (non-BQSR)

  #Options for resources:
  -q  | -queue  <STRING>       : set the queue to use. (default: normal)
  -t  | -thread    <INT>       : set the number of threads to use. (default: 4)
  -m  | -memory    <INT>       : set the size of memory to use. (default: 40 [40GB])

  #Options for other:
  -v  | -version               : show the version information.
  -h  | -help                  : show the help information.

Notes:
  (a) Each germplasm directory can only contain two files (dirName-1.fq.gz, dirName-2.fq.gz).
  (b) The bundled software is recommended.
  (c) Remember to install the Perl module 'String::Approx' before using NGSQC.
  (d) Remember to build a genomic library and create index (.fai/.dict) before using the script.
      e.g.
        \$ bwa index -p genome genome.fa
        \$ samtools faidx genome.fa
        \$ gatk CreateSequenceDictionary -R genome.fa -O genome.dict
####################################################################################################

__GUIDE__

#output version and help information
die "$VERSION\n" if $version;
die $usage if $help;

#set whether to execute each step
$jumpNGSQC = '#' if $jumpNGSQC;
$jumpBWA = '#' if $jumpBWA;
$jumpSortSam = '#' if $jumpSortSam;
$jumpMarkDuplicates = '#' if $jumpMarkDuplicates;
$jumpHaplotypeCaller = '#' if $jumpHaplotypeCaller;

#generate a scheduling script for each germplasm, meanwhile, generate a comprehensive scheduling script
open my $oTOTAL, '>', 'step1.sh';
my $lsfDir = "step1.lsf";
system "mkdir -p $lsfDir";
foreach my $dir (glob "$germplasmDir/*") {
    my ($germplasm) = $dir =~ /([^\/]+)\z/;
    open my $oEACH, '>', "$lsfDir/$germplasm.lsf";
    lsfInfo($oEACH, $germplasm, $germplasmDir, $outputDir, $softwareDir, $queue, $thread, $memory, $refGenomeLib, $refGenomeFa, $jumpNGSQC, $jumpBWA, $jumpSortSam, $jumpMarkDuplicates, $jumpHaplotypeCaller);
    close $oEACH;
    print $oTOTAL "bsub < $lsfDir/$germplasm.lsf\n";
}
close $oTOTAL;

#create a subroutine to write all command into lsf script
sub lsfInfo {
    my ($handle, $germplasm, $germplasmDir, $outputDir, $softwareDir, $queue, $thread, $memory, $refGenomeLib, $refGenomeFa, $jumpNGSQC, $jumpBWA, $jumpSortSam, $jumpMarkDuplicates, $jumpHaplotypeCaller) = @_;
    my $command = <<__COMMAND__;
#!/bin/bash

#BSUB -q $queue
#BSUB -n $thread
#BSUB -R "span[hosts=1] rusage[mem=${memory}GB]"
#BSUB -J step1.$germplasm
#BSUB -o step1.$germplasm.out
#BSUB -e step1.$germplasm.err

#create output directory
mkdir -p $outputDir

#run NGSQC
$jumpNGSQC mkdir -p $outputDir/NgsqcOut
$jumpNGSQC perl $softwareDir/NGSQCToolkit_v2.3.3/QC/IlluQC.pl -pe $germplasmDir/$germplasm/$germplasm-1.fq.gz $germplasmDir/$germplasm/$germplasm-2.fq.gz N A -l 70 -s 20 -o $outputDir/NgsqcOut/$germplasm

#run BWA
$jumpBWA mkdir -p $outputDir/BwaOut
$jumpBWA $softwareDir/bwa-0.7.13/bwa mem -t $thread -Y -a -M -R "\@RG\\tID:$germplasm\\tLB:$germplasm\\tSM:$germplasm\\tPL:ILLUMINA" $refGenomeLib $outputDir/NgsqcOut/$germplasm/$germplasm-1.fq.gz_filtered $outputDir/NgsqcOut/$germplasm/$germplasm-2.fq.gz_filtered 2> $outputDir/BwaOut/$germplasm.log | $softwareDir/samtools-0.1.18/samtools view -bS - > $outputDir/BwaOut/$germplasm.bam

#run SortSam
$jumpSortSam mkdir -p $outputDir/SortSamOut
$jumpSortSam $softwareDir/gatk-4.1.2.0/gatk --java-options "-Dsamjdk.compression_level=5 -Xms4000m" SortSam --INPUT $outputDir/BwaOut/$germplasm.bam --OUTPUT $outputDir/SortSamOut/$germplasm.sort.bam --SORT_ORDER coordinate --CREATE_INDEX false --CREATE_MD5_FILE false --VALIDATION_STRINGENCY SILENT

#run MarkDuplicates
$jumpMarkDuplicates mkdir -p $outputDir/MarkDuplicatesOut
$jumpMarkDuplicates $softwareDir/gatk-4.1.2.0/gatk --java-options "-Dsamjdk.compression_level=5 -Xms4000m" MarkDuplicates --INPUT $outputDir/SortSamOut/$germplasm.sort.bam --OUTPUT $outputDir/MarkDuplicatesOut/$germplasm.sort.mark.bam --METRICS_FILE $outputDir/MarkDuplicatesOut/$germplasm.metrics --VALIDATION_STRINGENCY SILENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --ASSUME_SORT_ORDER coordinate
$jumpMarkDuplicates $softwareDir/samtools-0.1.18/samtools index $outputDir/MarkDuplicatesOut/$germplasm.sort.mark.bam

#run HaplotypeCaller (non-BQSR)
$jumpHaplotypeCaller mkdir -p $outputDir/HaplotypeCallerOut_nonBQSR
$jumpHaplotypeCaller $softwareDir/gatk-4.1.2.0/gatk --java-options "-Xmx${memory}G" HaplotypeCaller -R $refGenomeFa -I $outputDir/MarkDuplicatesOut/$germplasm.sort.mark.bam -O $outputDir/HaplotypeCallerOut_nonBQSR/$germplasm.g.vcf -ERC GVCF -stand-call-conf 30 -mbq 20 --native-pair-hmm-threads 20

__COMMAND__
    print $handle $command;
}
