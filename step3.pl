#!/usr/bin/env perl

#date:   2025-07-22
#writer: Xiangjian Gou (QQ = 862137261)

#load modules
use strict;
use warnings;
use Getopt::Long;

#record version information
my $VERSION = 'variant (SNP/INDEL) calling, step3 v1.2 (2025-07-22)';

#set default options
my $germplasmDir = "/public/home/xjgou/rawReads";
my $refGenomeFa = "/public/home/xjgou/genome/genome.fa";
my $softwareDir = "/public/home/xjgou/tools";
my $outputDir = "output";
my $mbq = 20;
my $jumpBaseRecalibrator = '';
my $jumpApplyBQSR = '';
my $jumpHaplotypeCaller = '';
my $queue = "normal";
my $thread = 6;
my $memory = 30;
my $version;
my $help;

#get options from command line
GetOptions(
    'germplasmDir=s'           => \$germplasmDir,
    'refGenomeFa=s'            => \$refGenomeFa,
    'softwareDir=s'            => \$softwareDir,
    'outputDir=s'              => \$outputDir,
    'mbq|b=i'                  => \$mbq,
    'jumpBaseRecalibrator|jb+' => \$jumpBaseRecalibrator,
    'jumpApplyBQSR|ja+'        => \$jumpApplyBQSR,
    'jumpHaplotypeCaller|jh+'  => \$jumpHaplotypeCaller,
    'queue=s'                  => \$queue,
    'thread=i'                 => \$thread,
    'memory|m=i'               => \$memory,
    'version+'                 => \$version,
    'help+'                    => \$help,
);

#describe program information
my $usage = <<__GUIDE__;
####################################################################################################
Function: variants calling for each sample, including: BaseRecalibrator, ApplyBQSR, and HaplotypeCaller

Usage: perl step3.pl option1 <value1> option2 <value2> ... optionN <valueN>

Options:
  #Options for path:
  -g  | -germplasmDir <STRING> : the directory of storing all germplasms. (default: /public/home/xjgou/rawReads)
  -r  | -refGenomeFa  <STRING> : the reference genomic fasta file. (default: /public/home/xjgou/genome/genome.fa)
  -s  | -softwareDir  <STRING> : the directory of storing all software. (default: /public/home/xjgou/tools)
  -o  | -outputDir    <STRING> : set a directory for storing output information. (default: output)
  -b  | -mbq          <INT>    : minimum base quality required to consider a base for calling (default: 20)

  #Options for step:
  -jb | -jumpBaseRecalibrator  : no execute BaseRecalibrator
  -ja | -jumpApplyBQSR         : no execute ApplyBQSR
  -jh | -jumpHaplotypeCaller   : no execute HaplotypeCaller

  #Options for resources:
  -q  | -queue  <STRING>       : set the queue to use. (default: normal)
  -t  | -thread    <INT>       : set the number of threads to use. (default: 6)
  -m  | -memory    <INT>       : set the size of memory to use. (default: 30 [30GB])

  #Options for other:
  -v  | -version               : show the version information.
  -h  | -help                  : show the help information.
####################################################################################################

__GUIDE__

#output version and help information
die "$VERSION\n" if $version;
die $usage if $help;

#set whether to execute each step
$jumpBaseRecalibrator = '#' if $jumpBaseRecalibrator;
$jumpApplyBQSR = '#' if $jumpApplyBQSR;
$jumpHaplotypeCaller = '#' if $jumpHaplotypeCaller;

#generate a scheduling script for each germplasm, meanwhile, generate a comprehensive scheduling script
open my $oTOTAL, '>', 'step3.sh';
my $lsfDir = "step3.lsf";
system "mkdir -p $lsfDir";
foreach my $dir (glob "$germplasmDir/*") {
    my ($germplasm) = $dir =~ /([^\/]+)\z/;
    open my $oEACH, '>', "$lsfDir/$germplasm.lsf";
    lsfInfo($oEACH, $germplasm, $queue, $thread, $memory, $outputDir, $softwareDir, $refGenomeFa, $jumpBaseRecalibrator, $jumpApplyBQSR, $jumpHaplotypeCaller, $mbq);
    close $oEACH;
    print $oTOTAL "bsub < $lsfDir/$germplasm.lsf\n";
}
close $oTOTAL;

#create a subroutine to write all command into lsf script
sub lsfInfo {
    my ($handle, $germplasm, $queue, $thread, $memory, $outputDir, $softwareDir, $refGenomeFa, $jumpBaseRecalibrator, $jumpApplyBQSR, $jumpHaplotypeCaller, $mbq) = @_;
    my $command = <<__COMMAND__;
#!/bin/bash

#BSUB -q $queue
#BSUB -n $thread
#BSUB -J step3.$germplasm
#BSUB -o step3.$germplasm.out
#BSUB -e step3.$germplasm.err

#run BaseRecalibrator
$jumpBaseRecalibrator mkdir -p $outputDir/BaseRecalibratorOut
$jumpBaseRecalibrator $softwareDir/gatk-4.1.2.0/gatk --java-options "-Xms4000m" BaseRecalibrator -R $refGenomeFa -I $outputDir/MarkDuplicatesOut/$germplasm.sort.mark.bam -O $outputDir/BaseRecalibratorOut/$germplasm.table --known-sites $outputDir/knownsites/SNP.vcf --known-sites $outputDir/knownsites/INDEL.vcf

#run ApplyBQSR
$jumpApplyBQSR mkdir -p $outputDir/ApplyBQSROut
#$jumpApplyBQSR $softwareDir/gatk-4.1.2.0/gatk --java-options "-Xms4000m" ApplyBQSR -R $refGenomeFa -I $outputDir/MarkDuplicatesOut/$germplasm.sort.mark.bam -O $outputDir/ApplyBQSROut/$germplasm.bqsr.bam -bqsr $outputDir/BaseRecalibratorOut/$germplasm.table --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 -add-output-sam-program-record
$jumpApplyBQSR $softwareDir/gatk-4.1.2.0/gatk --java-options "-Xms4000m" ApplyBQSR -R $refGenomeFa -I $outputDir/MarkDuplicatesOut/$germplasm.sort.mark.bam -O $outputDir/ApplyBQSROut/$germplasm.bqsr.bam -bqsr $outputDir/BaseRecalibratorOut/$germplasm.table --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 -add-output-sam-program-record --create-output-bam-index false
$jumpApplyBQSR $softwareDir/samtools-0.1.18/samtools index $outputDir/ApplyBQSROut/$germplasm.bqsr.bam

#run HaplotypeCaller
$jumpHaplotypeCaller mkdir -p $outputDir/HaplotypeCallerOut
$jumpHaplotypeCaller $softwareDir/gatk-4.1.2.0/gatk --java-options "-Xmx${memory}G" HaplotypeCaller -R $refGenomeFa -I $outputDir/ApplyBQSROut/$germplasm.bqsr.bam -O $outputDir/HaplotypeCallerOut/$germplasm.g.vcf -ERC GVCF -stand-call-conf 30 -mbq $mbq --native-pair-hmm-threads 25

__COMMAND__
    print $handle $command;
}
