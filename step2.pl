#!/usr/bin/env perl

#date:   2021-11-15
#writer: Xiangjian Gou (QQ = 862137261)

#load modules
use strict;
use warnings;
use Cwd;
use Getopt::Long;

#record version information
my $VERSION = 'variant (SNP/INDEL) calling, step2 v1.0 (2021.11.15)';

#get current working directory
my $cwd = getcwd;

#set default options
my $germplasmDir = "/public/home/xjgou/base_composition/six_sample/rawReads";
my $refGenomeFa = "/public/home/xjgou/base_composition/six_sample/barleyv2/barleyv2.fa";
my $intersectionScript = "/public/home/xjgou/base_composition/six_sample/getIntersectionVariants.pl";
my $filterScript = "/public/home/xjgou/base_composition/six_sample/filterVcf.pl";
my $softwareDir = "/public/home/xjgou/tools";
my $outputDir = "output";
my $jumpSamtools = '';
my $jumpGatk = '';
my $jumpIntersection = '';
my $queue = "smp";
my $thread = 4;
my $memory = 100;
my $version;
my $help;

#get options from command line
GetOptions(
    'germplasmDir=s'        => \$germplasmDir,
    'refGenomeFa=s'         => \$refGenomeFa,
    'intersection=s'        => \$intersectionScript,
    'filterScript=s'        => \$filterScript,
    'softwareDir=s'         => \$softwareDir,
    'outputDir=s'           => \$outputDir,
    'jumpSamtools|js+'      => \$jumpSamtools,
    'jumpGatk|jg+'          => \$jumpGatk,
    'jumpIntersection|ji+'  => \$jumpIntersection,
    'queue=s'               => \$queue,
    'thread=i'              => \$thread,
    'memory=i'              => \$memory,
    'version+'              => \$version,
    'help+'                 => \$help,
);

#describe program information
my $usage = <<__GUIDE__;
####################################################################################################
Function: get the knownsites variants (SNPs+INDELs) for subsequent Base Quality Score Recalibrator

Usage: perl step2.pl option1 <value1> option2 <value2> ... optionN <valueN>

Options:
  #Options for path:
  -g  | -germplasmDir <STRING> : the directory of storing all germplasms. (default: /public/home/xjgou/base_composition/six_sample/rawReads)
  -r  | -refGenomeFa  <STRING> : the reference genomic fasta file. (default: /public/home/xjgou/base_composition/six_sample/barleyv2/barleyv2.fa)
  -i  | -intersection <STRING> : the script for intersection. (default: /public/home/xjgou/base_composition/six_sample/getIntersection.pl)
  -f  | -filterScript <STRING> : the script for filtering. (default: /public/home/xjgou/base_composition/six_sample/filterVcf.pl)
  -s  | -softwareDir  <STRING> : the directory of storing all software. (default: /public/home/xjgou/tools)
  -o  | -outputDir    <STRING> : set a directory for storing output information. (default: output)

  #Options for step:
  -js | -jumpSamtools          : no execute samtools
  -jg | -jumpGatk              : no execute GATK
  -ji | -jumpIntersection      : no execute intersection

  #Options for resources:
  -q  | -queue  <STRING>       : set the queue to use. (default: smp)
  -t  | -thread    <INT>       : set the number of threads to use. (default: 4)
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
$jumpSamtools = '#' if $jumpSamtools;
$jumpGatk = '#' if $jumpGatk;
$jumpIntersection = '#' if $jumpIntersection;

#get the all sequence names for subsequent GenomicsDBImport
open my $iFASTA, '<', $refGenomeFa or die "Error: cannot open file '$refGenomeFa': $!";
my @names;
while (<$iFASTA>) {
    if (/\A>/) {
        my ($name) = /\A>(\S+)/;
        push @names, $name;
    }
}
close $iFASTA;
my $chrSets = join " ", map {"-L $_"} @names;

#get the all sorted bam files and gvcf files for subsequent bcftools and GATK, respectively
my @sortedBamFiles;
my @gvcfFiles;
foreach my $dir (glob "$germplasmDir/*") {
    my ($germplasm) = $dir =~ /([^\/]+)\z/;
    push @sortedBamFiles, "$outputDir/SortSamOut/$germplasm.sort.bam";
    push @gvcfFiles, "$outputDir/HaplotypeCallerOut_nonBQSR/$germplasm.g.vcf";
}
my $sortedBamSet = join " ", @sortedBamFiles;
my $gvcfSets = join " ", map {"-V $_"} @gvcfFiles;

#generate a scheduling script for each germplasm, meanwhile, generate a comprehensive scheduling script
my $lsfDir = "step2.lsf";
system "mkdir -p $lsfDir";
open my $oTOTAL, '>', 'step2.sh';
print $oTOTAL "bsub < $lsfDir/1.lsf\n";
close $oTOTAL;
open my $oEACH, '>', "$lsfDir/1.lsf";
lsfInfo($oEACH, $queue, $thread, $memory, $outputDir, $softwareDir, $refGenomeFa, $sortedBamSet, $gvcfSets, $chrSets, $intersectionScript, $filterScript, $jumpSamtools, $jumpGatk, $jumpIntersection, $cwd);
close $oEACH;

#create a subroutine to write all command into lsf script
sub lsfInfo {
    my ($handle, $queue, $thread, $memory, $outputDir, $softwareDir, $refGenomeFa, $sortedBamSet, $gvcfSets, $chrSets, $intersectionScript, $filterScript, $jumpSamtools, $jumpGatk, $jumpIntersection, $cwd) = @_;
    my $command = <<__COMMAND__;
#!/bin/bash

#BSUB -q $queue
#BSUB -n $thread
#BSUB -R "span[hosts=1] rusage[mem=${memory}GB]"
#BSUB -J step2
#BSUB -o step2.out
#BSUB -e step2.err

#create output directory
mkdir -p $outputDir
mkdir -p $outputDir/knownsites

#samtools call
$jumpSamtools mkdir -p $outputDir/knownsites/samtools
$jumpSamtools $softwareDir/bcftools-1.9/bcftools mpileup -d 100000 -Ou -f $refGenomeFa $sortedBamSet --threads $thread | $softwareDir/bcftools-1.9/bcftools call -m -v -Ob -o $outputDir/knownsites/samtools/samtools.bcf
$jumpSamtools $softwareDir/bcftools-1.9/bcftools view $outputDir/knownsites/samtools/samtools.bcf > $outputDir/knownsites/samtools/samtools.vcf
$jumpSamtools perl $filterScript $outputDir/knownsites/samtools/samtools.vcf $outputDir/knownsites/samtools/samtools.refisnothet.vcf $outputDir/knownsites/samtools/samtools.refishet.vcf
$jumpSamtools $softwareDir/gatk-4.1.2.0/gatk IndexFeatureFile -F $outputDir/knownsites/samtools/samtools.refisnothet.vcf
$jumpSamtools $softwareDir/gatk-4.1.2.0/gatk VariantFiltration -V $outputDir/knownsites/samtools/samtools.refisnothet.vcf -O $outputDir/knownsites/samtools/samtools.filter.vcf --cluster 4 --window 10 --mask-extension 3 --filter-name "lowMQ" --filter "MQ < 40.0" --filter-name "lowDP" --filter "DP < 8.0" --filter-name "LowQual" --filter "QUAL < 20"
$jumpSamtools $softwareDir/bcftools-1.9/bcftools view -f PASS $outputDir/knownsites/samtools/samtools.filter.vcf > $outputDir/knownsites/samtools/tmp
$jumpSamtools mv $outputDir/knownsites/samtools/tmp $outputDir/knownsites/samtools/samtools.filter.vcf
$jumpSamtools rm $outputDir/knownsites/samtools/samtools.filter.vcf.idx
$jumpSamtools $softwareDir/gatk-4.1.2.0/gatk SelectVariants -O $outputDir/knownsites/samtools/samtools.SNP.vcf --variant $outputDir/knownsites/samtools/samtools.filter.vcf -select-type SNP
$jumpSamtools $softwareDir/gatk-4.1.2.0/gatk SelectVariants -O $outputDir/knownsites/samtools/samtools.INDEL.vcf --variant $outputDir/knownsites/samtools/samtools.filter.vcf -select-type INDEL

#GATK call
$jumpGatk mkdir -p $outputDir/knownsites/gatk
$jumpGatk $softwareDir/gatk-4.1.2.0/gatk --java-options "-Xmx${memory}g -Xms${memory}g" GenomicsDBImport --genomicsdb-workspace-path $outputDir/knownsites/gatk/db --batch-size 50 -R $refGenomeFa $gvcfSets $chrSets
$jumpGatk cd $outputDir/knownsites/gatk
$jumpGatk $softwareDir/gatk-4.1.2.0/gatk --java-options "-Xmx${memory}g -Xms${memory}g" GenotypeGVCFs -R $refGenomeFa -V gendb://db -O gatk.vcf -new-qual -G StandardAnnotation --use-new-qual-calculator
$jumpGatk cd $cwd
$jumpGatk $softwareDir/gatk-4.1.2.0/gatk VariantFiltration -V $outputDir/knownsites/gatk/gatk.vcf -O $outputDir/knownsites/gatk/gatk.HC.vcf --cluster 4 --window 10 --mask-extension 3 --filter-name "lowMQ" --filter "MQ < 40.0" --filter-name "lowDP" --filter "DP < 8.0" --filter-name "LowQual" --filter "QUAL < 20" --filter-name "lowQD" --filter "QD < 2.0" --filter-name "lowReadPosRankSum" --filter "ReadPosRankSum < -8.0" --filter-name "highFS" --filter "FS > 60.0" --filter-name "lowMQRankSum" --filter "MQRankSum < -12.5"
$jumpGatk $softwareDir/bcftools-1.9/bcftools view -f PASS $outputDir/knownsites/gatk/gatk.HC.vcf > $outputDir/knownsites/gatk/tmp
$jumpGatk mv $outputDir/knownsites/gatk/tmp $outputDir/knownsites/gatk/gatk.HC.vcf
$jumpGatk rm $outputDir/knownsites/gatk/gatk.HC.vcf.idx
$jumpGatk $softwareDir/gatk-4.1.2.0/gatk SelectVariants -O $outputDir/knownsites/gatk/gatk.SNP.vcf --variant $outputDir/knownsites/gatk/gatk.HC.vcf -select-type SNP
$jumpGatk $softwareDir/gatk-4.1.2.0/gatk SelectVariants -O $outputDir/knownsites/gatk/gatk.INDEL.vcf --variant $outputDir/knownsites/gatk/gatk.HC.vcf -select-type INDEL

#Keep the intersection of samtools and GATK
$jumpIntersection perl $intersectionScript $outputDir/knownsites/samtools $outputDir/knownsites/gatk $outputDir/knownsites
$jumpIntersection $softwareDir/gatk-4.1.2.0/gatk IndexFeatureFile -F $outputDir/knownsites/SNP.vcf
$jumpIntersection $softwareDir/gatk-4.1.2.0/gatk IndexFeatureFile -F $outputDir/knownsites/INDEL.vcf

__COMMAND__
    print $handle $command;
}
