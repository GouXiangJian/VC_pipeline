#!/usr/bin/env perl

#date:   2025-07-22
#writer: Xiangjian Gou (QQ = 862137261)

#load modules
use strict;
use warnings;
use Getopt::Long;

#record version information
my $VERSION = 'variant (SNP/INDEL) calling, step5 v1.2 (2025-07-22)';

#set default options
my $population = "population_genotype";
my $filterScript = "/public/home/xjgou/snpFilter.pl";
my $filterMAF = 0.05;
my $filterMR = 0.2;
my $refGenomeFa = "/public/home/xjgou/genome/genome.fa";
my $softwareDir = "/public/home/xjgou/tools";
my $outputDir = "output";
my $queue = "smp";
my $thread = 10;
my $memory = 100;
my $version;
my $help;

#get options from command line
GetOptions(
    'population=s'       => \$population,
    'filterScript|fs=s'  => \$filterScript,
    'filterMAF|ff=f'     => \$filterMAF,
    'filterMR|fr=f'      => \$filterMR,
    'refGenomeFa=s'      => \$refGenomeFa,
    'softwareDir=s'      => \$softwareDir,
    'outputDir=s'        => \$outputDir,
    'queue=s'            => \$queue,
    'thread=i'           => \$thread,
    'memory=i'           => \$memory,
    'version+'           => \$version,
    'help+'              => \$help,
);

#describe program information
my $usage = <<__GUIDE__;
####################################################################################################
Function: variant filtering, including: GatherVcfs, VariantFiltration, SelectVariants, and Tassel5

Usage: perl step5.pl option1 <value1> option2 <value2> ... optionN <valueN>

Options:
  #Options for filtering:
  -p  | -population   <STRING> : the population name. (default: population_genotype)
  -fs | -filterScript <STRING> : the filtering script. (default: /public/home/xjgou/snpFilter.pl)
  -ff | -filterMAF    <FLOAT>  : the minor allele frequency. (default: 0.05)
  -fr | -filterMR     <FLOAT>  : the missing rate. (default: 0.2)

  #Options for path:
  -r  | -refGenomeFa  <STRING> : the reference genomic fasta file. (default: /public/home/xjgou/genome/genome.fa)
  -s  | -softwareDir  <STRING> : the directory of storing all software. (default: /public/home/xjgou/tools)
  -o  | -outputDir    <STRING> : set a directory for storing output information. (default: output)

  #Options for resources:
  -q  | -queue        <STRING> : set the queue to use. (default: smp)
  -t  | -thread       <INT>    : set the number of threads to use. (default: 10)
  -m  | -memory       <INT>    : set the size of memory to use. (default: 100 [100GB])

  #Options for other:
  -v  | -version               : show the version information.
  -h  | -help                  : show the help information.
####################################################################################################

__GUIDE__

#output version and help information
die "$VERSION\n" if $version;
die $usage if $help;

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
my $input = join ' ', map {"-I $outputDir/GenotypeGVCFsOut/$_.vcf"} @names;

#generate a scheduling script for filtering, meanwhile, generate a comprehensive scheduling script
my $lsfDir = "step5.lsf";
system "mkdir -p $lsfDir";
open my $oTOTAL, '>', 'step5.sh';
print $oTOTAL "bsub < $lsfDir/1.lsf\n";
close $oTOTAL;
open my $oEACH, '>', "$lsfDir/1.lsf";
lsfInfo($oEACH, $queue, $thread, $memory, $input, $outputDir, $softwareDir, $population, $filterScript, $filterMAF, $filterMR);
close $oEACH;

#create a subroutine to write all command into lsf script
sub lsfInfo {
    my ($handle, $queue, $thread, $memory, $input, $outputDir, $softwareDir, $population, $filterScript, $filterMAF, $filterMR) = @_;
    my $command = <<__COMMAND__;
#!/bin/bash

#BSUB -q $queue
#BSUB -n $thread
#BSUB -J step5
#BSUB -o step5.out
#BSUB -e step5.err

#run GatherVcfs
$softwareDir/gatk-4.1.2.0/gatk GatherVcfs $input -O $outputDir/GenotypeGVCFsOut/$population.vcf

#run VariantFiltration
mkdir -p $outputDir/final
$softwareDir/gatk-4.1.2.0/gatk VariantFiltration -V $outputDir/GenotypeGVCFsOut/$population.vcf -O $outputDir/final/$population.HC.vcf --cluster 4 --window 10 --mask-extension 3 --filter-name "lowMQ" --filter "MQ < 40.0" --filter-name "lowDP" --filter "DP < 8.0" --filter-name "LowQual" --filter "QUAL < 20" --filter-name "lowQD" --filter "QD < 2.0" --filter-name "lowReadPosRankSum" --filter "ReadPosRankSum < -8.0" --filter-name "highFS" --filter "FS > 60.0" --filter-name "lowMQRankSum" --filter "MQRankSum < -12.5"

#run bcftools
$softwareDir/bcftools-1.9/bcftools view -f PASS $outputDir/final/$population.HC.vcf > $outputDir/final/tmp
mv $outputDir/final/tmp $outputDir/final/$population.HC.vcf
rm $outputDir/final/$population.HC.vcf.idx

#run SelectVariants
$softwareDir/gatk-4.1.2.0/gatk SelectVariants -O $outputDir/final/$population.SNP.vcf --variant $outputDir/final/$population.HC.vcf -select-type SNP
$softwareDir/gatk-4.1.2.0/gatk SelectVariants -O $outputDir/final/$population.INDEL.vcf --variant $outputDir/final/$population.HC.vcf -select-type INDEL

#run SortGenotypeFilePlugin
perl $softwareDir/tassel5/run_pipeline.pl -Xmx${memory}g -SortGenotypeFilePlugin -inputFile $outputDir/final/$population.SNP.vcf -outputFile $outputDir/final/$population.SNP.sort.vcf -fileType VCF
perl $softwareDir/tassel5/run_pipeline.pl -Xmx${memory}g -SortGenotypeFilePlugin -inputFile $outputDir/final/$population.INDEL.vcf -outputFile $outputDir/final/$population.INDEL.sort.vcf -fileType VCF

#get final hmphap file
perl $softwareDir/tassel5/run_pipeline.pl -Xmx${memory}g -fork1 -vcf $outputDir/final/$population.SNP.sort.vcf -export $outputDir/final/$population.SNP -exportType Hapmap -runfork1
perl $softwareDir/tassel5/run_pipeline.pl -Xmx${memory}g -fork1 -vcf $outputDir/final/$population.INDEL.sort.vcf -export $outputDir/final/$population.INDEL -exportType Hapmap -runfork1

#filter snp
perl $filterScript $outputDir/final/$population.SNP.hmp.txt $filterMAF $filterMR $outputDir/final/$population.SNP

__COMMAND__
    print $handle $command;
}
