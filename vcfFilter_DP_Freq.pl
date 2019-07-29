#!/usr/bin/perl -w
# ========================command line =======================
# ./vcfFilter_DP_Freq.pl --input=outVCF.vcf --list=sampleList.txt --DP=30 --AD=5 --FREQ=25
use strict;
use Switch;
use Getopt::Long;
use IO::File;
use File::Basename;

## Global variables
my %options = ();
$options{'input'} = '';
$options{'list'} = '';
$options{'DP'}    = 30;
$options{'AD'}    = 5;
$options{'FREQ'}    = 25;

my $optionsResult
	= GetOptions(\%options,
		"input=s",
		"list=s",
		"DP=i",
		"AD=i",
		"FREQ=i",
		);
		

my $in = $options{'input'};
my $fileName = $in;
$fileName =~ s/.vcf//g;
open(fileIn, "<$in") || die "Err fileIn:$!\n"; 
my @sampleList;
if(defined $options{'list'}) {
	open(fileList, "<".$options{'list'}) || die "fileList:$!\n"; 
	while(<fileList>) {
		chomp;
		my $name = basename($_);
		$name =~ s/\r//g;
		push @sampleList, $name;
	}
} 
my $DP = $options{'DP'};
my $AD = $options{'AD'};
my $FREQ = $options{'FREQ'};

open(fileOut, ">".$fileName."_final.vcf") || die "Err fileOut:$!\n";
open(fileDP, ">".$fileName."_DP.txt") || die "Err fileDP:$!\n";
open(fileFREQ, ">".$fileName."_FREQ.txt") || die "Err fileFREQ:$!\n";
# Print header
while (<fileIn>) {
	if(/^##/) {
		print fileOut $_;
	} elsif(/^#CHROM/) {
		if (scalar @sampleList == 0) {
			print fileOut "\n";
		} else {
			print fileOut "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
			print fileDP "Position";
			print fileFREQ "Position";
			foreach (@sampleList) {
				chomp;
				my $tmpName = $_;
				$tmpName =~ s/_final.bam//g;
				$tmpName =~ s/.vcf.gz//g;
				$tmpName =~ s/.vcf//g;
				print fileOut "\t$tmpName";
				print fileDP "\t$tmpName";
				print fileFREQ "\t$tmpName";
			}
			print fileOut "\n";
			print fileDP "\n";
			print fileFREQ "\n";
		}
	}	else {
		my $line = $_;
		chomp;
		my @input = split /\t/, $_;
		my $outString = join("\t", @input[0 .. 8]);
		my $outStringDP = $input[0].':'.$input[1].$input[3].'>'.$input[4];
		my $outStringFREQ = $input[0].':'.$input[1].$input[3].'>'.$input[4];
		my $outFlag = 0;
		foreach my $i (9 .. $#input) {
			if($input[$i] !~ /^\./) {
				my @score = split /:/, $input[$i];
				my $s1 = $score[3];
				my $s2 = $score[5];
				my $s3 = $score[6];
				$s3 =~ s/%//g;
				if ($s1 >= $DP && $s2 >= $AD && $s3 >= $FREQ) {
					$outFlag = 1;
				} elsif( $s3 < $FREQ) {
					substr($input[$i], 2, 1) = "0";
				}
				$outString = join("\t",  $outString, $input[$i]);
				$outStringDP = join("\t",  $outStringDP, $s1);
				$outStringFREQ = join("\t",  $outStringFREQ, $s3);
			} else {
				$outString = join("\t",  $outString, '.');
				$outStringDP = join("\t",  $outStringDP, 0);
				$outStringFREQ = join("\t",  $outStringFREQ, 0);
			}
		}
		if($outFlag == 1) {
			print fileOut $outString."\n";
			print fileDP $outStringDP."\n";
			print fileFREQ $outStringFREQ."\n";
		}	
	}
}
