#!/usr/bin/perl -w
# ========================command line =======================
# ./samCheckInsertNM_paired.pl --input=align_filtered.sam --output=align_rmContam.sam --insertSize=100 --NMSize=5 --maxQ=50 --sumQ=90
use warnings;
use strict;
use Switch;
use Getopt::Long;
use IO::File;
use List::Util qw[min max];

## Global variables
my %options = ();

$options{'input'} = '';
$options{'output'} = '';
$options{'insertSize'} = 100;
$options{'nmSize'}  = 5;
$options{'maxQ'}  = 50;
$options{'sumQ'}  = 90;

my $optionsResult
	= GetOptions(\%options,
		"input=s",
		"output=s",
		"insertSize=i",
		"nmSize=i",
		"maxQ=i",
		"sumQ=i",
		);
		
my $input = $options{'input'};
my $output = $options{'output'};
# my $output = $fileName."_NMfiltered.sam";
my $insertSize = $options{'insertSize'};
my $nmSize = $options{'nmSize'};
my $maxQ = $options{'maxQ'};
my $sumQ = $options{'sumQ'};
my $index1 = rindex($input, '/');
my $fileName = substr($input, $index1 + 1, -4);
my $qHistOut = $fileName."_qScoreHist.txt";
my %qVal;
open(samIn, '<'.$input) || die "Err samIn:$!\n";
open(samOut, '>'.$output) || die "Err samOut:$!\n";
open(qHist, '>'.$qHistOut) || die "Err qHist:$!\n";

while (<samIn>) {
	if(/^@/) {
	  	print samOut $_;
    }
	else {
		# first line of the same pID
		my $in1 = $_;
		chomp $in1;
		my @line1 = split /\t/,$in1;
		# second line of the same pID
		my $in2 = <samIn>;
		last unless defined $in2;
		chomp $in2;
		my @line2 = split /\t/,$in2;		
 		$qVal{$line1[4]} ++;
		$qVal{$line2[4]} ++;
	    # if both qScore < maxQ, skip this pair
		if(max($line1[4], $line2[4]) < $maxQ) {
			next;
		} 
		# if sum of qScore >= sumQ, consider insertSize and NM
		elsif(($line1[4] + $line2[4]) >= $sumQ) {
			my $flag1 = 1;
			my $flag2 = 1;
			if(abs($line1[8]) <= $insertSize) {
				my $nm1 = $line1[11];
				my $nm2 = $line2[11];
				if ($nm1 !~ /^NM/) {
					$flag1 = 0;
				}
				else {
					$nm1 =~ s/NM\:i\://g;
					$flag1 = 0 if ($nm1 >= 5);
				}
				if ($nm2 !~ /^NM/) {
					$flag2 = 0;
				}
				else {
					$nm2 =~ s/NM\:i\://g;
					$flag2 = 0 if ($nm2 >= 5);
				}
			}
		print samOut $in1."\n" if($flag1);
		print samOut $in2."\n" if($flag2);
		}
		# one read qScore good, one read qScore very poor
		else {
			if($line1[4] > $line2[4]) {
				print samOut $in1."\n"
			} 
			else {
				print samOut $in2."\n"
			}
		}
	}
}
foreach my $qNum (sort { $a <=> $b } keys %qVal) {
    print qHist "$qNum\t$qVal{$qNum}\n";
}
