#!/usr/bin/perl
#AF=10, Q=200
use warnings;
use strict;
my @file=`ls *.xls`;
foreach my $f(@file){
	chomp $f;
	open(FILE,$f)||die;
	open(LOW,">${f}_lowqual.txt")||die;
	open(HIGH,">${f}_highqual.txt")||die;
	while(<FILE>){
		chomp;
		if(/^Chrom/){
			print LOW $_."\n";
			print HIGH $_."\n";
			next;
		}
		my @line=split /\t/;
		if($line[4]=~/Homo|Heter/){
			if(($line[6]>=10) or ($line[7]>=200)){
				print HIGH join("\t",@line),"\n";
			}else{
				print LOW join("\t",@line),"\n";
			}
		}else{
			print HIGH join("\t",@line),"\n";
		}
	}
	close FILE;close LOW;close HIGH;
}
