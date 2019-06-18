#!/usr/bin/perl -w

use strict;
if(scalar @ARGV!=3){
	print "USAGE:perl gene_fre_v4.pl variants_fre.xls gene_fre.xls \$n(sample number)\n";
	exit;
}
open(VAR,"$ARGV[0]")||die;
open(GENE,">$ARGV[1]")||die;
my %gene;
while(<VAR>){
	chomp;
	my @line=split(/\t/,$_);
	my @sample=split(/\s/,$line[2]);
	foreach(@sample){
		if($_=~/H.*:(.*)/){
			push @{$gene{$line[1]}},$1;
		}
	}
}
print GENE "gene\tsample\tsam_num\tfrequency\n";
foreach my $k(sort (keys %gene)){
	my %count;
	@{$gene{$k}}=grep {++$count{$_}<2;} @{$gene{$k}};
	my $gene_fre=scalar @{$gene{$k}}/$ARGV[2];
	print GENE "$k\t@{$gene{$k}}\t",scalar @{$gene{$k}},"\t$gene_fre\n";
}
			
