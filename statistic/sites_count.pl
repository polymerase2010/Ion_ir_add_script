#!/usr/bin/perl

use warnings;
use strict;
my ($all,$mut,$nocall,$mut_hot,$mut_novel,$nomut,$rm,$rm_hot,$rm_novel,$last,$last_hot,$last_novel,$gene);
my @file=`ls  *fre*xls`;
open(STAT,">$ARGV[0]")||die;
for my $f(@file){
	chomp $f;
	my $c=(split(" ",`wc -l $f`))[0]-1;
	if($f=~/variants_fre_all\.xls/){
		$all=$c;
	}elsif($f=~/variants_fre_mut\.xls/){
		$mut=$c;
	}elsif($f=~/variants_fre_nocall\.xls/){
		$nocall=$c;
	}elsif($f=~/variants_fre_mut_hotspot\.xls/){
		$mut_hot=$c;
	}elsif($f=~/variants_fre_mut_novel\.xls/){
		$mut_novel=$c;
	}elsif($f=~/variants_fre_nomut\.xls/){
		$nomut=$c;
	}elsif($f=~/variants_fre_rmnocall\.xls/){
		$rm=$c;
	}elsif($f=~/variants_fre_rmnocall_hotspot\.xls/){
		$rm_hot=$c;
	}elsif($f=~/variants_fre_rmnocall_novel\.xls/){
		$rm_novel=$c;
	}elsif($f=~/variants_fre\.xls/){
		$last=$c;
	}elsif($f=~/variants_fre_hotspot\.xls/){
		$last_hot=$c;
	}elsif($f=~/variants_fre_novel\.xls/){
		$last_novel=$c;
	}elsif($f=~/gene_fre\.xls/){
		$gene=$c;
	}
}
print STAT "Filter\tMut\tMut_hotspot\tMut_novel\tNocall\tNomut\tAll\n";
print STAT "Raw:\t$mut\t$mut_hot\t$mut_novel\t$nocall\t$nomut\t$all\n";
print STAT "Rmnocall:\t$rm\t$rm_hot\t$rm_novel\t-\t-\t-\n";
print STAT "Rmkghigh:\t$last\t$last_hot\t$last_novel\t-\t-\t-\n";
print STAT "Gene:\t$gene\t-\t-\t-\t-\t-\t-\n";
close STAT;
