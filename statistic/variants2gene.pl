#!/usr/bin/perl -w

use strict;
open(FILE,"variants_rate.xls")||die;
my %gene;
while(<FILE>){
	chomp;
	my @line=split(/\t/,$_);
	push @{$gene{$line[1]}}, $line

       for(my $i=5; $i<scalar @line; $i++){
      push @{$gene{$1}},$line[$i];
    }
  }
}
foreach (keys %gene){
  print $_,"\t",join(",",@{$gene{$_}}),"\t",scalar @{$gene{$_}},"\n";
}
