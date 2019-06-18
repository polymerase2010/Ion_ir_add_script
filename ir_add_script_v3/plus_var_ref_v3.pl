#!/usr/bin/perl -w
use strict;
print "\nStart:add variant allele frequency.\tplus_var_ref_v3.pl\n";
open (A,"$ARGV[0]")||die;
open (B,">$ARGV[0].var_ref")||die;
my $ref_af;
my $var_af;
while (<A>){
	chomp;
	my @line=split/\t/;
	print "line$.\n";
	my $temp=(split(/,/,$line[35]))[0];
#	print "$temp\n";
	if ($temp=~/(\d.*)/){
		$ref_af =$1;
		$var_af = 1-$1;
	#	print "$ref_af\t$var_af\n";
	}
	if (/^Sample_name/) {
		print B "$_\t"."Var_af\n";
		next;
	}
	else {
		print B "$_\t"."$var_af\n";	
	}
}
close A;
close B;
print "End:add variant allele frequency.\tplus_var_ref_v3.pl\n";
