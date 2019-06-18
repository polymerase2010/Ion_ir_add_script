#!/usr/bin/perl 
#add ICGC data
print "\nStart:add ICGC lung cancer frequency in CN,KR and US.\n";
use strict;
open (A,"/results/langerdan/ICGC/LUSC-CN/mut/mut_var_CN.xls") || die $!;
open (B,"/results/langerdan/ICGC/LUSC-KR/mut/mut_var_KR.xls") || die $!;
open (C,"/results/langerdan/ICGC/LUSC-US/mut/mut_var_US.xls") || die $!;
open (D,$ARGV[0]) || die $!;
open (E,">$ARGV[0].icgc")  || die $!;
my %hash;my %hash1;my %hash2;my %hash3;
while (<A>) {
	chomp;
	my @line = split(/\t/,$_);
	$hash1{$line[1]} = $line[5];
}
while (<B>) {
	chomp;
	my @line = split(/\t/,$_);
	$hash2{$line[1]} = $line[5];
}
while (<C>) {
	chomp;
	my @line = split(/\t/,$_);
	$hash3{$line[1]} = $line[5];
}
my $sift;
my $poly;
while (<D>) {
	chomp;
	if (/^#.*/){
	print E "$_\n" ;
	next;
	}
	if (/^Locus/){
		print E "$_\tMAF_ICGC_CN\tMAF_ICGC_KR\tMAF_ICGC_US\n";
		my @title=split /\t/;
		for(0..$#title){
			if($title[$_] eq "SIFT"){$sift=$_;}
			if($title[$_] eq "PolyPhen"){$poly=$_;}
		}
		next;
	}
	my @line = split(/\t/,$_);
	print "line$.\n";
	if ( $line[$sift]  && $line[$sift] >= 0 && $line[$sift] <= 0.05){
		$line[$sift]= "$line[$sift]:Deleterious";print "$line[$sift]\n";
	}elsif ($line[$sift] > 0.05) {
		$line[$sift]= "$line[$sift]:Benign";
	}
	if ( $line[$poly] && $line[$poly]>=0 && $line[$poly]<=0.15){
		$line[$poly]= "$line[$poly]:Benign";
	}elsif($line[$poly]>0.15) {
		$line[$poly]= "$line[$poly]:Deleterious";
	}
        if(!exists $hash1{$line[0]}){$hash1{$line[0]}='.';}
        if(!exists $hash2{$line[0]}){$hash2{$line[0]}='.';}
        if(!exists $hash3{$line[0]}){$hash3{$line[0]}='.';}
	print E join("\t",@line,$hash1{$line[0]},$hash2{$line[0]},$hash3{$line[0]}),"\n";
#	@{$hash{$line[0]}} = join ("\t",$line[0],$line[1],$line[2],$line[3],$line[4],$line[5],$line[6],$line[7],$line[8],$line[9],$line[10],$line[11],$line[12],$line[13],$line[14],$line[15],$line[16],$line[17],$line[$poly],$line[19],$line[20],$line[21],$line[22],$line[23],$line[24],$line[25],$line[26],$line[27],$line[28],$line[29],$line[30],$line[31],$line[32],$line[33],$line[34],$line[35],$line[36],$line[37],$line[38]);
		
}
#foreach my $k (keys %hash){
#	if(!exists $hash1{$k}){$hash1{$k}='.';}
#	if(!exists $hash2{$k}){$hash2{$k}='.';}
#	if(!exists $hash3{$k}){$hash3{$k}='.';}
#	print E "@{$hash{$k}}\t$hash1{$k}\t$hash2{$k}\t$hash3{$k}\n";
#}	

close A;
close B;
close C;
close D;
print "End:add ICGC lung cancer frequency in CN,KR and US.\n";
