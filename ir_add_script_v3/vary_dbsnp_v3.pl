#!/usr/bin/perl 
#add dbsnp 
print "\nStart:add dbsnp.\tplus_dbsnp_af_v3.pl\n";
use strict;
open(A,$ARGV[0]) || die;
open(B,"/results/duzhao_test/scripts/dbsnp/dbsnp_chp2_colonlung_BRCA.vcf") || die;
open(NEW,">$ARGV[0].dbsnp")  || die;
my %hash; my %hash1;my %hash2;
while(<B>){
	chomp;
	next if /^#.*/;
	my @line = split/\t/;
	my @alt=split(",",$line[4]);
	print "alter:@alt\n";
	my $loc;
	foreach my $alter(@alt){
		$loc = "chr$line[0]:$line[1]:$line[3]:$alter";
		$hash2{$loc} = $line[2]." ".$hash2{$loc};
	}
	if (/CAF=(.*);/){
		my @fre=split(",",$1);
		print "frequency:@fre\n";
		for(my $i=1;$i<scalar @fre;$i++){	
			$loc = "chr$line[0]:$line[1]:$line[3]:$alt[$i-1]";
			$hash1{$loc} = $fre[$i]." ".$hash1{$loc};
			print "$loc\t$hash2{$loc}\t$hash1{$loc}\n";
		}
	}
}
while (<A>){
        chomp;
        if(/^Sample/){
 	       print NEW "$_\n";next;
        }
        if(/^Locus/){
        	print NEW "$_\t"."dbsnp_af"."\t"."dbsnp_rs"."\n";
        	next;
        }
#	print "line $.\n";
        my @line=split/\t/;
	my $ref=$line[47];
	my $alt=(split("/",$line[46]))[1];
	if(not exists $hash1{"$line[1]:$ref:$alt"}){
		$hash1{"$line[1]:$ref:$alt"}='.';
	}
        if(not exists $hash2{"$line[1]:$ref:$alt"}){
                $hash2{"$line[1]:$ref:$alt"}='.';
        }
	print "$line[1]:$ref:$alt\n";
	$line[15]=$hash2{"$line[1]:$ref:$alt"};
	$line[17]=$hash1{"$line[1]:$ref:$alt"};
	print NEW join("\t",@line),"\n";
}
close B;
close A;
print "End:add dbsnp.\tplus_dbsnp_af_v3.pl\n";
