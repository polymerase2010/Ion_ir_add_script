#!/usr/bin/perl
#筛除千人基因组中国汉族正常人高频位点
print "\nStart:remove the sites with high frequency(more than 0.05) in CHAF.(this script just adds CHAF.)\tfilter_vcf_kg_base_v3.pl\n";
use strict;
my %kg;
open(KG,"/results/duzhao_test/1000GenomesProject/new/1000g_CH_AF_chp2_and_colon_lung_and_BRCA1_2.txt")||die;
while(<KG>){
	chomp;
	my @line=split /\t/;
	$kg{$line[0]}=$line[1];
}
close KG;
$kg{"chr5:149433596:TG:GA"}=0.7;#add sites not in 1000g_CH_AF_chp2_and_colon_lung.txt
print "Hash of Chinese allele frequency in 1000Genomes Project has completed.\n";#千人基因组中国汉人位点及频率存入%kg

open(VCF,"$ARGV[0]")||die;
open(NEW,">$ARGV[0].kg")||die;
my $intag=24;#插入位置信息
while(<VCF>){
	chomp;
	if(/^#/){print NEW $_."\n";next;}
	my @line=split /\t/;
	print "line$.\n";
	my $mutation;
	if($ARGV[1] eq "IR"){
		if(/^Locus/){
			if(/No Call Reason/){$intag=25;}
			print NEW join("\t",@line[0..$intag],"CHAF",@line[($intag+1)..$#line]),"\n";
			next;
		}
		$mutation=$line[0].":".$line[2].":".(split(/\//,$line[1]))[1];
#		print "Mut $mutation\n";
	}elsif($ARGV[1] eq "Vanno"){
                $mutation="chr".$line[2].":".$line[3].":".$line[4].":".$line[5];
	}elsif($ARGV[1] eq "xls"){
		next unless($line[4] eq "Homozygous"||$line[4] eq "Heterozygous");
		$mutation=$line[0].":".$line[1].":".$line[2].":".$line[3];
	}else{
		print "Need parameter 'IR' or 'Vanno' or 'xls'!\n";
		last;
	}
	if(exists $kg{$mutation}){
		if($kg{$mutation}<=1){
			if($ARGV[1] eq "IR"){
				print NEW join("\t",@line[0..$intag],$kg{$mutation},@line[($intag+1)..$#line]),"\n";
			}else{
                                print NEW join("\t",@line,$kg{$mutation}),"\n";
			}
		}else{
			print $mutation."\t".$kg{$mutation}."\n";
		}
	}else{
		print NEW join("\t",@line[0..$intag],".",@line[($intag+1)..$#line]),"\n";
	}
}
close VCF;
close NEW;
print "End:remove the sites with high frequency(more than 0.05) in CHAF.(this script just adds CHAF.)\tfilter_vcf_kg_base_v3.pl\n";
