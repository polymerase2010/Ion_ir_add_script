#!/usr/bin/perl 

#use warnings;
use strict;
print "\nStart:adjust the order of each column.\tadjust_order_v3.pl\n";
my @order=(0,1,6,11,12,13,14,15,33,101,102,103,104,105,106,42,49,41,27,23,24,25,26,28,43,44,45,16,17,18,19,34,35,52,36,37,38,39,40,20,21,22,29,30,31,32,2,3,4,5,7,8,9,10,46,47,48,50,51,53,54);
my @column=("Sample_name","Locus","Genes","Strand","Transcript","Coding","Amino Acid Change","Variant Effect","ClinVar","A","B","C","D","E","F","dbsnp_rs","seq","dbsnp_af","CHAF","MAF","EMAF","AMAF","GMAF","UCSC Common SNPs","MAF_ICGC_CN","MAF_ICGC_KR","MAF_ICGC_US","PhyloP","SIFT","Grantham","PolyPhen","Allele Coverage","Allele Ratio","Var_ref","p-value","Phred QUAL Score","Coverage","Ref+/Ref-/Var+/Var-","Homopolymer Length","PFAM","dbSNP","DGV","COSMIC","OMIM","Gene Ontology","DrugBank","Genotype","Ref","Type","No Call Reason","Location","Length","Variant ID","Variant Name","CDS","EXIN","cds_start_end","cds_seq","cds_length","Novel or Hotspot");
open(RAW,"$ARGV[0]")||die;
open(NEW,">$ARGV[0].adjust")||die;
while(<RAW>){
	chomp;
	if(/^#/){
		print NEW $_,"\n";	next;
	}
	my @line=split /\t/;
	print "line$.\n";
	$line[101]=' ';
	$line[102]=' ';
	$line[103]=' ';
	$line[104]=' ';
	$line[105]=' ';
	$line[106]=' ';
	if(/^Sample_name/){	
		$line[101]='A';
		$line[102]='B';
		$line[103]='C';
		$line[104]='D';
		$line[105]='E';
		$line[106]='F';
		for my $i(@order){
			print NEW $line[$i],"\t";
		}
		print NEW "\n";
#		print NEW join("\t",$line[0],$line[1],$line[6],$line[11],$line[12],$line[13],$line[14],$line[15],$line[33],'','',$line[42],$line[49],$line[41],$line[27],$line[23],$line[24],$line[25],$line[26],$line[28],$line[43],$line[44],$line[45],$line[16],$line[17],$line[18],$line[19],$line[34],$line[35],$line[52],$line[36],$line[37],$line[38],$line[39],$line[40],$line[20],$line[21],$line[22],$line[29],$line[30],$line[31],$line[32],$line[2],$line[3],$line[4],$line[5],$line[7],$line[8],$line[9],$line[10],$line[46],$line[47],$line[48],$line[50],$line[51]),"\n";
	}else{	
		for my $i(@order){
#			if($line[$i] eq ''){print "$i\n";}
			print NEW $line[$i],"\t";
		}
		print NEW "\n";
#		 print NEW join("\t",$line[0],$line[1],$line[6],$line[11],$line[12],$line[13],$line[14],$line[15],$line[33],'','',$line[42],$line[49],$line[41],$line[27],$line[23],$line[24],$line[25],$line[26],$line[28],$line[43],$line[44],$line[45],$line[16],$line[17],$line[18],$line[19],$line[34],$line[35],$line[52],$line[36],$line[37],$line[38],$line[39],$line[40],$line[20],$line[21],$line[22],$line[29],$line[30],$line[31],$line[32],$line[2],$line[3],$line[4],$line[5],$line[7],$line[8],$line[9],$line[10],$line[46],$line[47],$line[48],$line[50],$line[51]),"\n";
	}
}
close RAW;close NEW;
print "End:adjust the order of each column.\tadjust_order_v3.pl\n";

