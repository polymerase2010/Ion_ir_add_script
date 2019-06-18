#!/usr/bin/perl -w

#去除检测率（非./.比例）小于20%的位点
use strict;
my (%var,%var_gene,%novel_hot);
#open(SUM,"Summary_of_identified_variants.txt")||die;
my @file=`ls *left_nor`;
foreach my $f(@file){
	chomp $f;
	open(FILE,$f)||die;
	my %uni_var;
	while(<FILE>){
		chomp;
		if(/^Chrom/){ 	next;	}
		my @line=split("\t",$_);
		$line[4]=~s/\s/_/g;
		if(exists $var{"$line[0]:$line[14]:$line[15]:$line[16]"}{$f} && $var{"$line[0]:$line[14]:$line[15]:$line[16]"}{$f}=~/Homozygous|Heterzygous/){
			next;
		}else{	
			$var{"$line[0]:$line[14]:$line[15]:$line[16]"}{$f}=$line[4];#将每个样本的突变放入hash数组里
			$novel_hot{"$line[0]:$line[14]:$line[15]:$line[16]"}.="$line[10]:$f|";
   			$var_gene{"$line[0]:$line[14]:$line[15]:$line[16]"}=$line[12];
		}
	}
	close FILE;
}
my %kg;
open(KG,"/results/duzhao_test/1000GenomesProject/new/1000g_CH_AF_chp2_and_colon_lung_and_BRCA1_2.txt")||die;
while(<KG>){
	chomp;
	next if /CH_AF/;
	my @line=split /\t/;
	my @loc=split(/:/,$line[0]);
	$kg{$line[0]}=$line[1];
}
$kg{"chr5:149433596:TG:GA"}=0.5;
$kg{"chr4:55141050:AGCCCAGATGGACATG:AGCCCGGATGGACATG"}=1.0000;
$kg{"chr5:112175769:CGG:CAG"}=0.7957;
$kg{"chr17:7579470:CGG:CGC"}=0.5745;
#$kg{"chr4:55962545"}=0.4784;
close KG;
print "Chines han allele frequency in KG hash completed.\n";
open(VAR,">variants_fre_all.xls")||die;
print VAR "variants\tgene\tsample_detected\tdetection\tdetection_rate\tvar\tvar_fre\tkg_CHAF\tHot_or_novel\n";
foreach my $site(sort keys %var){
	my @sample;
	my $sam_count=keys %{$var{$site}};#sample count of this sites
	my ($homo,$heter,$no_detect,$absent)=(0,0,0,0);
	foreach my $s(keys %{$var{$site}}){
		if($var{$site}{$s}=~/Homozygous/){
			$homo++;
		}elsif($var{$site}{$s}=~/Heterozygous/){
			$heter++;
		}elsif($var{$site}{$s}=~/No_Call/){
			$no_detect++;
		}elsif($var{$site}{$s}=~/Absent/){
			$absent++;
		}
		push @sample,$var{$site}{$s}.":".$s;
	}
	my $var_num=$homo*2+$heter;
	my $af;
	if($homo+$heter+$absent==0){
		$af=100;
	}else{
		$af=$var_num/(2*($homo+$heter+$absent));#MAF计算，剔除未检测位点。$n是样本总数，$no_dec是没有检测到该位点的样本
	}
	my $detect_fre=($homo+$heter+$absent)/$sam_count;
	print VAR $site."\t".$var_gene{$site}."\t"."@sample"."\t".($homo+$heter+$absent)."|".$sam_count."\t".$detect_fre."\t".$var_num."|".(2*($homo+$heter+$absent))."\t";
	printf VAR "%0.4f",$af;
	if(exists $kg{$site}){#添加CHAF数据
		print VAR "\t$kg{$site}\t";
	}else{
		print VAR "\tNA\t";
	}
	if($novel_hot{"$site"}=~/Hotspot/){
		print VAR "Hotspot\n";
	}else{
		print VAR "Novel\n";
	}
}
close VAR;
