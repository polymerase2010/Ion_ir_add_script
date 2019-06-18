#!/usr/bin/perl -w
#add the sequence at both sides of the mutation
print "\nStart:add the sequence at both sides of the mutation.\tadd_bothside_sequence_v3.pl\n";
use strict;
#提取人类基因组序列
open(HG,"/results/referenceLibrary/tmap-f3/hg19/hg19.fasta")||die;
my %hg;
my $seq="";
my $header;
while(<HG>){
	chomp;
	my $line=$_;
	if(/^>/){
		if(length $seq>0){
			$hg{$header}=$seq;
		}
		$seq="";
		$header=substr($line,1);#去掉>符号
#		print "$header\n";
	}else{
		$seq=$seq.$line;
	}
}
$hg{$header}=$seq;
print "Hash of hg19 sequence ready.\n";
close HG;
open(BED,"$ARGV[0]")||die;
open(NEW,">$ARGV[0].seq")||die;
#print NEW "Chrom\tPosition\tRef\tVariant\tAllele Call\tFrequency\tQuality\tType\tGene ID\tRegion Name\tCoverage\tStrand Bias\tSequence\n";
while(<BED>){
	chomp;
	if(/^#/){print NEW "$_\n";next;}
	if(/^Locus/){
		print NEW "$_\tseq\tcds_seq\n";
		next;
	}
	my @line=split /\t/;
	print "line$.\n";
	my @loc=split(/:/,$line[0]);
	my $part;
#	my $part=substr($hg{$line[0]},$line[1]-200,350);
#	print NEW join("\t",@line,$part,substr($part,199,length $line[2])),"\n";
	if($line[12]=~/del/){
		$part=substr($hg{$loc[0]},$loc[1]-9,8).lc(substr($line[2],1)).substr($hg{$loc[0]},$loc[1]+$line[7]-1,8);
	}elsif($line[12]=~/ins/){
		$part=substr($hg{$loc[0]},$loc[1]-10,8).lc($line[2]).substr($hg{$loc[0]},$loc[1]+$line[7]-1,8);
	}else{
		$part=substr($hg{$loc[0]},$loc[1]-9,8).lc($line[2]).substr($hg{$loc[0]},$loc[1]+$line[7]-1,8);
        }
	my $part2=substr($hg{$loc[0]},$loc[1]-9,17);
	my $part3;#cds sequence
	my @cds=split(/ /,$line[47]);
#	print "cds_start_end=@cds\n";
	if($line[10] eq '+' || $line[10] eq '+, -'){
		for(my $i=0;$i<scalar@cds;$i+=2){
			$part3.=substr($hg{$loc[0]},$cds[$i]-1,$cds[$i+1]-$cds[$i]+1);
		}
		print "sense strand lentgh=",length $part3,"\n";
	}elsif($line[10] eq '-'){
		for(my $i=scalar@cds-1;$i>0;$i-=2){
                        $part3.=substr($hg{$loc[0]},$cds[$i]-1,$cds[$i-1]-$cds[$i]+1);
                }
		print "antisense strand length=",length $part3,"\n";
		$part3=&seq_reverse($part3);
	}
#	print "cds_seq=$part3\n";
#	print "$part\n$part2\n\n";
	if($line[10] eq '-'){
	}
	print NEW join("\t",@line,$part,$part3),"\n";
}
close BED;
close NEW;

sub seq_reverse(){
	my ($seq)=@_;
	my $reverse;
	my %rever=('A'=>'T','T'=>'A','G'=>'C','C'=>'G','a'=>'t','t'=>'a','g'=>'c','c'=>'g');
	for(my $i=(length $seq)-1;$i>=0;$i--){
		$reverse.=$rever{(substr($seq,$i,1))};
	}
#	print "re:$seq\n$reverse\n";
	return $reverse;
}
print "End::add the sequence at both sides of the mutation.\tadd_bothside_sequence_v3.pl\n";
