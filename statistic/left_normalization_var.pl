#!/usr/bin/perl 
#left-normalization for variants
#USAGE:perl a.tsv
use strict;
my @file=`ls *highqual.txt`;
for my $f(@file){
	chomp $f;
	open(FILE,$f)||die;
	open(NEW,">$f.left_nor")||die;
	while(<FILE>){
		chomp;
		if(/^#|Chrom/){
			print NEW "$_\n";next;
		}
		my ($rightloc,$leftloc,$newref,$newalt)=(0,0,"","");
		my @line=split /\t/;
		print "old\told:\t$line[0]:$line[14]:$line[15]:$line[16]\n";
		for(my $i=0;$i<&min(length $line[15],length $line[16]);$i++){#remove the same part from right side.
			if(substr($line[15],((length $line[15])-1-$i),1) ne substr($line[16],((length $line[16])-1-$i),1)){
				$rightloc=$i;last;
			}
		}
		print "LAST:right:$rightloc\t";
		$newref=substr($line[15],0,((length $line[15])-$rightloc));
		$newalt=substr($line[16],0,((length $line[16])-$rightloc));
		print "right normalization:\t$line[0]:$line[14]:$newref:$newalt\n";
		for(my $j=0;$j<&min(length $newalt,length $newref);$j++){#remove the same part from left side.
#			print "J:$j\n";
			$leftloc=$j;
			if(substr($newalt,$j,1) ne substr($newref,$j,1)){
				$leftloc=$j-1;last;
                       	}
		}
		if(length $newalt eq length $newref){
			$leftloc+=1;
		}
		print "LAST:left:$leftloc\t";
		$newref=substr($newref,$leftloc);
		$newalt=substr($newalt,$leftloc);
		$line[14]=$line[14]+$leftloc;
		print "left normalization:\t$line[0]:$line[14]:$newref:$newalt\n";
		$line[15]=$newref;
		$line[16]=$newalt;
		print NEW join("\t",@line),"\n";
	}
	close FILE;close NEW;
}	
#	print "$leftloc,$rightloc,$newref,$newalt\n";


sub min(){
	my($a,$b)=@_;
	if($a>$b){
		return $b;
	}else{
		return $a;
	}
}
print "END!\n";
