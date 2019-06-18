#!/usr/bin/perl 
#left-normolization for variants
#USAGE:perl a.tsv
print "\nStart:left-normolization for variants of vcf.\tlocus_modify_v3.pl\n";
use strict;
open(IR,$ARGV[0])||die;
open(NEW,">$ARGV[0].alt")||die;
while(<IR>){
	chomp;
	if(/^#|Locus/){
		print NEW "$_\n";next;
	}
	$_=~s///g;
	my ($rightloc,$leftloc,$newref,$newalt)=(0,0,"","");
	my @line=split /\t/;
	print "line$.\n";
	if($#line<38){	@line=(@line,'.','.');	}
	my @type=split(/\//,$line[1]);#$type[1] eq $alt, $line[2] eq $ref
	if(($line[2] ne $type[0]) && ($type[0] ne $type[1])){#liangzhongtubian
		print NEW join("\t",@line),"\n";next;
        }
	if($line[3] eq 'SNV' ){ #print "$line[3]\n";
		if(length $line[2]==1){print NEW "$_\n";next;}
		for(my $i=0;$i<=length $line[2];$i++){
			if(substr($type[1],$i,1) ne substr($line[2],$i,1)){
				$leftloc=$i;
				$newref=substr($line[2],$i,1);
				$newalt=substr($type[1],$i,1);
				my @chr=split(/:/,$line[0]);
				$line[0]=$chr[0].":".($chr[1]+$i);
				last;
			}
		}
	}elsif($line[3] eq "INDEL"){#print "$line[3]\n";
#		if($line[2]>$type[1]){#del
			for(my $i=0;$i<&min(length $type[1],length $line[2]);$i++){#remove the same part from right side.
#				print "$line[1]\ti:$i\n";
				if(substr($type[1],((length $type[1])-1-$i),1) ne substr($line[2],((length $line[2])-1-$i),1)){
					$rightloc=$i;last;
				}
			}
			print "$line[1]:$line[2]\t";
			print "LAST:right:$rightloc\t";
			$newref=substr($line[2],0,((length $line[2])-$rightloc));
			$newalt=substr($type[1],0,((length $type[1])-$rightloc));
			print "right normolization:$newref/$newalt\n";
			for(my $j=0;$j<&min(length $newalt,length $newref);$j++){#remove the same part from left side.
#				print "J:$j\n";
				$leftloc=$j;
				if(substr($newalt,$j,1) ne substr($newref,$j,1)){
	                                $leftloc=$j-1;last;
                       		 }
			}
			print "LAST:left:$leftloc\t";
			$newref=substr($newref,$leftloc);
			$newalt=substr($newalt,$leftloc);
			print "left normolization:$newref/$newalt\n";
			my @chr=split(/:/,$line[0]);
			$line[0]=$chr[0].":".($chr[1]+$leftloc);
#                }
	}else{	print "ELSE TYPE=$line[3]\n";
		print NEW "$_\n";print "$line[1]\n";next;
	}		
#	print "$leftloc,$rightloc,$newref,$newalt\n";
	if($type[0] eq $type[1]){#homozytous
		$line[1]="$newalt/$newalt";
	}else{#heterozygous
		$line[1]="$newref/$newalt";
	}
	$line[2]=$newref;
	print NEW join("\t",@line),"\n";
}


sub min(){
	my($a,$b)=@_;
	if($a>$b){
		return $b;
	}else{
		return $a;
	}
}
print "End:left-normolization for variants of vcf.\tlocus_modify_v3.pl\n";
