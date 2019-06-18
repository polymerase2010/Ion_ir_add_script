#!/user/bin/perl -w
#mark a mutain as novel or hotspot
print "\nStart:Mark a mutation as novel or hotspot.\tplus_novel_v3.pl\n";
use strict;
open (A,"$ARGV[0]")||die;
open (B,">$ARGV[0].novel")||die;
while (<A>){
        chomp;
        if(/^##.*/){
 	       print B "$_\n";next;
        }
        if(/^Sample_name/){
        	print B "$_\t"."Novel or Hotspot"."\n";
        	next;
        }
        my @line=split /\t/;
	print "line$.\n";
	if($line[9] ne ''){
		print "$line[0]\t$line[1]\t$line[2]\t$line[12]\tcosmic=$line[9]\tHotspot\n";
		print B join("\t",@line,'Hotspot'),"\n";
	}else{
		print B join("\t",@line,'novel'),"\n";
	}
}
close A;
close B;
print "End:Mark a mutation as novel or hotspot.\tplus_novel_v3.pl\n";
