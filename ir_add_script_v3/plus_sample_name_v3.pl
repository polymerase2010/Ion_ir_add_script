#!/user/bin/perl -w
#add th sample name at the first column.
print "\nStart:add the sample name at the first column.\tplus_sample_name_v3.pl\n";
use strict;
open (IR,$ARGV[0])||die;
open (NEW,">$ARGV[0].name")||die;
my @sam=split(/\./,$ARGV[0]);
#print "$ARGV[0]\t$sam[0]\n";
my $cds_length;
while (<IR>){
	chomp;
	next if (/##/);
	my @line=split /\t/;
	print "line$.\n";
	$cds_length=length $line[-1];

	if (/Locus/){
		print NEW "Sample_name\t"."$_\t"."cds_length\n";
	}
	else{	
		print NEW "$sam[0]\t"."$_\t"."$cds_length\n";
	}
}
close IR;
close NEW;
print "End:add the sample name at the first column.\tplus_sample_name_v3.pl\n";
