#!/user/bin/perl -w
#filter REF or Nocall sites of vcf
print "\nStart:filter REF or Nocall sites of vcf.\tfilter_v3.pl\n";
use strict;
my $clumn=0;
open (IR,$ARGV[0]) || die;
open (IR_filter,">$ARGV[0].filter")||die;
while (<IR>){
	chomp;
	next if (/^##/);
	if (/^Locus/){
		my @line=split /\t/;
		if($line[8] eq 'Info'){
			splice(@line,8,1);
			$clumn=1;
		}
		print IR_filter join("\t",@line),"\n";
		next;
	}
	print "line$.\n";
	my @line = split/\t/;
	if($clumn==1){
		splice(@line,8,1);
	}
	if ($line[3]=~ "REF"||$line[3]=~"NOCALL"){
		next;	
#		print IR_filter "$_\n";
	}
	if($line[5]=~/FBXW7/){$line[10]='-';}
	if($line[5]=~/EGFR/){$line[10]='+';}
	if($line[5]=~/CSF1R/){$line[10]='-';}
	if($line[5] eq 'CSF1R,HMGXB3' && $line[11] eq 'NM_014983.2, NM_005211.3'){
		$line[5]='CSF1R';$line[11]='NM_005211.3';
	}
	if($line[12]=~/,/){	
		print "Maybe 2 different mutations: $line[1]\t$line[2]\t$line[12]\n";
		my @c=split(",",$line[12]);
		my @genotype=split("/",$line[1]);
		if($genotype[0] eq $genotype[1] || $genotype[0] eq $line[2]){
			print "Only 1 mutation\t$line[1]\t$line[2]\n";
			next;
		}
		my @type=split(",",$line[3]);
		my @location=split(",",$line[6]);
		my @nm=split(",",$line[11]);
		my @p=split(",",$line[13]);
		my @var_eff=split(",",$line[14]);
		foreach my $i(0..$#c){
			$line[1]=$line[2]."/".$genotype[$i];
			if(undef $type[$i]){
				$line[3]=$type[0];
			}else{
				$line[3]=$type[$i];
			}
			print "Divide the mutation to: $line[3]\t$type[$i]\n";
			$line[6]=$location[$i];
			$line[11]=$nm[$i];
			$line[12]=$c[$i];
			$line[13]=$p[$i];
			$line[14]=$var_eff[$i];
			foreach my $j(1,2,6,11,12,13,14){
				$line[$j]=~s/\s//g;
			}
			print IR_filter join("\t",@line),"\n";
		}
	next;
	}
#	my @var=split(/\//,$line[1]);
#	print "$var[1]\n";
	print IR_filter join("\t",@line),"\n";
	
}
close IR;	
close IR_filter;
print "End:filter REF or Nocall sites of vcf.\tfilter_v3.pl\n";
