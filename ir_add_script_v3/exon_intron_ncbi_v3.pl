#!/usr/bin/perl -w
#add intron or exon information and modify the 'Variant Effect' column
print "Start:add intron or exon information and modify the 'Variant Effect' column.\texon_intron_ncbi_v3.pl\n";
use strict;
open(GFF,"/results/duzhao_test/Analysis/IR-reporter_Analysis_data/gff_to_intron_extron/ref_GRCh37.p5_top_level.gff3")||die;
my ($chr,$gname,$trscript,$start,$end);
my (%exon,%cds,%rna);
while(<GFF>){
	chomp;
	if(/chromosome=(\d+);/){
		$chr=$1;
		next;
	}
	my @line=split /\t/;
	next if (scalar @line<8);
	if($line[2]=~/mRNA|exon/ && $line[8]=~/;gene=(BRCA1|BRCA2|ABL1|AKT1|ALK|APC|ATM|BRAF|CHD1|CDKN2A|CSF1R|CTNNB1|EGFR|ERBB2|ERBB4|EZH2|FBXW7|FGFR1|FGFR2|FGFR3|FLT3|GNA11|GNAQ|GNAS|HNF1A|HRAS|IDH1|IDH2|JAK2|JAK3|KDR|KIT|KRAS|MET|MLH1|MPL|NOTCH1|NPM1|NRAS|PDGFRA|PIK3CA|PTEN|PTPN11|RB1|RET|SMAD4|SMARCB1|SMO|SRC|STK11|TP53|VHL|DDR2|MAP2K1);.*transcript_id=(.*)\./){#50 gene,CHPv2
#	if($line[8]=~/;gene=(APC|ERBB2|ERBB4);.*transcript_id=(.*)\./){
#		print "$_\tchr$chr\n";
		$gname=$1;
		$trscript=$2;
#		print "aaa$line[2]\n";
		if ($line[2] eq "mRNA" && ($line[8]=~/ID=(rna\d+)/)){
			$rna{$1}=$trscript;	next;
		}
		if($line[6] eq "-" && $line[8]=~/ID=id(\d+);/){#antisense strand
			push @{$exon{"$gname:$trscript"}},($line[4],$line[3]);
		}else{#sense strand
			push @{$exon{"$gname:$trscript"}},($line[3],$line[4]);
		}
	}elsif($line[2]=~/CDS/ && $line[8]=~/^ID=(cds\d+);.*Parent=(rna\d+);/){
		if(exists $rna{$2}){
			if($line[6] eq "-"){
				push @{$cds{$gname.":".$rna{$2}}},($line[4],$line[3]);
			}else{
				push @{$cds{$gname.":".$rna{$2}}},($line[3],$line[4]);
			}
		}
	}
}
close GFF;
foreach(sort keys %cds){
#	print "$_\t",scalar @{$cds{$_}},"\t@{$cds{$_}}\n";
}
open(VCF,"$ARGV[0]")||die;
open(NEW,">$ARGV[0].inex")||die;
my %mat=('A'=>'T','T'=>'A','G'=>'C','C'=>'G');
while(<VCF>){
	chomp;
	if (/^#.*/){
        print NEW "$_\n";
        next;
        }
	if (/^Locus/){
       		print NEW "$_\tCDS\tEXIN\tcds_start_end\n";
		next;
	}
	my @line=split /\t/;
	print "line$.\tMutation site:$line[0]\t$line[1]\t$line[2]\t$line[3]\n";
	my $loc=(split(/:/,$line[0]))[1];
	my $nm=(split(/\./,$line[11]))[0];#NM_00124.2=>NM_00124
	$line[5]=~s/.*(BRCA1|BRCA2|ABL1|AKT1|ALK|APC|ATM|BRAF|CHD1|CDKN2A|CSF1R|CTNNB1|EGFR|ERBB2|ERBB4|EZH2|FBXW7|FGFR1|FGFR2|FGFR3|FLT3|GNA11|GNAQ|GNAS|HNF1A|HRAS|IDH1|IDH2|JAK2|JAK3|KDR|KIT|KRAS|MET|MLH1|MPL|NOTCH1|NPM1|NRAS|PDGFRA|PIK3CA|PTEN|PTPN11|RB1|RET|SMAD4|SMARCB1|SMO|SRC|STK11|TP53|VHL|DDR2|MAP2K1).*/$1/;
	if(exists $exon{"$line[5]:$nm"}){
		my @cds_posi=@{$cds{"$line[5]:$nm"}};
		my @exon_posi=@{$exon{"$line[5]:$nm"}};
		my $cds_length;
		my $cds_num='.';
		my ($ref,$alt)=($line[2],(split(/\//,$line[1]))[1]);
		if(($loc>$exon_posi[0] && $loc<$cds_posi[0]) || ($loc<$exon_posi[0] && $loc>$cds_posi[0])){#5-utr,sense || antisense
#			print "5-utr,sense||antisense:$line[6]\n";
			$line[14]='5-utr';	
		}elsif(($loc<$exon_posi[$#exon_posi] && $loc>$cds_posi[$#cds_posi]) || ($loc>$exon_posi[$#exon_posi] && $loc<$cds_posi[$#cds_posi])){#3-utr,sense||antisense
			$line[14]='3-utr';	
		}
		for(my $i=0;$i<$#cds_posi;$i++){#CDS
			if($i%2==0){$cds_length+=abs($cds_posi[$i+1]-$cds_posi[$i])+1;}
			if($cds_posi[$i]<=$loc && $loc<=$cds_posi[$i+1]){#sense
				my $up_dis=$loc-$cds_posi[$i];
        	                my $down_dis=$cds_posi[$i+1]-$loc;
				if($i%2!=0){#intron
					print "Sense:intron:$line[6]\t$line[12]\tREF_length=",length($ref),"\tALT_length=",length($alt),"\t";
					if($up_dis<=$down_dis){#close to upstream splicesite
						if(length($ref) > length $alt){#deletion
							if(length($ref)==2){
								$line[12]="c.".$cds_length."+".($up_dis+1)."del".substr($ref,1);
							}else{
								$line[12]="c.".$cds_length."+".($up_dis+1)."_".$cds_length."+".($up_dis+1+(length($ref)-2))."del".substr($ref,1);
							}
						}elsif(length($ref) < length $alt){#insertion
							$line[12]="c.".$cds_length."+".($up_dis)."_".$cds_length."+".($up_dis+1)."ins".substr($alt,1);
						}else{#SNV,MNV
							$line[12]="c.".$cds_length."+".$up_dis.$ref.">".$alt;
						}
					}else{#close to downstream splicesite
						if(length($ref) > length $alt){#deletion
							if(length($ref)==2){
								$line[12]="c.".($cds_length+1)."-".($down_dis+1)."del".substr($ref,1);
							}else{
	                                                        $line[12]="c.".($cds_length+1)."-".($down_dis+1+(length($ref)-2))."_".($cds_length+1)."-".($down_dis+1)."del".substr($ref,1);
							}
                                                }elsif(length($ref) < length $alt){#insertion
							$line[12]="c.".($cds_length+1)."-".($down_dis)."_".($cds_length+1)."-".($down_dis-1)."ins".substr($alt,1);
						}else{#SNV,MNV
							$line[12]="c.".($cds_length+1)."-".$down_dis.$ref.">".$alt;
						}
					}
					print "modify:$line[12]\n";
				}else{
					print "Sense:exon:$line[6]\t$line[12]\tREF_length=",length($ref),"\tALT_length=",length($alt),"\t";
					if($up_dis<=$down_dis || $up_dis>$down_dis){#both upstream and downstream are the same
                                                if(length($ref) > length $alt){#deletion
							if($line[6]=~/splicesite|exonic|utr/){#splicesite need add c.xxx
								if(length($ref)==2){
									$line[12]="c.".($cds_length-$down_dis+1)."del".substr($ref,1);
								}else{
									$line[12]="c.".($cds_length-$down_dis+1)."_".($cds_length-$down_dis+1+(length($ref)-2))."del".substr($ref,1);
								}
							}
							$cds_num="CDS".int($i/2+1).",".($cds_length-($cds_posi[$i+1]-$cds_posi[$i]))."-".$cds_length.",(".($up_dis+1).",".($down_dis-(length($ref)-1)).")";
                                                }elsif(length($ref) < length $alt){#insertion
                                                        if($line[6]=~/splicesite|exonic|utr/){#splicesite need add c.xxx
								$line[12]="c.".($cds_length-$down_dis)."_".($cds_length-$down_dis+1)."ins".substr($alt,1);
							}
							$cds_num="CDS".int($i/2+1).",".($cds_length-($cds_posi[$i+1]-$cds_posi[$i]))."-".$cds_length.",(".$up_dis.",".($down_dis-1).")";
                                                }else{#SNV,MNV
                                                        if($line[6]=~/splicesite|exonic|utr/){#splicesite need add c.xxx
								$line[12]="c.".($cds_length-$down_dis).$ref.">".$alt;
							}
							$cds_num="CDS".int($i/2+1).",".($cds_length-($cds_posi[$i+1]-$cds_posi[$i]))."-".$cds_length.",(".$up_dis.",".($down_dis-(length($ref)-1)).")";
                                                }
                                        }
					print "modify:$line[12]\n";
				}
				last;
			}elsif($cds_posi[$i]>=$loc && $loc>=$cds_posi[$i+1]){#antisense
				my $up_dis=$cds_posi[$i]-$loc;
				my $down_dis=$loc-$cds_posi[$i+1];
				if($i%2!=0){#intron
					print "Antisense:intron:$line[6]\t$line[12]\tREF_length=",length($ref),"\tALT_length=",length($alt),"\t";
					if($up_dis<=$down_dis){#close to upstream splicesite
						if(length($ref) > length $alt){#deletion
                                                        if(length($ref)==2){
								$line[12]="c.".$cds_length."+".($up_dis-1)."del".&rev(substr($ref,1));
							}else{
								$line[12]="c.".$cds_length."+".($up_dis-1-(length($ref)-2))."_".$cds_length."+".($up_dis-1)."del".&rev(substr($ref,1));
							}
                                                }elsif(length($ref) < length $alt){#insertion
							$line[12]="c.".$cds_length."+".($up_dis-1)."_".$cds_length."+".$up_dis."ins".&rev(substr($alt,1));
						}else{#SNV MNV
							$line[12]="c.".$cds_length."+".$up_dis.&rev($ref).">".&rev($alt);
						}
                                        }else{#close to downstream splicesite
						if(length($ref) > length $alt){#deltion
							if(length($ref)==2){
								$line[12]="c.".($cds_length+1)."-".($down_dis+1)."del".&rev(substr($ref,1));
							}else{
								$line[12]="c.".($cds_length+1)."-".($down_dis+1+(length($ref)-2))."_".($cds_length+1)."-".($down_dis+1)."del".&rev(substr($ref,1));#print "REF:$ref\t",substr($ref,1),"\n";
							}
                                                }elsif(length($ref) < length $alt){#insertion
                                                        $line[12]="c.".($cds_length+1)."-".($down_dis+1)."_".($cds_length+1)."-".($down_dis)."ins".&rev(substr($alt,1));#print "ALT:$alt\t",substr($alt,1),"\n";
                                                }else{#SNV MNV
							$line[12]="c.".($cds_length+1)."-".$down_dis.&rev($ref).">".&rev($alt);
						}
                                        }
					print "modify:$line[12]\n";
				}else{
					print "Antisense:exon:$line[6]\t$line[12]\tREF_length=",length($ref),"\tALT_length=",length($alt),"\t";
                                        if($up_dis<=$down_dis || $up_dis>$down_dis){#both upstream and downstream are the same
                                                if(length($ref) > length $alt){#deletion
                                                        if($line[6]=~/splicesite|exonic|utr/){#splicesite need add c.xxx
							#	print "length ref=",length($ref),"\t";
                                                                if(length($ref)==2){
									$line[12]="c.".($cds_length-$down_dis-1)."del".&rev(substr($ref,1));
								}else{
									$line[12]="c.".($cds_length-$down_dis-1-(length($ref)-2))."_".($cds_length-$down_dis-1)."del".&rev(substr($ref,1));
                                                                }
                                                        }
                                                        $cds_num="CDS".int($i/2+1).",".($cds_length-($cds_posi[$i]-$cds_posi[$i+1]))."-".$cds_length.",(".($up_dis+1).",".($down_dis+1-length($ref)-1).")";
                                                }elsif(length($ref) < length $alt){#insertion
                                                        if($line[6]=~/splicesite|exonic|utr/){#splicesite need add c.xxx
                                                                $line[12]="c.".($cds_length-$down_dis-1)."_".($cds_length-$down_dis)."ins".&rev(substr($alt,1));
                                                        }
                                                        $cds_num="CDS".int($i/2+1).",".($cds_length-($cds_posi[$i]-$cds_posi[$i+1]))."-".$cds_length.",(".$up_dis.",".$down_dis.")";
                                                }else{#SNV,MNV
                                                        if($line[6]=~/splicesite|exonic|utr/){#splicesite need add c.xxx
                                                                $line[12]="c.".($cds_length-$down_dis+(length($ref)-1)).&rev($ref).">".&rev($alt);
                                                        }
                                                        $cds_num="CDS".int($i/2+1).",".($cds_length-($cds_posi[$i]-$cds_posi[$i+1]))."-".$cds_length.",(".$up_dis.",".($down_dis-(length($ref)-1)).")";
                                                }
                                        }
					print "modify:$line[12]\n";
				}
                                last;
                        }elsif($cds_posi[0]<$cds_posi[1]){
				if($loc<$cds_posi[0]){#sense 5_utr
                                        if(length($ref) > length $alt){#deletion
                                                $line[12]="c.-".($cds_posi[0]-$loc)."del".substr($ref,1);
                                        }elsif(length($ref) < length $alt){#insertion
                                                $line[12]="c.-".($cds_posi[0]-$loc)."ins".substr($alt,1);
                                        }else{
                                                $line[12]="c.-".($cds_posi[0]-$loc).$ref.">".$alt;
                                        }
				}elsif($loc>$cds_posi[$#cds_posi]){#sense 3_utr
                                        if(length($ref) > length $alt){#deletion
                                                $line[12]="c.*".($loc-$cds_posi[$#cds_posi])."del".substr($ref,1);
                                        }elsif(length($ref) < length $alt){#insertion
                                                $line[12]="c.*".($loc-$cds_posi[$#cds_posi])."ins".substr($alt,1);
                                        }else{
                                                $line[12]="c.*".($loc-$cds_posi[$#cds_posi]).$ref.">".$alt;
                                        }
				}
			}elsif($cds_posi[0]>$cds_posi[1]){
                               if($loc>$cds_posi[0]){#antisense 5_utr
                                        if(length($ref) > length $alt){#deletion
                                                $line[12]="c.-".($loc-$cds_posi[0])."del".&rev(substr($ref,1));
                                        }elsif(length($ref) < length $alt){#insertion
                                                $line[12]="c.-".($loc-$cds_posi[0])."ins".&rev(substr($alt,1));
                                        }else{
                                                $line[12]="c.-".($loc-$cds_posi[0]).&rev($ref).">".&rev($alt);
                                        }
                               }elsif($loc<$cds_posi[$#cds_posi]){#antisense 3_utr
                                        if(length($ref) > length $alt){#deletion
                                                $line[12]="c.*".($cds_posi[$#cds_posi]-$loc)."del".&rev(substr($ref,1));
                                        }elsif(length($ref) < length $alt){#insertion
                                                $line[12]="c.*".($cds_posi[$#cds_posi]-$loc)."ins".&rev(substr($alt,1));
                                        }else{
                                                $line[12]="c.*".($cds_posi[$#cds_posi]-$loc).&rev($ref).">".&rev($alt);
                                        }
                               }
			}
		}
                if($line[12]=~/\d*_\d*(ins|del)(\w*)/&&(length($2)%3==0)){#exon-indel,such as 1496_1497insACG
			print "exon indel\t$line[12]\n";
                        $line[14]='CDS-indel';
                }elsif($line[12]=~/\d*(\+|-)(\d*)_\d*(\+|-)(\d*)/){#intron,such as 1496+25_1497+18insACG
                        print "intron indel\t$line[12]\t$2\t$4\n";
                        if($2>20 && $4>20){
                                $line[14]='intron';
                        }else{
                                $line[14]='splice';
                        }
                }elsif($line[12]=~/\d*(\+|-)(\d*)/){#intron,such as 1496+41insACG
                        print "intron indel\t$line[12]\t$2\n";
			if($2>20){
                                $line[14]='intron';
                        }else{
                                $line[14]='splice';
                        }
                }
		my $i;
		for($i=0;$i<$#exon_posi;$i++){#exon
                        if($exon_posi[$i]<=$loc && $loc<=$exon_posi[$i+1]){#sense
                                my $up_dis=$loc-$exon_posi[$i];
                                my $down_dis=$exon_posi[$i+1]-$loc;
                                if($i%2==0){
                                        print NEW join("\t",@line,$cds_num,"EX".int($i/2+1)),"\t@cds_posi\n";
                                }else{
                                        if($up_dis<=$down_dis){
                                                print NEW join("\t",@line,$cds_num,"IN".int($i/2+1)).",UP".$up_dis."|down".$down_dis,"\t@cds_posi\n";
                                        }else{
                                                print NEW join("\t",@line,$cds_num,"IN".int($i/2+1)).",up".$up_dis."|DOWN".$down_dis,"\t@cds_posi\n";
                                        }
                                }
                                last;
                        }elsif($exon_posi[$i]>=$loc && $loc>=$exon_posi[$i+1]){#antisense
                                my $up_dis=$exon_posi[$i]-$loc;
                                my $down_dis=$loc-$exon_posi[$i+1];
                                if($i%2==0){
                                        print NEW join("\t",@line,$cds_num,"EX".int($i/2+1)),"\t@cds_posi\n";
                                }else{
                                        if($up_dis<=$down_dis){
                                                print NEW join("\t",@line,$cds_num,"IN".int($i/2+1)).",UP".$up_dis."|down".$down_dis,"\t@cds_posi\n";
                                        }else{
                                                print NEW join("\t",@line,$cds_num,"IN".int($i/2+1)).",up".$up_dis."|DOWN".$down_dis,"\t@cds_posi\n";
                                        }
                                }
                                last;
                        }
                }
		if($i>=$#exon_posi){
			print NEW join("\t",@line,".\t.\t.\n");
		}
	}else{
		print "not exists $line[5]:$nm\n"; 
		print NEW join("\t",@line,".\t.\t.\n");
	}
}
close VCF;
close NEW;
=pod
sub indel(){
	my ($ref,$alt)=@_;
	my $indel;
#	print"$ref,$alt\n";
	if(length($ref)<length $alt){
		for(my $i=0;$i<=length($ref); $i++){
			if(substr($ref,$i,1) ne substr($alt,$i,1)){
				$indel=substr($alt,$i);
				last;
			}
		}
	}
#	print "indel:$indel\n";
	return $indel;
}
=cut
sub rev(){
	my ($seq)=@_;
#	print "parameter:@_\t$seq\n";
	my $anti='';
	my %mat=('A'=>'T','T'=>'A','G'=>'C','C'=>'G');
	for(my $i=(length $seq)-1;$i>=0;$i--){
		$anti.=$mat{substr($seq,$i,1)};
	}
#	print "$seq\t$anti\n";
	return $anti;
}
print "End:add intron or exon information and modify the 'Variant Effect' column.\texon_intron_ncbi_v3.pl\n";
