#!/usr/bin/perl
#modify the aa information of column 'Amino Acid Change'.
print "\nStart:modify the aa information of column 'Amino Acid Change'.\tbase2aa_v3.pl\n";
use warnings; 
use strict;
open(IR,$ARGV[0])||die;
open(NEW,">$ARGV[0].prot")||die;
while(<IR>){
	chomp;
	my $dna="N";
#	print "##########################################################\n";
	if(!/^chr/){print NEW "$_\n";next;}
	my @line=split /\t/;
	if($line[11] eq 'NM_001163213.1' && $line[12] eq 'c.2342_2343insC'){
		$line[13]='p.Ser782Leufs*37';
		print NEW join("\t",@line),"\n";
		next;
	}
	if($line[11] eq 'NM_005228.3' && $line[12] eq 'c.1441delT'){
		$line[13]='p.Phe481Leufs*11';
                print NEW join("\t",@line),"\n";
                next;
        }
	if($line[1] eq 'chr10:89685281' && $line[12] eq 'c.176C>G'){
		$line[13]='p.Ser58*';
                print NEW join("\t",@line),"\n";
                next;
        }
	if($line[14] eq 'synonymous' && $line[12]=~/(\d+)/){
	        if($line[10] eq '-'){
			$dna=substr(&seq_reverse($line[48]),6+($1%3),3);
		}else{
			$dna=substr($line[48],6+($1%3),3);
		}
		my $protein='';
		my $codon;
		for(my $i=0; $i<(length($dna)-2);$i+=3) {
			$codon=substr($dna,$i,3);
			$protein.=&codon2aa($codon);
		}
		my $p=int(($1-1)/3)+1;
		print "synonymous:\tphase=",$1%3,"\tstart=$1\tseq=$line[48]\tDNA=$dna\tProtein=$protein\n";
		$line[13]="p.$protein$p$protein";
	}elsif($line[14]=~/frameshift/ && $line[12]=~/(\d+)_(\d+)(ins|del)(.*)/){
		print "$line[12]\t$1\t$2\t$3\t$4\n",$line[48],"\n";
		if($3 eq 'ins'){
			$dna=substr($line[49],0,$1).lc($4).substr($line[49],$1);
			print "frameshift_ins\tstart=$1\tend=$2\t$3\tinseq=$4\n";
		}elsif($3 eq 'del'){
			$dna=substr($line[49],0,$1-1).substr($line[49],$2);
			print "frameshift_del\tstart=$1\tend=$2\t$3\tdelseq=$4\n";
		}
	}elsif($line[14]=~/frameshift/ && $line[12]=~/(\d+)(del)(.*)/){
		print "frameshift_del_one_base:$line[12]\tstart=$1\t$2\tdelseq=$3\n";
		$dna=substr($line[49],0,$1-1).substr($line[49],$1);
	}
	my $var_posit=0;
	my $stop_posit=0;
	if($line[13]=~/(\d+)fs/){
		$var_posit=$1;
		print "Raw p.XXX:$line[13]\tvar-position:$var_posit\n";
	}
	my $protein='';
	my $codon;
	my ($prot_alt,$i);
	print "line=$.\tdna=$dna\n";
	if($dna eq "N"){
		print NEW join("\t",@line),"\n";;
		print "line$.=no p.XXX or no modification of p.XXX \n";
		next;
	}
	for($i=0; $i<(length($dna)-2);$i+=3) {
	#	print "$i\tdna=$dna\n";
	        $codon=substr($dna,$i,3);
        	my $aa=&codon2aa($codon);
		$protein.=$aa;
#		print "raw:\t$line[13]\n";
		if(int($i/3)+1==$var_posit){
			$prot_alt=$aa;
			print "Alt:\tvar_position=$var_posit\tprotein_alt=$prot_alt\t";
		}#突变位置
		if($aa eq '_'){
			$stop_posit=$i/3+1;
			print "protein_stop_position=$stop_posit\tprotein_stop=$aa\n";
			last;
		}#终止氨基酸位置
	}
#	print "$line[13]*$stop_posit\n$protein\n";
	if($line[14]=~/frameshift/ && $i>=length($dna)-2){
		print "Stop codon over the protein!\n";
	}
	if($line[13]=~/p\.([A-Za-z]+)(\d+)fs/){
		print "Add protein stop position:c.XXX=$line[13]\tref_protein=$1\tposition=$2\talt_protein=$prot_alt\tstop_position=$stop_posit\n";
		if( $stop_posit==$var_posit){
			$line[13]=$line[13]."*";
		}else{
			my $end=$stop_posit-$2+1;
			$line[13]=~s/fs/${prot_alt}fs*${end}/g;
		}
		print "modify c.XXX=$line[13]\n";
	}
	print NEW join("\t",@line),"\n";
}
close IR;
close NEW;
#*****************************************************************************************#
# codon2aa #
# A subroutine to translate a DNA 3-character codon to an amino acid
# Version 3, using hash lookup
sub seq_reverse(){
        my ($seq)=@_;
        my $reverse;
        my %rever=('A'=>'T','T'=>'A','G'=>'C','C'=>'G','a'=>'t','t'=>'a','g'=>'c','c'=>'g');
        for(my $i=(length $seq)-1;$i>=0;$i--){
                $reverse.=$rever{(substr($seq,$i,1))};
        }
        print "Raw seq=$seq\treverse seq=$reverse\n";
        return $reverse;
}

sub codon2aa {
	my($codon) = @_; $codon = uc $codon;
#uc=uppercase;lc=lowercase
#也就是大小写转换，uc表示将所有的小写 转换为大写
#lc将所有的大写转换为小写 
	my(%genetic_code) = ('TCA' => 'S', # Serine 
			'TCC' => 'S', # Serine 
			'TCG' => 'S', # Serine 
			'TCT' => 'S', # Serine 
			'TTC' => 'F', # Phenylalanine 
			'TTT' => 'F', # Phenylalanine 
			'TTA' => 'L', # Leucine 
			'TTG' => 'L', # Leucine 
			'TAC' => 'Y', # Tyrosine 
			'TAT' => 'Y', # Tyrosine 
			'TAA' => '_', # Stop
			'TAG' => '_', # Stop 
			'TGC' => 'C', # Cysteine 
			'TGT' => 'C', # Cysteine 
			'TGA' => '_', # Stop 
			'TGG' => 'W', # Tryptophan 
			'CTA' => 'L', # Leucine 
			'CTC' => 'L', # Leucine 
			'CTG' => 'L', # Leucine 
			'CTT' => 'L', # Leucine 
			'CCA' => 'P', # Proline 
			'CCC' => 'P', # Proline 
			'CCG' => 'P', # Proline 
			'CCT' => 'P', # Proline 
			'CAC' => 'H', # Histidine 
			'CAT' => 'H', # Histidine 
			'CAA' => 'Q', # Glutamine 
			'CAG' => 'Q', # Glutamine 
			'CGA' => 'R', # Arginine 
			'CGC' => 'R', # Arginine 
			'CGG' => 'R', # Arginine 
			'CGT' => 'R', # Arginine 
			'ATA' => 'I', # Isoleucine 
			'ATC' => 'I', # Isoleucine 
			'ATT' => 'I', # Isoleucine 
			'ATG' => 'M', # Methionine 
			'ACA' => 'T', # Threonine 
			'ACC' => 'T', # Threonine 
			'ACG' => 'T', # Threonine 
			'ACT' => 'T', # Threonine 
			'AAC' => 'N', # Asparagine 
			'AAT' => 'N', # Asparagine 
			'AAA' => 'K', # Lysine 
			'AAG' => 'K', # Lysine 
			'AGC' => 'S', # Serine 
			'AGT' => 'S', # Serine 
			'AGA' => 'R', # Arginine 
			'AGG' => 'R', # Arginine 
			'GTA' => 'V', # Valine 
			'GTC' => 'V', # Valine 
			'GTG' => 'V', # Valine 
			'GTT' => 'V', # Valine 
			'GCA' => 'A', # Alanine 
			'GCC' => 'A', # Alanine 
			'GCG' => 'A', # Alanine 
			'GCT' => 'A', # Alanine 
			'GAC' => 'D', # Aspartic Acid 
			'GAT' => 'D', # Aspartic Acid 
			'GAA' => 'E', # Glutamic Acid 
			'GAG' => 'E', # Glutamic Acid 
			'GGA' => 'G', # Glycine 
			'GGC' => 'G', # Glycine 
			'GGG' => 'G', # Glycine 
			'GGT' => 'G', # Glycine 
			);
	        my(%genetic_code3) = ('TCA' => 'Ser', # Serine
                        'TCC' => 'Ser', # Serine
                        'TCG' => 'Ser', # Serine
                        'TCT' => 'Ser', # Serine
                        'TTC' => 'Phe', # Phenylalanine
                        'TTT' => 'Phe', # Phenylalanine
                        'TTA' => 'Leu', # Leucine
                        'TTG' => 'Leu', # Leucine
                        'TAC' => 'Tyr', # Tyrosine
                        'TAT' => 'Tyr', # Tyrosine
                        'TAA' => '_', # Stop
                        'TAG' => '_', # Stop
                        'TGC' => 'Cys', # Cysteine
                        'TGT' => 'Cys', # Cysteine
                        'TGA' => '_', # Stop
                        'TGG' => 'Trp', # Tryptophan
                        'CTA' => 'Leu', # Leucine
                        'CTC' => 'Leu', # Leucine
                        'CTG' => 'Leu', # Leucine
                        'CTT' => 'Leu', # Leucine
                        'CCA' => 'Pro', # Proline
                        'CCC' => 'Pro', # Proline
                        'CCG' => 'Pro', # Proline
                        'CCT' => 'Pro', # Proline
                        'CAC' => 'His', # Histidine
                        'CAT' => 'His', # Histidine
                        'CAA' => 'Gln', # Glutamine
                        'CAG' => 'Gln', # Glutamine
                        'CGA' => 'Arg', # Arginine
                        'CGC' => 'Arg', # Arginine
                        'CGG' => 'Arg', # Arginine
                        'CGT' => 'Arg', # Arginine
                        'ATA' => 'Ile', # Isoleucine
                        'ATC' => 'Ile', # Isoleucine
                        'ATT' => 'Ile', # Isoleucine
                        'ATG' => 'Met', # Methionine
                        'ACA' => 'Thr', # Threonine
			'ACC' => 'Thr', # Threonine
                        'ACG' => 'Thr', # Threonine
                        'ACT' => 'Thr', # Threonine
                        'AAC' => 'Asn', # Asparagine
                        'AAT' => 'Asn', # Asparagine
                        'AAA' => 'Lys', # Lysine
                        'AAG' => 'Lys', # Lysine
                        'AGC' => 'Ser', # Serine
                        'AGT' => 'Ser', # Serine
                        'AGA' => 'Arg', # Arginine
                        'AGG' => 'Arg', # Arginine
                        'GTA' => 'Val', # Valine
                        'GTC' => 'Val', # Valine
                        'GTG' => 'Val', # Valine
                        'GTT' => 'Val', # Valine
                        'GCA' => 'Ala', # Alanine
                        'GCC' => 'Ala', # Alanine
                        'GCG' => 'Ala', # Alanine
                        'GCT' => 'Ala', # Alanine
                        'GAC' => 'Asp', # Aspartic Acid
                        'GAT' => 'Asp', # Aspartic Acid
                        'GAA' => 'Glu', # Glutamic Acid
                        'GAG' => 'Glu', # Glutamic Acid
                        'GGA' => 'Gly', # Glycine
                        'GGC' => 'Gly', # Glycine
                        'GGG' => 'Gly', # Glycine
                        'GGT' => 'Gly', # Glycine
                        );
	if(exists $genetic_code3{$codon}) { 
		return $genetic_code3{$codon}; 
	} else { 
		print STDERR "Bad codon \"$codon\"!!\n"; exit; 
	}
}
#***************************************************************************************** #
print "End:modify the aa information of column 'Amino Acid Change'.\tbase2aa_v3.pl\n";
