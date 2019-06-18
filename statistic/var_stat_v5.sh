#!/bin/bash -w

################################
#author=duzhao
#version=1.0
#date 2016-11-01
####################
cp /statistic/* ./ # copy scripts to working dir.
echo "variants statistic. Version 5.0. $1=sample_number $2=Type(such as lung_cancer)."
perl new_af_qual_v4.pl
perl left_normalization_var.pl
perl variants_fre_v5.pl
head -n 1 variants_fre_all.xls>variants_fre_all.xls.head
sed '1d' variants_fre_all.xls|sort -t "	" -k 7 -bnr  -o variants_fre_all.xls.nohead
cat variants_fre_all.xls.head variants_fre_all.xls.nohead >variants_fre_all.xls
rm -f variants_fre_all.xls.head variants_fre_all.xls.nohead
awk -F"\t" '{if($1 == "variants"||$7 == 0) {print}}' variants_fre_all.xls>variants_fre_nomut.xls
awk -F"\t" '{if($1 == "variants"||$7 == 100) {print}}' variants_fre_all.xls>variants_fre_nocall.xls
awk -F"\t" '{if($1 == "variants"||($7>0 && $7<=1)) {print}}' variants_fre_all.xls>variants_fre_mut.xls
awk -F"\t" '{if($1 == "variants"||$9 == "Novel") {print}}' variants_fre_mut.xls>variants_fre_mut_novel.xls
awk -F"\t" '{if($1 == "variants"||$9 == "Hotspot") {print}}' variants_fre_mut.xls>variants_fre_mut_hotspot.xls
awk -F"\t" '{if($1 == "variants"||$5>0.2) {print}}' variants_fre_mut.xls>variants_fre_rmnocall.xls
awk -F"\t" '{if($1 == "variants"||$9 == "Novel"){print}}' variants_fre_rmnocall.xls>variants_fre_rmnocall_novel.xls
awk -F"\t" '{if($1 == "variants"||$9 == "Hotspot"){print}}' variants_fre_rmnocall.xls>variants_fre_rmnocall_hotspot.xls
awk -F"\t" '{if($1 == "variants"||$8== "NA" || $8<=0.05){print}}' variants_fre_rmnocall.xls>variants_fre.xls
awk -F"\t" '{if($1 == "variants"||$9 == "Novel"){print}}' variants_fre.xls>variants_fre_novel.xls
awk -F"\t" '{if($1 == "variants"||$9 == "Hotspot"){print}}' variants_fre.xls>variants_fre_hotspot.xls
perl /results/duzhao_test/scripts/statistic/gene_fre_v4.pl variants_fre.xls gene_fre.xls $1
head -n 1 gene_fre.xls> gene_fre.xls.head
sed '1d' gene_fre.xls|sort -t "	" -k 4 -bnr -o gene_fre.xls.nohead
cat gene_fre.xls.head gene_fre.xls.nohead>gene_fre.xls
rm -f gene_fre.xls.head gene_fre.xls.nohead
mv variants_fre_all.xls $2_variants_fre_all.xls
mv variants_fre_mut.xls $2_variants_fre_mut.xls
mv variants_fre_mut_novel.xls $2_variants_fre_mut_novel.xls
mv variants_fre_mut_hotspot.xls $2_variants_fre_mut_hotspot.xls
mv variants_fre_nomut.xls $2_variants_fre_nomut.xls
mv variants_fre_nocall.xls $2_variants_fre_nocall.xls
mv variants_fre_rmnocall.xls $2_variants_fre_rmnocall.xls
mv variants_fre_rmnocall_novel.xls $2_variants_fre_rmnocall_novel.xls
mv variants_fre_rmnocall_hotspot.xls $2_variants_fre_rmnocall_hotspot.xls
mv variants_fre.xls $2_variants_fre.xls
mv variants_fre_novel.xls $2_variants_fre_novel.xls
mv variants_fre_hotspot.xls $2_variants_fre_hotspot.xls
mv gene_fre.xls $2_gene_fre.xls
#wc -l $2_*.xls >$2_count_stat.txt
perl sites_count.pl $2_stat.txt

rm *.pl *.sh #remove scripts.
