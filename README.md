# Ion_ir_add_script
Annotation for data from Ion torrent suite software, and statistics for gene variants.


#cp ~/ir_add_v3/* $wording_dir # copy scripts to working dir.

##1.注释文件*.tsv批量处理
#*.tsv file is variant calling results from ion torrent suite(v4.12).
for i in *.tsv
do
        echo "#############################################################"
        echo "#############################################################"
        echo "#############################################################"
        echo "$i analysis begins."
        perl filter_v3.pl $i
        perl locus_modify_v3.pl $i.filter
        perl filter_vcf_kg_base_v3.pl ${i}.filter.alt IR
        perl plus_dbsnp_af_v3.pl ${i}.filter.alt.kg
        perl add_icgc_score2word_v3.pl ${i}.filter.alt.kg.dbsnp
        perl exon_intron_ncbi_v3.pl  ${i}.filter.alt.kg.dbsnp.icgc
        perl add_bothside_sequence_v3.pl ${i}.filter.alt.kg.dbsnp.icgc.inex
        perl base2aa_v3.pl ${i}.filter.alt.kg.dbsnp.icgc.inex.seq
        perl plus_sample_name_v3.pl ${i}.filter.alt.kg.dbsnp.icgc.inex.seq.prot
        perl plus_var_ref_v3.pl ${i}.filter.alt.kg.dbsnp.icgc.inex.seq.prot.name
        perl plus_novel_v3.pl ${i}.filter.alt.kg.dbsnp.icgc.inex.seq.prot.name.var_ref
        perl adjust_order_v3.pl ${i}.filter.alt.kg.dbsnp.icgc.inex.seq.prot.name.var_ref.novel

        new=${i%%.*}
        echo ${new}_ir_add.xls
        cp ${i}.filter.alt.kg.dbsnp.icgc.inex.seq.prot.name.var_ref.novel.adjust ${new}_ir_add.xls
        rm -rf $i.filter*
done

rm *.pl *.sh #remove scripts.

