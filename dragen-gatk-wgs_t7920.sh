Sample=5006224
thread=56
compression_level=5
java_opt=-Xmx64g
#REF=/media/hanjinu/PM883/db/refs/hg38_masked/hg38.analysisSet.fa
#REF=/media/hanjinu/PM883/db/refs/hg38_alt_masked/GRCh38_full_analysis_set_plus_decoy_hla.fa
#INTERVAL=/media/hanjinu/PM883/db/refs/interval_list/Twist_Exome/Twist_ComprehensiveExome_hg38.interval_list
#INTERVAL=/media/hanjinu/PM883/db/refs/interval_list/whole.exome.hg38.gene.interval_list
#INTERVAL=/media/hanjinu/PM883/db/refs/interval_list/IDT_exome/xgen-exome-hyb-panel-v2.hg38.interval_list
REF=/media/hanjinu/PM883/db/refs/hg38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna


#gatk SplitSamByNumberOfReads \
#-I ${Sample}.unmapped.bam \
#-O $PWD \
#--SPLIT_TO_N_FILES 12 \
#--TMP_DIR $PWD \84

#dragen-os --build-hash-table true --ht-reference  --output-directory $PWD

#gatk ComposeSTRTableFile \
#-R GRCh38_full_analysis_set_plus_decoy_hla.fa \
#-O hg38.str.zip

dragen-os -r /media/hanjinu/PM883/db/refs/hg38_no_alt/ \
--num-threads ${thread} -1 ${Sample}_R1.fastq.gz -2 ${Sample}_R2.fastq.gz | \
samtools view -@ ${thread} -Sb | samtools sort -n -@ ${thread} > ${Sample}.mapped.bam


#markduplicatesSpark (#--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \ for patterned flow cell model)
gatk --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" MarkDuplicatesSpark \
-I ${Sample}.mapped.bam \
-O ${Sample}.dedup.bam \
-M ${Sample}.mark_dup_metrics.txt \
--conf 'spark.executor.cores=${thread}' \
--optical-duplicate-pixel-distance 2500 \
--tmp-dir $PWD \
--verbosity ERROR


#rm -rf ${Sample}.mapped.bam
#rm -rf ${Sample}.mapped.bai

#samtools index -@ {thread} ${Sample}.markduplicates.bam


#gatk CalibrateDragstrModel (BETA)
gatk CalibrateDragstrModel \
-I ${Sample}.dedup.bam \
-O ${Sample}.CalibrateDragstrModel \
-R ${REF} \
-str /media/hanjinu/PM883/db/refs/hg38_no_alt/hg38.str.zip \
--tmp-dir $PWD 

#Haplotypecaller
gatk HaplotypeCaller \
-R ${REF} \
-I ${Sample}.dedup.bam \
-O ${Sample}.vcf.gz \
-bamout ${Sample}.bamout.bam \
--tmp-dir $PWD \
--dragen-mode true \
--dragstr-params-path ${Sample}.CalibrateDragstrModel

gatk VariantFiltration \
-V ${Sample}.vcf.gz \
--filter-expression "QUAL < 10.4139" \
--filter-name "DRAGENHardQual" \
-O ${Sample}.filtered.vcf.gz

