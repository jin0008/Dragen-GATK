#!/bin/bash
# arg1: number of threads
# to run: 
# chmod +x dragen.sh
# <path>/trim.sh <number of threads>
# Example: ./dragen.sh 40

compression_level=5
java_opt=-Xmx128g
REF=/media/hanjinu/PM883/db/refs/hg38_broad/Homo_sapiens_assembly38.fasta

for f in *_R1.fastq.gz # for each n

do
    n=${f%%_R1.fastq.gz} # strip part of file name
    dragen-os -r /media/hanjinu/PM883/db/refs/hg38_broad/ \
    --num-threads $1 -1 ${n}_R1.fastq.gz -2 ${n}_R2.fastq.gz | \
    samtools view -@ $1 -Sb | samtools sort -n -@ $1 > ${n}.mapped.bam

    #markduplicatesSpark (#--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \ for patterned flow cell model)
    gatk --java-options "${java_opt}" MarkDuplicatesSpark \
    -I ${n}.mapped.bam \
    -O ${n}.dedup.bam \
    -M ${n}.mark_dup_metrics.txt \
    --spark-master local[$1] \
    --optical-duplicate-pixel-distance 2500 \
    --tmp-dir $PWD \
    --verbosity ERROR

    rm -rf ${n}.mapped.bam
    #rm -rf ${n}.mapped.bai

    #samtools index -@ 1 ${n}.markduplicates.bam

    #gatk CalibrateDragstrModel
    gatk --java-options "${java_opt}" CalibrateDragstrModel \
    -I ${n}.dedup.bam \
    -O ${n}.CalibrateDragstrModel \
    -R ${REF} \
    -str /media/hanjinu/PM883/db/refs/hg38_broad/hg38.str.zip \
    --threads $1 \
    --tmp-dir $PWD

    #Haplotypecaller
    # gatk --java-options "${java_opt}" HaplotypeCaller \
    # -R ${REF} \
    # -I ${n}.dedup.bam \
    # -O ${n}.g.vcf.gz \
    # -bamout ${n}.bamout.bam
    # --tmp-dir $PWD \
    # --dragen-mode true \
    # -ERC GVCF \
    # --native-pair-hmm-threads 8 \
    # --dragstr-params-path ${n}.CalibrateDragstrModel


    #Haplotypecaller (M>Y>22>19>X>21>18>15>14>20>17>13>12>11>9>8>7>5>6>16>3>4>)
    for i in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y","M"}
    do
        gatk --java-options "-Xmx16G" HaplotypeCaller -R /media/hanjinu/PM883/db/refs/hg38_broad/Homo_sapiens_assembly38.fasta -L chr${i} -I ${n}.dedup.bam -bamout ${n}_chr${i}.bamout.bam -O ${n}_chr${i}.g.vcf.gz --tmp-dir $PWD --dragen-mode true -ERC GVCF --dragstr-params-path ${n}.CalibrateDragstrModel &
    done

    wait

    #GatherVcfs
    gatk GatherVcfs \
    -I ${n}_chr1.g.vcf.gz \
    -I ${n}_chr2.g.vcf.gz \
    -I ${n}_chr3.g.vcf.gz \
    -I ${n}_chr4.g.vcf.gz \
    -I ${n}_chr5.g.vcf.gz \
    -I ${n}_chr6.g.vcf.gz \
    -I ${n}_chr7.g.vcf.gz \
    -I ${n}_chr8.g.vcf.gz \
    -I ${n}_chr9.g.vcf.gz \
    -I ${n}_chr10.g.vcf.gz \
    -I ${n}_chr11.g.vcf.gz \
    -I ${n}_chr12.g.vcf.gz \
    -I ${n}_chr13.g.vcf.gz \
    -I ${n}_chr14.g.vcf.gz \
    -I ${n}_chr15.g.vcf.gz \
    -I ${n}_chr16.g.vcf.gz \
    -I ${n}_chr17.g.vcf.gz \
    -I ${n}_chr18.g.vcf.gz \
    -I ${n}_chr19.g.vcf.gz \
    -I ${n}_chr20.g.vcf.gz \
    -I ${n}_chr21.g.vcf.gz \
    -I ${n}_chr22.g.vcf.gz \
    -I ${n}_chrX.g.vcf.gz \
    -I ${n}_chrY.g.vcf.gz \
    -I ${n}_chrM.g.vcf.gz \
    -O ${n}.g.vcf.gz

    rm -rf ${n}_chr*.g.vcf.gz
    rm -rf ${n}_chr*.g.vcf.gz.tbi
    #rm -rf ${n}.mapped.bam

    samtools merge -@ 1 ${n}.bamout.bam *chr*.bamout.bam
    samtools index -@ 1 ${n}.bamout.bam
    rm -rf ${n}_chr*.bamout.bam
    rm -rf ${n}_chr*.bamout.bai

    gatk IndexFeatureFile \
    -I ${n}.g.vcf.gz \
    --tmp-dir $PWD

    #n=URD_SH_005_3

    #create tdf files
    igvtools count ${n}.dedup.bam ${n}.dedup.bam.tdf hg38

done

