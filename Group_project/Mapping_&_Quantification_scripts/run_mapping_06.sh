referenceFile=$SCRATCH/np_data/reference.fa
referenceTranscriptome=$SCRATCH/np_data/transcriptome_reference.fa
annotationFile=$SCRATCH/np_data/annotation.gtf

experimentFolder="/cluster/scratch/ochsneto/np_data/np_data/1_STA_ProjektEvan/NP06-RNA-20191003-dKD1/NP06-RNA-20191003-dKD1/20191003_0817_MN31533_FAK97601_6b2f9ce7"


#Join all the fastq files in one pipe
cat $experimentFolder/fastq_pass/*.fastq > $experimentFolder/NP06-dKD1.fastq 

mkdir "$experimentFolder/alignments"

echo "Performing genome alignment on NP06-dKD1"
#Mapping with minimap2
minimap2 -ax splice -uf -k14 $referenceFile  $experimentFolder/NP06-dKD1.fastq > $experimentFolder/alignments/genomic-aln.sam
# Convert sam to bam
samtools view -bS  $experimentFolder/alignments/genomic-aln.sam >  $experimentFolder/alignments/genomic-aln.bam
# Sort bam
samtools sort -o  $experimentFolder/alignments/sorted-genomic-aln.bam  $experimentFolder/alignments/genomic-aln.bam 
# Index the sorted bam
samtools index $experimentFolder/alignments/sorted-genomic-aln.bam
# Create bam with primary alignment only
samtools view -b -h -F 2308 $experimentFolder/alignments/sorted-genomic-aln.bam > $experimentFolder/alignments/primary-genomic-aln.bam
# Convert to bed12 (this is the input for FLAIR)
bedtools bamtobed -bed12 -i $experimentFolder/alignments/primary-genomic-aln.bam > $experimentFolder/alignments/primary-genomic-aln.bed12

echo "Performing transcriptome alignment on NP06-dKD1"
minimap2 -ax map-ont -N 15 $referenceTranscriptome $experimentFolder/NP06-dKD1.fastq > $experimentFolder/alignments/transcriptomic-aln.sam
# Convert sam to bam
samtools view -bS $experimentFolder/alignments/transcriptomic-aln.sam > $experimentFolder/alignments/transcriptomic-aln.bam
# Sort bam
samtools sort -o $experimentFolder/alignments/sorted-transcriptomic-aln.bam $experimentFolder/alignments/transcriptomic-aln.bam
# Index the sorted bam
samtools index $experimentFolder/alignments/sorted-transcriptomic-aln.bam
# Create bam with primary alignment and secondary only
samtools view -h -F 2052 $experimentFolder/alignments/sorted-transcriptomic-aln.bam > $experimentFolder/alignments/primary-transcriptomic-aln.bam
# Conver to bed12
bedtools bamtobed -bed12 -i $experimentFolder/alignments/primary-transcriptomic-aln.bam > $experimentFolder/alignments/primary-transcriptomic-aln.bed12

mkdir "$experimentFolder/quantifications"

echo "Performing gene quantification"
# Quantifty genes with featureCounts
featureCounts -L -a $annotationFile -o $experimentFolder/quantifications/primary-genomic-quant $experimentFolder/alignments/genomic-aln.bam --primary
    # All alignments for differential expression analysis
featureCounts -L -a $annotationFile -o $experimentFolder/quantifications/genomic-quant $experimentFolder/alignments/genomic-aln.bam

echo "Performing transcript quantification"
# Quantifying transcripts with NanoCount
NanoCount -i $experimentFolder/alignments/transcriptomic-aln.bam -p align_score -a -b $experimentFolder/alignments/transcriptomic-aln-nanocount-filtered.bam --extra_tx_info -o $experimentFolder/quantifications/transcriptomic-quant.tsv

echo "Finished alignment and quantification for this sample"

