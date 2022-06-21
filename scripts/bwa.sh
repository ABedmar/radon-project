module load bwa java picard samtools skewer/0.2.2

folder=$1
mkdir $folder/sam $folder/bam_raw $folder/bam_sort/ $folder/bam/ $folder/QC/ $folder/tmp/
#bed=$(ls $folder/INPUTS/*.bed)
genome=/slgpfs/projects/idib57/data/genome/Rattus/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa

picard="java -Xmx132G -jar /apps/PICARD/2.24.0/bin/picard.jar"

for bb in $folder/fastq/*/*_1.fastq.gz; do fq1+=($(echo $(basename ${bb}))) ; done
#for bb in $folder/fastq/*R2*; do fq2+=($(echo $(basename ${bb}))) ; done

for fastq in ${fq1[@]}; do
nameHeader=$(basename ${fastq/_1.fastq.gz//})
FILE=$folder/bam/$nameHeader.bam
#if [ ! -f "$FILE" ]; then 
f1=$folder/fastq/$fastq
f2=$folder/fastq/${fastq/R1/R2}

if [ ! -f $folder/bam_nodup/$nameHeader.bam ]; then

echo "

skewer -o $folder/fastq/$nameHeader $f1 $f2

bwa mem -c 1000 -t 12 $genome $folder/fastq/$nameHeader"-trimmed-pair1.fastq" $folder/fastq/$nameHeader"-trimmed-pair2.fastq" -R '@RG\tID:'$nameHeader'\tLB:'$nameHeader'\tSM:'$nameHeader'\tPL:ILLUMINA'  > $folder/sam/$nameHeader.sam

$picard SamFormatConverter -I $folder/sam/$nameHeader.sam  -O $folder/bam_raw/$nameHeader.bam  -VALIDATION_STRINGENCY SILENT -QUIET true

$picard AddOrReplaceReadGroups -I $folder/bam_raw/$nameHeader.bam -O $folder/bam/$nameHeader.bam -RGID 4 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM 20 -VALIDATION_STRINGENCY SILENT -SO coordinate

$picard MarkDuplicates -I $folder/bam/$nameHeader.bam -O $folder/bam_nodup/$nameHeader.bam -M $folder/QC/$nameHeader.marked_dup_metrics.txt --REMOVE_DUPLICATES true

mv $folder/bam_nodup/$nameHeader.bam $folder/bam/$nameHeader.bam

samtools index $folder/bam/$nameHeader.bam

" >> /home/idib57/idib57798/scripts/sbatch/sbatch.cpt1 ;

sbatch < /home/idib57/idib57798/scripts/sbatch/sbatch.cpt1 ;

sed -i '11,$d' /home/idib57/idib57798/scripts/sbatch/sbatch.cpt1
fi
done


#$picard SortSam  -SO coordinate -I $folder/bam_raw/$nameHeader.bam -O $folder/bam_sort/$nameHeader.bam -VALIDATION_STRINGENCY SILENT -QUIET true --TMP_DIR $folder/tmp/$nameHeader/
