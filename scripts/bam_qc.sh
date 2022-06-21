module load intel java fastqc samtools R qualimap star python/3.6.5 gcc/9.2.0 sortmerna skewer/0.2.2

folder=/slgpfs/projects/idib57/GA_002_21_Radon
gtf=/slgpfs/projects/idib57/data/genome/Rattus/Rattus_norvegicus.Rnor_6.0.104.gtf


bamfiles=(BAM/RnD1-134_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD1-162_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD1-195_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD1-207_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD12-135_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD12-145_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD1-214_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD12-169_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD12-178_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD12-184_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD12-185-1_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD12-185-2_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD12-205_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD12-212_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD12-229_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD12-231_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD12-24_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD1-228_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD12-49_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD12-51_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD1-36_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD1-81_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD3-11_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD3-168_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD3-203_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD3-212_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD3-232_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD3-233_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD3-66_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD3-92_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD6-132_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD6-152_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD6-16_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD6-179_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnD6-204_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnPr109_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnPr140_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnPr146_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnPr214_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnSE38_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnSE3_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnSE80_R1_val_1.fqAligned.sortedByCoord.out.bam BAM/RnSE90_R1_val_1.fqAligned.sortedByCoord.out.bam)
for bam in ${bamfiles[@]}; do

nameHeader=$(basename ${bam/_R1_val_1.fqAligned.sortedByCoord.out.bam/})

echo " 
samtools index $folder/$bam

samtools idxstats $folder/$bam > $folder/QC/$nameHeader.idxstats.txt

fastqc $folder/$bam -o $folder/QC/fastqc/

qualimap bamqc -nt 48 -bam $folder/$bam -c -gd RAT - RGSC6 -gff $gtf -outdir $folder/QC/$nameHeader -sd --java-mem-size=8G

geneBody_coverage.py -i $folder/$bam -r /slgpfs/projects/idib57/data/genome/Rattus/RGSC6.0_rn6.exome.nochr.bed -o $folder/QC/RSeQC/$nameHeader
read_distribution.py -i $folder/$bam -r /slgpfs/projects/idib57/data/genome/Rattus/RGSC6.0_rn6.exome.nochr.bed > $folder/QC/RSeQC/rd_$nameHeader
" >> /slgpfs/projects/idib57/GA_002_21_Radon/sbatch.12h.cpt1 ;

sbatch /slgpfs/projects/idib57/GA_002_21_Radon/sbatch.12h.cpt1 ;

sed -i '10,$d' /slgpfs/projects/idib57/GA_002_21_Radon/sbatch.12h.cpt1 ;

done
