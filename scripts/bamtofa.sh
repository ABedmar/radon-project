module load intel java samtools

for bam in BAM/*.bam; do

nameHeader=$(basename ${bam/_R1_val_1.fqAligned.sortedByCoord.out.bam/})

echo " 

samtools fastq -1 fastq1/$nameHeader"_R1.fq" -2 fastq1/$nameHeader"_R2.fq" -0 fastq1/$nameHeader"_of" -s fastq1/$nameHeader"_sing" -n -F 0x900 $bam

" >> /slgpfs/projects/idib57/GA_002_21_Radon/sbatch.12h.cpt1 ;
sbatch /slgpfs/projects/idib57/GA_002_21_Radon/sbatch.12h.cpt1 ;
sed -i '10,$d' /slgpfs/projects/idib57/GA_002_21_Radon/sbatch.12h.cpt1 ;
done