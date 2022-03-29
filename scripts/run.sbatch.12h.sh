input=$1

echo "$input

" >> /slgpfs/projects/idib57/GA_002_21_Radon/sbatch.12h.cpt1 ;

sbatch < /slgpfs/projects/idib57/GA_002_21_Radon/sbatch.12h.cpt1 ;

sed -i '10,$d' /slgpfs/projects/idib57/GA_002_21_Radon/sbatch.12h.cpt1
