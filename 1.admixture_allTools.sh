#download vcf parser to pick only the centre SNP
#wget https://raw.githubusercontent.com/CoBiG2/RAD_Tools/6648d1ce1bc1e4c2d2e4256abdefdf53dc079b8c/vcf_parser.py

#recode vcf
vcftools --gzvcf dingo_200_dog_wolf_default_unique_722.3.vcf.gz --maf 0.005 --max-missing 0.60 --recode --out dingo_200_unique
mv dingo_200_unique.recode.vcf dingo_200_unique.vcf
python3 vcf_parser.py --center-snp -vcf dingo_200_unique.vcf
cp annotation_678.txt dingo_200.ind
structure_threader run -i dingos7k_nm.vcf -o ./results_dingos7k_nm -als /share/ScratchGeneral/nenbar/miniconda3/envs/R4.0/bin/alstructure_wrapper.R -K 10 -t 3 --ind dingo_200.ind

vcftools --vcf dingos7k_nm.vcf --maf 0.005 --max-missing 0.60 --recode --out dingo7k_nm_unique
mv dingo7k_nm_unique.recode.vcf dingo7k_nm_unique.vcf
python3 vcf_parser.py --center-snp -vcf dingo7k_nm_unique.vcf
cp dingos7k_nm.annotation.txt dingo7k_nm.ind


########## STRUCTURE


ncores=24
for i in {8..1}; do
  strucLine="structure_threader run -i dingos7k_nm.str -o ./results_str_dingos7k_nm -st /share/ScratchGeneral/nenbar/miniconda3/envs/R4.0/bin/structure -Klist $i -t 12 -R 10 --extra_opts \"-D $i\" --ind dingo7k_nm.ind"
  qsub -b y -cwd -j y -N struc1.$i -R y -pe smp $ncores -V $strucLine
done


for i in {8..2}; do
  strucLine="structure_threader run -i canids9.str -o ./results_str_canids9 -st /share/ScratchGeneral/nenbar/miniconda3/envs/R4.0/bin/structure -Klist $i -t 12 -R 10 --extra_opts \"-D $i\" --ind canids9.populations.str.txt"
  qsub -b y -cwd -j y -N struc$i -R y -pe smp $ncores -V $strucLine
done


########## FASTSTRUCTURE

ncores=24
for i in {3..5}; do
  strucLine="structure_threader run -i dingos7k_nm.str -o ./results_fs_dingos7k_nm -fs /share/ScratchGeneral/nenbar/miniconda3/envs/R4.0/bin/fastStructure -Klist $i -t 12 --ind dingo7k_nm.ind"
  qsub -b y -cwd -j y -N fsD.$i -R y -pe smp $ncores -V $strucLine
done


#run for all canids, but 3 first
for i in {3..7}; do
  strucLine="structure_threader run -i canids9.str -o ./results_fs_canids9 -fs /share/ScratchGeneral/nenbar/miniconda3/envs/R4.0/bin/fastStructure -Klist $i -t 12 --ind canids9.populations.str.txt"
  qsub -b y -cwd -j y -N fsC.$i -R y -pe smp $ncores -V $strucLine
done


########## ALstructure 

ncores=24
for i in {3..5}; do
  strucLine="structure_threader run -i dingo7k_nm_uniqueCenterSNP.vcf -o ./results_als_dingos7k_nm -als /share/ScratchGeneral/nenbar/miniconda3/envs/R4.0/bin/alstructure_wrapper.R -Klist $i -t 12 --ind dingo7k_nm.ind"
  qsub -b y -cwd -j y -N strucALSd.$i -R y -pe smp $ncores -V $strucLine
done;

ncores=24
for i in {3..7}; do
  #strucLine="structure_threader run -i canids9.vcf -o ./results_als_canids9 -als /share/ScratchGeneral/nenbar/miniconda3/envs/R4.0/bin/alstructure_wrapper.R -K 7 -t 12 --ind canids9.annotation.txt"
  strucLine="structure_threader run -i canids9.recode.vcf -o ./results_als_canids9 -als /share/ScratchGeneral/nenbar/miniconda3/envs/R4.0/bin/alstructure_wrapper.R -Klist $i -t 12 --ind canids9.annotation.txt"
  qsub -b y -cwd -j y -N strucALSc.$i -R y -pe smp $ncores -V $strucLine
done;

structure_threader run -i dingo_200_uniqueCenterSNP.vcf -o ./results_dingo_200_unique -als /share/ScratchGeneral/nenbar/miniconda3/envs/R4.0/bin/alstructure_wrapper.R -K 5 -t 3 --ind dingo_200.ind



ncores=32
strucLine="structure_threader run -i dingos7k_nm.str -o ./results_dingos7k_nm -st /share/ScratchGeneral/nenbar/miniconda3/envs/R4.0/bin/structure -Klist 3 -t 16 -R 10 --ind dingo7k_nm.ind"
qsub -b y -cwd -j y -N struc3 -R y -pe smp $ncores -V $strucLine

structure_threader run -i dingo_200_uniqueCenterSNP.vcf -o ./results_dingo_200_unique -als /share/ScratchGeneral/nenbar/miniconda3/envs/R4.0/bin/alstructure_wrapper.R -K 10 -t 3 --ind dingo_200.ind

dingos7k_nm.annotation.txt

structure_threader run -i dingo_200_uniqueCenterSNP.vcf -o ./results_structure_dingo_200_unique -st /share/ScratchGeneral/nenbar/miniconda3/envs/R4.0/bin/structure -K 10 -t 3 --ind dingo_200.ind




for i in {1..10}
do
 admixture --cv dingo_200_unique.bed $i > log${i}.out &
done
grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > dingo_200_unique.cv.error
grep "CV" *out | awk '{print $3,$4}' | cut -c 4,7-20 > dingo_200_unique.cv.error
awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > dingo_200_unique.cv.error
awk '{split($1,name,"."); print $1,name[2]}' dingo_200_unique.nosex > dingo_200_unique.list



for i in {11..15}
do
 admixture --cv dingos7k_nm.bed $i > log${i}.out &
done
grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > dingos7k_nm.cv.error
grep "CV" *out | awk '{print $3,$4}' | cut -c 4,7-20 > dingos7k_nm.cv.error
awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > dingos7k_nm.cv.error
awk '{split($1,name,"."); print $1,name[2]}' dingos7k_nm.nosex > dingos7k_nm.list


for i in {1..10}
do
 admixture --cv canids9.bed $i > log${i}.out &
done
grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > canids9.cv.error
grep "CV" *out | awk '{print $3,$4}' | cut -c 4,7-20 > canids9.cv.error
awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > canids9.cv.error
awk '{split($1,name,"."); print $1,name[2]}' canids9.nosex > canids9.list



