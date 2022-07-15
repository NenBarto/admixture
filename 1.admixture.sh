#mkdir ~/str_analyses
#cd ~/str_analyses
#wget https://raw.githubusercontent.com/radcamp/radcamp.github.io/master/Lisbon2020/05_POPULATION_STRUCTURE_files/oyster.vcf.gz
#wget https://raw.githubusercontent.com/radcamp/radcamp.github.io/master/Lisbon2020/05_POPULATION_STRUCTURE_files/oyster.indfile
#
#conda install -c bioconda vcftools
#centos6.10/pethum/vcftools/gcc-4.4.6/0.1.10
#vcftools --gzvcf oyster.vcf.gz --maf 0.005 --max-missing 0.60 --recode --out oysterMAF005MM60
#mv oysterMAF005MM60.recode.vcf oysterMAF005MM60.vcf
#
#wget https://raw.githubusercontent.com/CoBiG2/RAD_Tools/6648d1ce1bc1e4c2d2e4256abdefdf53dc079b8c/vcf_parser.py
#python3 vcf_parser.py --center-snp -vcf oysterMAF005MM60.vcf
#
#structure_threader run -i oysterMAF005MM60CenterSNP.vcf -o ./results_oysterMAF005MM60CenterSNP -als /share/ScratchGeneral/nenbar/miniconda3/envs/R4.0/bin/alstructure_wrapper.R -K 10 -t 3 --ind oyster.indfile

#recode vcf
vcftools --gzvcf dingo_200_dog_wolf_default_unique_722.3.vcf.gz --maf 0.005 --max-missing 0.60 --recode --out dingo_200_unique
mv dingo_200_unique.recode.vcf dingo_200_unique.vcf
python3 vcf_parser.py --center-snp -vcf dingo_200_unique.vcf
cp annotation_678.txt dingo_200.ind
structure_threader run -i dingo_200_uniqueCenterSNP.vcf -o ./results_dingo_200_unique -als /share/ScratchGeneral/nenbar/miniconda3/envs/R4.0/bin/alstructure_wrapper.R -K 10 -t 3 --ind dingo_200.ind


#run structure analysis

library(RColorBrewer)
library(ggplot2)
cols<-brewer.pal(7,"Set1")

data<-read.csv("./results_dingo_200_unique/alstr_K6")
df<-read.table("dingo_200.ind")
colnames(df)<-c("sample","pop","idx")
df$type<-"dingo"
df$type[grepl("dog",df$pop)]<-"dog"
df$type[grepl("wolf|Fox",df$pop)]<-"wolf"
df$type[grepl("mixed|hybrid|backcross",df$pop)]<-"mixed"
df$type[grepl("vill",df$pop)]<-"village"
df$type[grepl("New_Guinea_SD",df$pop)]<-"New_Guinea_SD"
df$type[grepl("Mallee",df$pop)]<-"Mallee"

data<-data[,-1]
pdf("ALStructure_678.pdf",width=12,height=8)
pca<-princomp(data)
#plot(pca$loading,pch=19, cex=2,col=cols)
plot(pca$scores,pch=19, cex=2)
text(pca$scores, df$pop,pos = 1)
dev.off()

df$score1<-pca$scores[,1]
df$score2<-pca$scores[,2]


pdf("ALStructure_678_ggplot.pdf",width=12,height=8)
p<-ggplot(df,aes(score1,score2,group=type))
p<-p+geom_point(aes(color=type))
p<-p+scale_color_manual(values=cols)
p
dev.off()


structure_threader run -i dingo_200_uniqueCenterSNP.vcf -o ./results_structure_dingo_200_unique -st /share/ScratchGeneral/nenbar/miniconda3/envs/R4.0/bin/structure -K 10 -t 3 --ind dingo_200.ind




for i in {3..10}
do
 admixture --cv dingo_200_unique.bed $i > log${i}.out &
done
grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > dingo_200_unique.cv.error
grep "CV" *out | awk '{print $3,$4}' | cut -c 4,7-20 > dingo_200_unique.cv.error
awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > dingo_200_unique.cv.error
awk '{split($1,name,"."); print $1,name[2]}' dingo_200_unique.nosex > dingo_200_unique.list

