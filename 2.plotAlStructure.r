
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(scales)
library(forcats)
library(plyr)



# Assign the first argument to prefix
prefix="canids9"

# Get individual names in the correct order
labels<-read.table(paste0(prefix,"/",prefix,".annotation.txt"))
labels$V2[labels$V2=="Captive"]="captive"
labels$V2[labels$V2=="Wolves"]="wolves"
labels$V2[labels$V2=="NG_Singing_Dog"]="NG_singing_dog"

# Name the columns
names(labels)<-c("ind","pop")

labels$ind[labels$ind=="AmHairlessTerr01"]<-"AmericanHairlessTerrier01"
labels$ind[labels$ind=="AustCattleDog01"]<-"AustralianCattleDog01"
labels$ind[labels$ind=="AustCattleDog02"]<-"AustralianCattleDog02"
labels$ind[labels$ind=="BerneseMntnDog01"]<-"BerneseMountainDog01"
labels$ind[labels$ind=="BlackRussianTerr01"]<-"BlackRussianTerrier01"
labels$ind[labels$ind=="Ovcharka01"]<-"CaucasianOvcharka01"
labels$ind[labels$ind=="ChineseIndiDog01"]<-"ChineseIndigenousDog01"
labels$ind[labels$ind=="CockerSpanielAm01"]<-"CockerSpanielAmerican01"
labels$ind[labels$ind=="EnglishSpringer01"]<-"EnglishSpringerSpaniel01"
labels$ind[labels$ind=="Entlebucherhund01"]<-"EntlebucherSennenhund01"
labels$ind[labels$ind=="FlatcoatRetvr01"]<-"FlatcoatedRetriever01"
labels$ind[labels$ind=="IndiDogNigeria01"]<-"IndigenousDogNigeria01"
labels$ind[labels$ind=="IndiDogVietnam01"]<-"IndigenousDogVietnam01"
labels$ind[labels$ind=="IrishWaterSpanl01"]<-"IrishWaterSpaniel01"
labels$ind[labels$ind=="IstrianShorthair01"]<-"IstrianShorthairedHound01"
labels$ind[labels$ind=="JackRussell02"]<-"JackRussellTerrier02"
labels$ind[labels$ind=="Labrador13"]<-"LabradorRetriever13"
labels$ind[labels$ind=="Lagotto02"]<-"LagottoRomagnolo02"
labels$ind[labels$ind=="NewGuineaSD01"]<-"NewGuineaSingingDog01"
labels$ind[labels$ind=="NovaScotiaRetrv01"]<-"NovaScotiaDuckTollingRetriever01"
labels$ind[labels$ind=="PortuguesePod01"]<-"PortuguesePodengo01"
labels$ind[labels$ind=="PortugueseWater01"]<-"PortugueseWaterDog01"
labels$ind[labels$ind=="RhodesianRidgebk01"]<-"RhodesianRidgeback01"
labels$ind[labels$ind=="ScottishDeerhnd01"]<-"ScottishDeerhound01"
labels$ind[labels$ind=="VillDog_Aust01"]<-"VillDog_Australia01"
labels$ind<-gsub("CFA_","CFA.",labels$ind)

# Add a column with population indices to order the barplots
# Use the order of populations provided as the fourth argument (list separated by commas)
labels$n<-factor(labels$pop)
levels(labels$n)<-c(1:length(levels(labels$n)))
labels$n<-as.integer(as.character(labels$n))

# read in the different admixture output files
minK=3
maxK=7

results<-list()
oldLabels<-read.table("backup/canids9.annotation.txt",header=T)


for(i in minK:maxK){
  data<-read.csv(paste0(prefix,"/","results_als_canids9/alstr_K",i),header=T)
  data<-data[,-1]
  data$id<-oldLabels$anno.id
  data$pop<-oldLabels$anno.pop
  dataM<-melt(data)
  dataM$K<-i
  results[[paste0("K",i)]]<-dataM
}

resultsDF<-do.call("rbind",results)

#first convert the numbers to Peter's ids
anno<-read.table("../annotation/DCan22combined_annotation.txt",header=T,sep="\t")
resultsDF$newID<-resultsDF$id
resultsDF$newID<-mapvalues(resultsDF$id,anno$targetid,anno$New_6734_id)
resultsDFS<-resultsDF[resultsDF$newID %in% labels$ind,]
resultsDFS$popNew<-mapvalues(resultsDFS$newID,labels$ind,labels$pop)


resultsDFS$popNew<-factor(resultsDFS$popNew,levels=c("wolves","dogs","village_dog","NG_singing_dog","hybrid","backcross","captive","alpine","desert","Mallee"))

cols<-c(brewer.pal(9,"Set1"),"gray")
pdf(paste0("figures/",prefix,"_ALStructure_K3-7.pdf"),width=30,height=20)
p<-ggplot(resultsDFS,aes(newID,value,fill=variable,group=popNew))
p<-p+geom_bar(position="stack", stat="identity",colour="NA",width=1)
p<-p+facet_grid(K~popNew,space="free_x",scales="free")
p<-p+scale_fill_manual(values=cols)
p<-p+theme(axis.text.x=element_text(angle = -90, hjust = 0,size=1))
p<-p+theme(axis.text.y=element_text(size=24))
p<-p+theme(axis.title.x=element_text(size=24))
p<-p+theme(axis.title.y=element_text(size=24))
p<-p+
  theme(strip.text.x = element_text(size = 12,angle = 90),legend.position="none",panel.spacing=unit(.02, "lines"),
        strip.text.y = element_text(size = 24),panel.border = element_rect(color = "NA", fill = NA, size = 1), 
        strip.background = element_rect(color = "NA", size = 24))
p<-p+scale_y_continuous(labels=percent)
p<-p+ylab("Proportion")+xlab("Population")
p
dev.off()



# Assign the first argument to prefix
prefix="dingos7k_nm"

# Get individual names in the correct order
labels<-read.table(paste0(prefix,"/",prefix,".annotation.txt"))

# Name the columns
names(labels)<-c("ind","pop")

# Add a column with population indices to order the barplots
# Use the order of populations provided as the fourth argument (list separated by commas)
labels$n<-factor(labels$pop)
levels(labels$n)<-c(1:length(levels(labels$n)))
labels$n<-as.integer(as.character(labels$n))

# read in the different admixture output files
minK=3
maxK=5

results<-list()
for(i in minK:maxK){
  data<-read.csv(paste0(prefix,"/","results_als_dingos7k_nm/alstr_K",i),header=T)
  data<-data[,-1]
  data$id<-labels$ind
  data$pop<-labels$pop
  dataM<-melt(data)
  dataM$K<-i
  results[[paste0("K",i)]]<-dataM
}


resultsDF<-do.call("rbind",results)
resultsDF$pop<-factor(resultsDF$pop,levels=c("dog","alpine","desert","Mallee"))

#resultsDF$variable <- relevel(resultsDF$variable, c('V1'))
resultsDF$NewVariable<-resultsDF$variable
#resultsDF$NewVariable[resultsDF$variable=="V1"]<-"V3"
#resultsDF$NewVariable[resultsDF$variable=="V2"]<-"V2"
#resultsDF$NewVariable[resultsDF$variable=="V3"]<-"V1"
#resultsDF$NewVariable<-factor(resultsDF$NewVariable,levels=c("V2","V1","V3"))


cols<-c(brewer.pal(6,"Set1"))
pdf(paste0("figures/",prefix,"_ALStructure_K",minK,"-",maxK,".pdf"),width=30,height=10)
p<-ggplot(resultsDF,aes(id,value,fill=NewVariable,group=pop,order = -as.numeric(NewVariable)))
p<-p+geom_bar(stat="identity",colour="NA",width=1,position = position_fill(reverse = TRUE))
p<-p+facet_grid(K~pop,space="free_x",scales="free")
p<-p+scale_fill_manual(values=cols)
p<-p+theme(axis.text.x=element_text(angle = -90, hjust = 0,size=1))
p<-p+theme(axis.text.y=element_text(size=24))
p<-p+theme(axis.title.x=element_text(size=24))
p<-p+theme(axis.title.y=element_text(size=24))
p<-p+
  theme(strip.text.x = element_text(size = 24),legend.position="none",panel.spacing=unit(.02, "lines"),
        strip.text.y = element_text(size = 24),panel.border = element_rect(color = "NA", fill = NA, size = 1), 
        strip.background = element_rect(color = "NA", size = 24))
p<-p+scale_y_continuous(labels=percent)
p<-p+ylab("Proportion")+xlab("Population")
p
dev.off()


colnames(resultsDF)[3]<-"ancestry"
resultsDF$ancestry<-gsub("V","",resultsDF$ancestry)
pdf(paste0("../figures/",prefix,"_boxplot2_ALStructure_dartDogDingo_K",i,".pdf"),width=12,height=5)
p<-ggplot(resultsDF,aes(pop,value,group=ancestry))
p<-p+geom_boxplot(aes(fill=ancestry))
p<-p+facet_grid(~pop,space="free",scales="free")
p<-p+scale_fill_manual(values=cols)
p<-p+theme(axis.text.x=element_text(size=10))
p<-p+
  theme(strip.text.x = element_text(size = 10),
        strip.background = element_rect(color = "NA", size = 1),axis.ticks = element_blank())
p<-p+scale_y_continuous(labels=percent)
p<-p+ylab("Proportion")+xlab("Population")
p
dev.off()

