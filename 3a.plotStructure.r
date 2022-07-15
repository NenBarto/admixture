
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(dplyr)
library(plyr)
library(scales)

# Assign the first argument to prefix
prefix="canids9"

# Get individual names in the correct order
labelsSTR<-read.table(paste0(prefix,"/",prefix,".annotation.txt"))

# Name the columns
names(labelsSTR)<-c("ind","pop")

# Add a column with population indices to order the barplots
# Use the order of populations provided as the fourth argument (list separated by commas)
labelsSTR$n<-factor(labelsSTR$pop)
levels(labelsSTR$n)<-c(1:length(levels(labelsSTR$n)))
labelsSTR$n<-as.integer(as.character(labelsSTR$n))

# read in the different admixture output files
minK=2
maxK=8

parseStructure<-function(structurefile){
  input<-readLines(file(structurefile))
  #input<-input[(grep("Inferred ancestry of individuals",input)+1):length(input)]
  #input<-input[2:(grep("Estimated Allele Frequencies",input)-3)]
  foo<-strsplit(sapply(strsplit(input,":"),function(x) x[2])," ")
  foo.mat<-sapply(foo,function(x) x[3:length(foo[[2]])])
  return(t(foo.mat))
}


results<-list()
for(i in minK:maxK){

  #inFile<-paste0(prefix,"/results_str_",prefix,"/str_K",i,"_rep",j,"_f")
  inFile<-paste0(prefix,"/clumpak.",prefix,".files/ClumppIndFile.",i,".output")
  data<-parseStructure(inFile)
  data<-as.data.frame(data)
  data$id<-labelsSTR$ind
  data$pop<-labelsSTR$pop
  dataM<-melt(data,id.vars=c("id","pop"))
  dataM$value<-as.numeric(dataM$value)
  dataM$K<-i
  results[[paste0("K",i)]]<-dataM
}

resultsDF<-do.call("rbind",results)
resultsDF$tool<-"structure"
resultsDF<-resultsDF[!is.na(resultsDF$pop),]

dataM<-resultsDF

cols<-c(brewer.pal(9,"Set1"),"gray")

dataM$pop[dataM$pop=="NG_Singing_Dog"]<-"NGSD"
dataM$pop[dataM$pop=="Wolves"]<-"wolves"
dataM$pop[dataM$pop=="Captive"]<-"captive"

dataM$pop<-factor(dataM$pop,levels=c("wolves","dogs","village_dog","hybrid","backcross","NGSD","captive","alpine","desert","Mallee",""))

pdf(paste0("figures/",prefix,"_structure_K2-8.pdf"),width=30,height=25)
p<-ggplot(dataM,aes(id,value,fill=variable,group=pop))
p<-p+geom_bar(position="stack", stat="identity",colour="NA",width=1)
p<-p+facet_grid(K~pop,space="free_x",scales="free")
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
  data<-read.table(paste0(prefix,"/results_fs_dingos7k_nm/fS_run_K.",i,".meanQ"))
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
pdf(paste0("figures/",prefix,"_fastStructure_K",minK,"-",maxK,".pdf"),width=30,height=10)
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





















colnames(dataM)[3]<-"ancestry"
dataM$ancestry<-gsub("V","",dataM$ancestry)
pdf(paste0("figures/",prefix,"_boxplot2_structure_K",i,".pdf"),width=12,height=5)
p<-ggplot(dataM,aes(pop,value,group=ancestry))
p<-p+geom_boxplot(aes(fill=ancestry))
p<-p+facet_grid(tool~pop,space="free",scales="free")
p<-p+scale_fill_manual(values=cols)
p<-p+theme(axis.text.x=element_text(size=10))
p<-p+
  theme(strip.text.x = element_text(size = 10),
        strip.background = element_rect(color = "NA", size = 1),axis.ticks = element_blank())
p<-p+scale_y_continuous(labels=percent)
p<-p+ylab("Proportion")+xlab("Population")
p
dev.off()


dataMS<-dataM[dataM$ancestry %in% c(1,2),]
dataMS<-dataMS[!dataMS$pop %in% "Mallee",]

pdf(paste0("figures/",prefix,"_boxplot2_dogvsdingo_structure_K",i,".pdf"),width=12,height=5)
p<-ggplot(dataMS,aes(pop,value,group=ancestry))
p<-p+geom_boxplot(aes(fill=ancestry))
p<-p+facet_grid(tool~pop,space="free",scales="free")
p<-p+scale_fill_manual(values=cols)
p<-p+theme(axis.text.x=element_text(size=10))
p<-p+
  theme(strip.text.x = element_text(size = 10),
        strip.background = element_rect(color = "NA", size = 1),axis.ticks = element_blank())
p<-p+scale_y_continuous(labels=percent)
p<-p+ylab("Proportion")+xlab("Population")

p
dev.off()









