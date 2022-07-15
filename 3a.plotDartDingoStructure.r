
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(dplyr)
library(plyr)
library(scales)

# Assign the first argument to prefix
prefix="dingos7k"

# Get individual names in the correct order
labels<-read.table(paste0(prefix,"/",prefix,".populations.str.txt"))

# Name the columns
names(labels)<-c("ind","pop")

# Add a column with population indices to order the barplots
# Use the order of populations provided as the fourth argument (list separated by commas)
labels$n<-factor(labels$pop)
levels(labels$n)<-c(1:length(levels(labels$n)))
labels$n<-as.integer(as.character(labels$n))

# read in the different admixture output files
minK=3
maxK=3

totalResults<-list()


###### Structure 

parseStructure<-function(structurefile){
  input<-readLines(file(structurefile))
  input<-input[(grep("Inferred ancestry of individuals",input)+1):length(input)]
  input<-input[2:(grep("Estimated Allele Frequencies",input)-3)]
  foo<-strsplit(sapply(strsplit(input,":"),function(x) x[2])," ")
  foo.mat<-sapply(foo,function(x) x[3:length(foo[[2]])])
  return(t(foo.mat))
}


results<-list()
for(i in 3){
  inFile<-paste0("dingos7k_nm/results_dingos7k_nm/str_K",i,"_rep16_f")
  data<-parseStructure(inFile)
  data<-as.data.frame(data)
  data$id<-labels$ind
  data$pop<-labels$pop
  dataM<-melt(data,id.vars=c("id","pop"))
  dataM$value<-as.numeric(dataM$value)
  dataM$K<-i
  results[[paste0("K",i)]]<-dataM
}

resultsDF<-do.call("rbind",results)
resultsDF$tool<-"structure"
resultsDF<-resultsDF[!is.na(resultsDF$pop),]
totalResults[["structure"]]<-resultsDF



cols<-c(brewer.pal(4,"Set1"),"gray")

dataM<-totalResults[["structure"]]

dataM$pop<-factor(dataM$pop,levels=c("dog","alpine","desert","Mallee",""))

pdf(paste0("figures/",prefix,"_structure_dartDogDingo_K",i,".pdf"),width=30,height=5)
p<-ggplot(dataM,aes(id,value,fill=variable,group=pop))
p<-p+geom_bar(position="stack", stat="identity",colour="NA",width=1)
p<-p+facet_grid(tool~pop,space="free",scales="free")
p<-p+scale_fill_manual(values=cols)
p<-p+theme(axis.text.x=element_text(angle = -90, hjust = 0,size=2))
p<-p+
  theme(strip.text.x = element_text(size = 8),legend.position="none",panel.spacing=unit(.02, "lines"),
        panel.border = element_rect(color = "NA", fill = NA, size = 1), 
        strip.background = element_rect(color = "NA", size = 1))
p<-p+scale_y_continuous(labels=percent)
p<-p+ylab("Proportion")+xlab("Population")
p
dev.off()



colnames(dataM)[3]<-"ancestry"
dataM$ancestry<-gsub("V","",dataM$ancestry)
pdf(paste0("figures/",prefix,"_boxplot2_structure_dartDogDingo_K",i,".pdf"),width=12,height=5)
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


dataMS<-dataM[dataM$ancestry %in% c(2,3),]
dataMS<-dataMS[!dataMS$pop %in% "Mallee",]

pdf(paste0("figures/",prefix,"_boxplot2_dogvsdingo_structure_dartDogDingo_K",i,".pdf"),width=12,height=5)
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











