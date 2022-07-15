
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(dplyr)
library(plyr)


# Assign the first argument to prefix
prefix="canids9"

# Get individual names in the correct order
labels<-read.table(paste0(prefix,".annotation.txt"))

# Name the columns
names(labels)<-c("ind","pop")

# Add a column with population indices to order the barplots
# Use the order of populations provided as the fourth argument (list separated by commas)
labels$n<-factor(labels$pop)
levels(labels$n)<-c(1:length(levels(labels$n)))
labels$n<-as.integer(as.character(labels$n))

# read in the different admixture output files
minK=4
maxK=4

totalResults<-list()


###### Admixture 
results<-list()
for(i in minK:maxK){
  data<-read.table(paste0("admixture/",prefix,".",i,".Q"))
  data$id<-labels$ind
  data$pop<-labels$pop
  dataM<-melt(data)
  dataM$K<-i
  results[[paste0("K",i)]]<-dataM
}

resultsDF<-do.call("rbind",results)
resultsDF$tool<-"Admixture"

totalResults[["Admixture"]]<-resultsDF
#get data from other sources

###### fastStructure 
labelsSTR<-read.table(paste0(prefix,".populations.str.txt"))
colnames(labelsSTR)<-c("ind","pop")
results<-list()
for(i in minK:maxK){
  data<-read.table(paste0("results_canids9/fS_run_K.",i,".meanP"))
  data$id<-labelsSTR$ind
  data$pop<-labelsSTR$pop
  dataM<-melt(data)
  dataM$K<-i
  results[[paste0("K",i)]]<-dataM
}

resultsDF<-do.call("rbind",results)
resultsDF$tool<-"fastStructure"
totalResults[["fastStructure"]]<-resultsDF

###### ALStructure 
labelsALStructure<-read.table(paste0(prefix,".annotation.txt"))
colnames(labelsALStructure)<-c("ind","pop")

results<-list()
for(i in minK:maxK){
  data<-read.csv(paste0("results_canids9/alstr_K",i),header=T)
  data<-data[,-1]
  data$id<-labelsALStructure$ind
  data$pop<-labelsALStructure$pop
  dataM<-melt(data)
  dataM$K<-i
  results[[paste0("K",i)]]<-dataM
}

resultsDF<-do.call("rbind",results)
resultsDF$tool<-"ALStructure"
totalResults[["ALStructure"]]<-resultsDF



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
for(i in minK:maxK){
  inFile<-paste0("results_canids9/str_K",i,"_rep15_f")
  data<-parseStructure(inFile)
  data<-as.data.frame(data)
  data$id<-labelsSTR$ind
  data$pop<-labelsSTR$pop
  dataM<-melt(data,id.vars=c("id","pop"))
  dataM$K<-i
  results[[paste0("K",i)]]<-dataM
}

resultsDF<-do.call("rbind",results)
resultsDF$tool<-"structure"
totalResults[["structure"]]<-resultsDF

df<-do.call("rbind",totalResults)

annotation<-read.table("DCan22combined_annotation.txt",sep="\t",header=T)

df$newID<-df$id

df$id<-mapvalues(df$newID,annotation$targetid,annotation$id)
df<-df[df$id %in% data$id,]

#now plot the data


cols<-c(brewer.pal(9,"Set1"),"gray")
pdf(paste0("figures/",prefix,"_allTools.pdf"),width=30,height=10)
p<-ggplot(df,aes(id,value,fill=variable,group=pop))
p<-p+geom_bar(position="stack", stat="identity",colour="NA",width=1)
p<-p+facet_grid(tool~pop,space="free",scales="free")
p<-p+scale_fill_manual(values=cols)
p<-p+theme(axis.text.x=element_text(angle = -90, hjust = 0,size=2))
p<-p+
  theme(strip.text.x = element_text(size = 8),legend.position="none",panel.spacing=unit(.02, "lines"),
        panel.border = element_rect(color = "NA", fill = NA, size = 1), 
        strip.background = element_rect(color = "NA", size = 1))
p
dev.off()



