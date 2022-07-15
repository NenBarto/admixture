
library(RColorBrewer)
library(ggplot2)
library(reshape2)

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

cols<-c(brewer.pal(9,"Set1"),"gray")
pdf(paste0("figures/",prefix,".pdf"),width=30,height=10)
p<-ggplot(resultsDF,aes(id,value,fill=variable,group=pop))
p<-p+geom_bar(position="stack", stat="identity",colour="NA",width=1)
p<-p+facet_grid(K~pop,space="free_x",scales="free")
p<-p+scale_fill_manual(values=cols)
p<-p+theme(axis.text.x=element_text(angle = -90, hjust = 0,size=2))
p<-p+
  theme(strip.text.x = element_text(size = 8),legend.position="none",panel.spacing=unit(.02, "lines"),
        panel.border = element_rect(color = "NA", fill = NA, size = 1), 
        strip.background = element_rect(color = "NA", size = 1))
p
dev.off()

#get data from other sources

# read in the different admixture output files
minK=4
maxK=4

results<-list()
for(i in minK:maxK){
  data<-read.table(paste0("results_canids9/fS_run_K.",i,".meanQ"))
  data$id<-labels$ind
  data$pop<-labels$pop
  dataM<-melt(data)
  dataM$K<-i
  results[[paste0("K",i)]]<-dataM
}