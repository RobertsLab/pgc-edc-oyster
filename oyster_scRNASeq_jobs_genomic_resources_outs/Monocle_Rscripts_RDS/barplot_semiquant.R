library(data.table)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(readr)
library(grid)
#read in a single cov file
semiq<-fread("~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/semiquant.data.txt")
head(semiq)
#read in stdev's
#seems like this is better so I don't have to learn how to do this (since I'll never do it again)

#summarize mean and sd of each life-stage/gene combo for plotting
semiq2<- semiq %>% group_by(life_stage,gene) %>% summarize (mean_ex=mean(reletive_exp), sd_ex=sd(reletive_exp))
head(semiq2)

#barplot with error bars, thanks to https://statdoe.com/barplot-for-two-factors-in-r/
ggplot(semiq2, aes(x = factor(gene), y = mean_ex, fill = life_stage)) + 
geom_bar(stat = "identity", position = "dodge") +
scale_fill_manual(values = c("#0a4158", "#4b8378","#ff9636","#e4d7d0")) +
geom_errorbar(aes(ymin=mean_ex-sd_ex, ymax=mean_ex+sd_ex), position = position_dodge(0.9), width = 0.25, colour="black") +
theme_minimal() +
theme(legend.position="bottom", axis.text.y = element_text(size = 20)) 

ggsave("semiq_barplot.png",
       dpi = 750,
       height = 5,
       width = 6,
       bg = "white")
  