#!/usr/bin/env Rscript
library(ggplot2)
library(ggpubr)
library(readr)
library(plyr)
library(dplyr)
library(parallel)
library(stringr)
library(pbapply)
library(ggpointdensity)
library(reshape2)

####fig1b-d
heatmap <- function(df_long,start,mid,mid_value) {
  heatmap_fig<-ggplot(df_long, aes(x = Col, y = Samples, fill = value)) + geom_tile() + 
    scale_fill_gradientn(colors = c("blue","white", "red"), breaks = c(start,mid, 1),         
                         labels = c(start,mid, "1"),limits = c(start, 1),          
                         values = c(0,mid_value,1),na.value = "grey80") + theme_minimal() +   
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5,size=8),
          legend.direction = "horizontal",legend.position = "top",legend.margin = margin(t = 10, r = 5, b = 0.01, l = 5),
          legend.key.width = unit(0.6, "cm"),legend.key.height = unit(0.25, "cm"),
          axis.title = element_text(size=8),axis.text = element_text(size=5),
          legend.text = element_text(size=8),legend.title = element_text(size=8)) + 
    coord_fixed() +  labs(x = "Samples", y = "Samples")+
    ggtitle("Pearson's correlation coefficient for the\n number of reads mapped to each genes")
  return(heatmap_fig)
}

load("/mnt/data6/disk/huangyanying/update_M_heart/MR_counts_heatmap.Rdata")##add data
heatmap_M_new<-heatmap(df_long_M,0.4,0.75,0.6)
heatmap_R_new<-heatmap(df_long_R,0.5,0.85,0.7)
heatmap_C_new<-heatmap(df_long_C,0.3,0.65,0.5)
fig1b_d<-ggarrange(heatmap_M_new,heatmap_R_new,heatmap_C_new,ncol=3)


###fig1b inside 
load("/mnt/data6/disk/huangyanying/mouse_plot/04_pea/02conbine_pea/new_all_M.Rdata")
filter_2reps <-data_M %>% dplyr::select("genename",`M-heart-D1`, `M-heart-D2`) %>%na.omit()
filter_2reps$`M-heart-D1` <- log10(filter_2reps$`M-heart-D1` + 1)  
filter_2reps$`M-heart-D2` <- log10(filter_2reps$`M-heart-D2` + 1)
p_value=round(cor.test(filter_2reps$`M-heart-D1`,filter_2reps$`M-heart-D2`,method="p")$estimate,3)

fig1b_inside<-ggplot(filter_2reps,aes(x=`M-heart-D1`,y=`M-heart-D2`))+geom_point(size=0.01)+  
  scale_x_continuous(breaks = c(0, 1, 3, 5, 7),  
                     labels = c("0", "10¹", "10³", "10⁵", "10⁷"),limits = c(0, 7))+    
  scale_y_continuous(breaks = c(0, 1, 3, 5, 7), 
                     labels = c("0", "10¹", "10³", "10⁵", "10⁷"),limits = c(0, 7))+ 
  labs(title = "",x = "Number of reads", y = "Number of reads")+
  theme(plot.background = element_blank(),panel.background = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
        plot.title = element_text(hjust = 0.5),axis.ticks = element_line(linewidth = 0.3),
        panel.grid = element_blank(),axis.ticks.length = unit(0.04, "cm"),
        axis.title = element_text(size=5.5),axis.text = element_text(size=5.5),
        legend.text = element_text(size=5.5),legend.title = element_text(size=5.5))


###fig 1e
###gene
##横坐标为average numbers of reads of the gene，纵坐标为两个样本之间每个基因的FSR的pearson相关性值
pearson1<-read.table("/mnt/data6/disk/huangyanying/mouse_plot/01_rawdata/filtersite/C/C-brain-N_union.txt",header = T)
density1<-read.table("/mnt/data6/disk/huangyanying/mouse_plot/FSR02/filter_density/C-brain-N_freq1.txt",header = T)
colnames(density1)=c("genename","freq","genelength","exon","intron","density")
pearson2=merge(pearson1,density1,by="genename",all = T)
reps_pearson2<-pearson2[order(pearson2$density),]
reps_pearson2_1<-na.omit(reps_pearson2)
fig1e<-ggplot(reps_pearson2_1,aes(density,rep1_rep2))+
  scale_x_log10(limits=c(1,1000),breaks = c(1,10,100,1000),labels=as.character(c(1,10,100,1000)))+
  scale_y_continuous(limits = c(0,1),breaks = c(0,0.25,0.50,0.75,1),labels=as.character(c(0,0.25,0.50,0.75,1)))+
  labs(x="Reads coverage of the genes",
       y="Pearson’s correlation of RT stop counts between \nbiological repeats of NAI-N3 treated Cavia brain")+
  theme_bw()+geom_pointdensity(size=0.8,alpha=0.6,adjust = 0.5)+
  geom_point(data = reps_pearson2_1[4000, ],color = "red", size = 1,alpha=1) + #例子标红（U5 snRNA）
  theme(plot.title = element_text(hjust = 0.5),panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),axis.title = element_text(size=8),
        axis.text = element_text(size=8),legend.text = element_text(size=8),
        legend.title = element_text(size=8),legend.position="none")


###fig 1f
#######C_brain_U5_RT_stopsites_reps
C_brain_U5_rep1<-read.table("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig/C_brain_U5/C_brain_U5_rep1.txt",
                            header = F,stringsAsFactors = F)
C_brain_U5_rep2<-read.table("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig/C_brain_U5/C_brain_U5_rep2.txt",
                            header = F,stringsAsFactors = F)
U5_add<-function(RT_site,rep_name){
  rep1<-as.data.frame(table(RT_site))
  rep1$RT_site<-as.numeric(as.character(rep1$RT_site))
  add1<-data.frame(RT_site=setdiff(1:116,rep1$RT_site),Freq=0)
  rep1_add<-rbind(rep1,add1)
  rep1_add<-rep1_add%>%dplyr::mutate(type=rep_name)
  return(rep1_add)
}
U5_C_brain_rep1_add<-U5_add(C_brain_U5_rep1$V2,"rep1")
U5_C_brain_rep2_add<-U5_add(C_brain_U5_rep2$V2,"rep2")
U5_C_brain_rep1_rep2<-rbind(U5_C_brain_rep1_add,U5_C_brain_rep2_add)

fig_1f<-ggplot(U5_C_brain_rep1_rep2,aes(x=RT_site,y=Freq,fill=type))+geom_bar(stat = "identity")+
  facet_grid(type~.) +ylab("RT stops")+xlab("Nucleotide positions(nt) of U5 snRNA")+
  scale_x_continuous(breaks = c(1,20,40,60,80,100,116))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.title = element_text(size=8),axis.text = element_text(size=8),
        legend.text = element_text(size=8),legend.title = element_text(size=8),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size=8),legend.position = "none")



