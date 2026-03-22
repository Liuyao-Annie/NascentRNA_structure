#!/usr/bin/env Rscript
library(ggplot2)
library(ggpubr)
library(readr)
library(plyr)
library(dplyr)
library(parallel)
library(stringr)
library(pbapply)

###fig3a_exon_data_bar_plot
changecounts<-pbapply::pblapply(cl=3,c("M","R","C"),function(spe){
  load(sprintf("/mnt/data6/disk/liuyao/nascentRNA_data/filtersite_speOverlapexon2/%s_v1_DC_sum.Rdata",spe))
  select_D_counts<-sum_DC_counts[,c("BL_Dcounts","BH_Dcounts","LH_Dcounts","BL_BH_Dcounts",
                                    "BL_LH_Dcounts","BH_LH_Dcounts","BLH_Dcounts","BLH_1pair_changecounts","BLH_1pair_allcounts")]
  select_D_counts<-select_D_counts%>%mutate(changeRate_1pair=BLH_1pair_changecounts/BLH_1pair_allcounts)
})%>%rbind.fill()
changeRate_1pair<-changecounts[,8:10]
changecounts<-changecounts[,1:7]
rownames(changecounts)<-c("mouse","rat","cavia")
changecounts$BL_BH_Dcounts<-changecounts$BL_BH_Dcounts-changecounts$BLH_Dcounts
changecounts$BL_LH_Dcounts<-changecounts$BL_LH_Dcounts-changecounts$BLH_Dcounts
changecounts$BH_LH_Dcounts<-changecounts$BH_LH_Dcounts-changecounts$BLH_Dcounts
changecounts$species<-c("mouse","rat","cavia")
melted_data <- reshape2::melt(changecounts, id.vars = "species", variable.name = "group", value.name = "counts")
melted_data$species<-factor(melted_data$species,levels = c("mouse","rat","cavia"))

bar_plot<-ggplot(melted_data,aes(x = group, y = counts, fill = species)) +
  geom_col(width = 0.6,position= position_dodge(width=0.75))+
  scale_x_discrete(limits =c("BL_Dcounts","BH_Dcounts","LH_Dcounts","BL_BH_Dcounts",
                             "BL_LH_Dcounts","BH_LH_Dcounts","BLH_Dcounts"), 
                   labels = c("Brain-liver", "Brain-heart","Heart-liver",
                              "Brain-specific","Liver-specific", "Heart-specific","Three-way differences"))+
  scale_y_continuous(limits =c(0,1e6),breaks = c(0,2e5,4e5,6e5,8e5,1e6),
                     labels = parse(text = c("0","2%*%10^5","4%*%10^5","6%*%10^5","8%*%10^5","1%*%10^6")))+
  labs(x = NULL, y = "Number of sites", fill = "Species") +
  theme(plot.background = element_rect(fill = "white"),panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = NA),
        plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.title = element_text(size=8),
        axis.text = element_text(size=8),legend.text = element_text(size=8),
        axis.text.x = element_text(angle = 15, hjust = 1),
        legend.title = element_text(size=8),legend.position = "none") 
##加右上角表格
library(patchwork)
library(gridExtra)
library(gtable)
library(grid)
library(patchwork)
legend_df <- data.frame(
  Species = c("Mouse", "Rat", "Cavia"),
  Fraction = round(changeRate_1pair$changeRate_1pair*100,2))
legend_df$Color <- c("#F8766D", "#00BA38", "#619CFF")
names(legend_df)[2] <- "Fraction of nucleotides with\ninter-organ FSR differences(%)"
legend_table <- tableGrob(legend_df[, 1:2],rows = NULL,theme = ttheme_default(base_size = 8))
core_id <- which(grepl("^core-fg", legend_table$layout$name))
rows <- legend_table$layout$t[core_id] 
cols <- legend_table$layout$l[core_id] 

for (i in seq_along(legend_df$Species)) {
  idx <- core_id[rows == (i + 1) & cols == 1]
  if (length(idx)) {old_grob <- legend_table$grobs[[idx]]
  legend_table$grobs[[idx]] <-
    gTree(children = gList(rectGrob(gp = gpar(fill = legend_df$Color[i], col = NA)),
                           textGrob(label = old_grob$label,x = old_grob$x, y = old_grob$y,
                                    just = old_grob$just,gp = old_grob$gp)))
  }
}
##把图例 inset 到主图
fig3a<-bar_plot+inset_element(legend_table,left = 0.7, bottom = 0.6,right = 0.8, top =0.6)
ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_new/fig3a_exondata_inset.pdf",
       plot = fig3a,width = 12, height = 6,units = "cm")


###fig3b
load("/mnt/data6/disk/liuyao/nascentRNA_data/filtersite_speOverlapexon2/mrc_betweenspecies_summary.Rdata")
####用organ-specific的counts总和为1按比例画bar图
df<-mrc_withinspecies_summary[1:3,]
df_proportions_matrix <- apply(df[, -1], MARGIN = 2, function(x) x / sum(x, na.rm = TRUE)*100)
df_proportions <- cbind(df[, 1], as.data.frame(df_proportions_matrix))
colnames(df_proportions)[1] <- colnames(df)[1]
Species<- c(rep(c("M_R", "M_C", "R_C","M_R_C"), each = 3))
Organs  <- c(rep(c("Brain","Heart", "Liver"), times = 4))
percent <- c(df_proportions$mr_counts[1:3],df_proportions$mc_counts[1:3],
             df_proportions$rc_counts[1:3],df_proportions$mrc_counts[1:3])
Data      <- data.frame(Species, Organs, percent)
Data$Organs <- factor(Data$Organs,levels = c("Brain", "Liver","Heart"))
Data$Species <- factor(Data$Species,levels = c("M_R","M_C","R_C","M_R_C"))
Data$percent<-round(Data$percent,1)
colors <- c( "#FFD7C4", "#71BBB2","#FFA09B")
fig3b<-ggplot(Data,aes(x = Species, y = percent, fill = Organs)) +
  geom_bar(stat = "identity", position = "stack",width = 0.8) +
  geom_text(aes(label = percent, group = Organs),size = 1.5,
            position = position_stack(vjust = 0.5))+
  labs(x="Overlap between species",
       y="Percentage of overlapping organ-specific\n structures between species(%)")+
  scale_fill_manual(values = colors) +theme_bw() +
  theme(panel.grid.major = element_blank(),legend.text = element_text(size = 7),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.text = element_text(size = 7),axis.title = element_text(size = 7))
ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_0.5/fig3b_2.pdf", 
       plot = fig3b,width = 6.5, height = 5.6,units = "cm")



##fig3c
allsample_average_FSR<-read.table("/mnt/data6/disk/liuyao/nascentRNA_data/FSR/allsample_average_FSR.txt",
                                  header = T,stringsAsFactors = F)
allsample_average_FSR<-allsample_average_FSR[,1:2]
allsample_average_FSR2<-allsample_average_FSR%>%dplyr::mutate(Row=str_split(sample,pattern = "_",simplify = T)[,1],
                                                              Col=str_split(sample,pattern = "_",simplify = T)[,2])
allsample_average_FSR2<-allsample_average_FSR2[,c(3,4,2)]
colnames(allsample_average_FSR2) <- c("Row", "Col", "Average_FSR")
allsample_average_FSR2$Row <- gsub("M", "Mouse", allsample_average_FSR2$Row) 
allsample_average_FSR2$Row <- gsub("R", "Rat", allsample_average_FSR2$Row) 
allsample_average_FSR2$Row <- gsub("C", "Cavia", allsample_average_FSR2$Row)
allsample_average_FSR2$Col <- factor(allsample_average_FSR2$Col, levels = c("brain", "liver", "heart"))
allsample_average_FSR2$Row <- factor(allsample_average_FSR2$Row, levels = c("Rat", "Cavia","Mouse" ))

fig3c<-ggplot(allsample_average_FSR2, aes(x = Col, y = Row, fill = Average_FSR)) +
  geom_tile() +scale_fill_gradientn(colors = c("#DEECFF","#988FF7", "blue"),
                                    breaks = c(seq(0.27,0.35,by=0.02)),
                                    labels = c(seq(0.27,0.35,by=0.02)),
                                    limits = c(0.27, 0.35),values = c(0,0.75,1))+
  geom_text(aes(label = Average_FSR), color = "black", size = 2.5) +theme_minimal() +
  theme(axis.text.x = element_text(),plot.title = element_text(hjust = 0.5,size=8),
        legend.direction = "vertical",legend.position = "right",
        legend.margin = margin(t = 0, r = 0, b = 0, l = -0.5, unit = "cm"),
        legend.key.width = unit(0.3, "cm"),legend.justification = c(0,0.5), 
        legend.key.height = unit(1, "cm"),axis.title = element_text(size=8),
        axis.text = element_text(size=8),legend.text = element_text(size=8),
        legend.title = element_text(size=8, hjust = 1, vjust = 3)) + 
  coord_fixed() +  labs(x = "Organs", y = "Species")+ggtitle("")
ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_0.5/fig3c_average_fsr.pdf", 
       plot = fig3c,width =8 , height =8, units = "cm")

##fig3d
allsample_average_FSR_gini<-read.table("/mnt/data6/disk/liuyao/nascentRNA_data/FSR/allsample_average_FSR_gini.txt",
                                       header = T,stringsAsFactors = F)
allsample_average_FSR1_gini<-allsample_average_FSR_gini[,1:2]
allsample_average_FSR2_gini<-allsample_average_FSR1_gini%>%dplyr::mutate(Row=str_split(sample,pattern = "_",simplify = T)[,1],
                                                                         Col=str_split(sample,pattern = "_",simplify = T)[,2])
allsample_average_FSR2_gini<-allsample_average_FSR2_gini[,c(3,4,2)]
colnames(allsample_average_FSR2_gini) <- c("Row", "Col", "Average_FSR_gini")
allsample_average_FSR2_gini$Row <- gsub("M", "Mouse", allsample_average_FSR2_gini$Row) 
allsample_average_FSR2_gini$Row <- gsub("R", "Rat", allsample_average_FSR2_gini$Row) 
allsample_average_FSR2_gini$Row <- gsub("C", "Cavia", allsample_average_FSR2_gini$Row)
allsample_average_FSR2_gini$Col <- factor(allsample_average_FSR2_gini$Col, levels = c("brain", "liver", "heart"))
allsample_average_FSR2_gini$Row <- factor(allsample_average_FSR2_gini$Row, levels = c("Rat", "Cavia","Mouse" ))

fig3d<-ggplot(allsample_average_FSR2_gini, aes(x = Col, y = Row, fill = Average_FSR_gini)) +
  geom_tile() +scale_fill_gradientn(colors = c("#DEECFF","#988FF7", "blue"),
                                    breaks = c(seq(0.58,0.68,by=0.02)),   
                                    labels = c(seq(0.58,0.68,by=0.02)), 
                                    limits = c(0.58,0.68),values = c(0,0.5,1))+
  geom_text(aes(label = Average_FSR_gini), color = "black", size = 2.5) +theme_minimal() +
  theme(axis.text.x = element_text(),plot.title = element_text(hjust = 0.5,size=8),
        legend.direction = "vertical",legend.position = "right",
        legend.margin = margin(t = 0, r = 0, b = 0, l = -0.5, unit = "cm"),
        legend.key.width = unit(0.3, "cm"),legend.justification = c(0,0.5),  
        legend.key.height = unit(0.9, "cm"),axis.title = element_text(size=8),
        axis.text = element_text(size=8),legend.text = element_text(size=8),
        legend.title = element_text(size=8, hjust = 1, vjust = 3)) + 
  coord_fixed() +  labs(x = "Organs", y = "Species")+ggtitle("")
ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_new/fig3b_average_fsr_gini.pdf", 
       plot = fig3d,width =9 , height =9 , units = "cm")



###fig3e
MR_conservedSite<-read.table("/mnt/data6/disk/liuyao/nascentRNA_data/filtersite_speOverlapexon2/MR_conservedSite.txt")
MR_tissues_oddratio<-pbapply::pblapply(cl=3,c("brain","liver","heart"),function(tissue){
  name=paste(tissue,"specific",sep = "_")
  a_data<-MR_conservedSite%>%dplyr::filter(V7==name)
  b_data<-MR_conservedSite%>%dplyr::filter(V5==name)%>%dplyr::filter(V6!=name)
  c_data<-MR_conservedSite%>%dplyr::filter(V5!=name)%>%dplyr::filter(V6==name)
  d_data<-MR_conservedSite%>%dplyr::filter(V5!=name)%>%dplyr::filter(V6!=name)
  a=nrow(a_data);b=nrow(b_data);c=nrow(c_data);d=nrow(d_data)
  data<-data.frame(Pairwises="mouse_rat",Organs=tissue,a=a,b=b,c=c,d=d,
                   odd_ratio=fisher.test(matrix(c(a,b,c,d),ncol=2,byrow = "T"))$estimate,
                   p_value=fisher.test(matrix(c(a,b,c,d),ncol=2,byrow = "T"))$p.value,
                   conf95CI_1=fisher.test(matrix(c(a,b,c,d),ncol=2,byrow = "T"))$conf.int[1],
                   conf95CI_2=fisher.test(matrix(c(a,b,c,d),ncol=2,byrow = "T"))$conf.int[2],
                   Relative_Risk=(a/(a+b))/(c/(c+d)))
})%>%rbind.fill()

MC_conservedSite<-read.table("/mnt/data6/disk/liuyao/nascentRNA_data/filtersite_speOverlapexon2/MC_conservedSite.txt")
MC_tissues_oddratio<-pbapply::pblapply(cl=3,c("brain","liver","heart"),function(tissue){
  name=paste(tissue,"specific",sep = "_")
  a_data<-MC_conservedSite%>%dplyr::filter(V7==name)
  b_data<-MC_conservedSite%>%dplyr::filter(V5==name)%>%dplyr::filter(V6!=name)
  c_data<-MC_conservedSite%>%dplyr::filter(V5!=name)%>%dplyr::filter(V6==name)
  d_data<-MC_conservedSite%>%dplyr::filter(V5!=name)%>%dplyr::filter(V6!=name)
  a=nrow(a_data);b=nrow(b_data);c=nrow(c_data);d=nrow(d_data)
  data<-data.frame(Pairwises="mouse_cavia",Organs=tissue,a=a,b=b,c=c,d=d,
                   odd_ratio=fisher.test(matrix(c(a,b,c,d),ncol=2,byrow = "T"))$estimate,
                   p_value=fisher.test(matrix(c(a,b,c,d),ncol=2,byrow = "T"))$p.value,
                   conf95CI_1=fisher.test(matrix(c(a,b,c,d),ncol=2,byrow = "T"))$conf.int[1],
                   conf95CI_2=fisher.test(matrix(c(a,b,c,d),ncol=2,byrow = "T"))$conf.int[2],
                   Relative_Risk=(a/(a+b))/(c/(c+d)))
})%>%rbind.fill()


RC_conservedSite<-read.table("/mnt/data6/disk/liuyao/nascentRNA_data/filtersite_speOverlapexon2/RC_conservedSite.txt")
RC_tissues_oddratio<-pbapply::pblapply(cl=3,c("brain","liver","heart"),function(tissue){
  name=paste(tissue,"specific",sep = "_")
  a_data<-RC_conservedSite%>%dplyr::filter(V7==name)
  b_data<-RC_conservedSite%>%dplyr::filter(V5==name)%>%dplyr::filter(V6!=name)
  c_data<-RC_conservedSite%>%dplyr::filter(V5!=name)%>%dplyr::filter(V6==name)
  d_data<-RC_conservedSite%>%dplyr::filter(V5!=name)%>%dplyr::filter(V6!=name)
  a=nrow(a_data);b=nrow(b_data);c=nrow(c_data);d=nrow(d_data)
  data<-data.frame(Pairwises="rat_cavia",Organs=tissue,a=a,b=b,c=c,d=d,
                   odd_ratio=fisher.test(matrix(c(a,b,c,d),ncol=2,byrow = "T"))$estimate,
                   p_value=fisher.test(matrix(c(a,b,c,d),ncol=2,byrow = "T"))$p.value,
                   conf95CI_1=fisher.test(matrix(c(a,b,c,d),ncol=2,byrow = "T"))$conf.int[1],
                   conf95CI_2=fisher.test(matrix(c(a,b,c,d),ncol=2,byrow = "T"))$conf.int[2],
                   Relative_Risk=(a/(a+b))/(c/(c+d)))
})%>%rbind.fill()

mr_mc_oddratio<-rbind(MR_tissues_oddratio,MC_tissues_oddratio)
mr_mc_rc_oddratio<-rbind(mr_mc_oddratio,RC_tissues_oddratio)
mr_mc_rc_oddratio$log2OR<-log2(mr_mc_rc_oddratio$odd_ratio)
mr_mc_rc_oddratio$Organs<-factor(mr_mc_rc_oddratio$Organs,levels = c("brain","liver","heart"))
mr_mc_rc_oddratio$Pairwises<-factor(mr_mc_rc_oddratio$Pairwises,levels = c("mouse_rat","mouse_cavia","rat_cavia"))

mr_mc_rc_oddratio_adj <- mr_mc_rc_oddratio %>%
  mutate(odds_ratio_adj = odd_ratio - 1)

fig3e<-ggplot(mr_mc_rc_oddratio_adj,aes(x = Organs, y = odds_ratio_adj, fill = Pairwises)) +
  geom_col(width = 0.8, position= position_dodge(width=0.9))+
  geom_hline(yintercept = 0, linetype = "dashed",color = "grey",linewidth = 0.5)+
  scale_y_continuous(limits = c(-0.2,0.7), breaks = seq(-0.2,0.6,0.2),labels = seq(0.8,1.6,0.2))+
  coord_flip() +labs(x = "", y = "Conservation of organ-specificity\n(Odds Ratio)", fill = "Species pair") +
  scale_fill_manual(values = c("#F7D8B7","#CF9CB4","#7985B3"))+
  theme(plot.background = element_rect(fill = "white"),panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", linewidth = 0.8, fill = NA),
        plot.title = element_text(hjust = 0.5),panel.grid= element_blank(),
        axis.title = element_text(size=8),axis.text = element_text(size=8),
        legend.text = element_text(size=8),legend.title = element_text(size=8),
        axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) 

ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_0.5/fig3e_odfrom1.pdf", 
       plot = fig3e,width = 8.7, height = 5.2,units = "cm")


###fig3f
load("/mnt/data6/disk/liuyao/nascentRNA_tissue/exon_windows_p/windows_10_50/filter_w10_5_one_exon.Rdata")
window_Fsite=as.numeric(str_split(filter_w10_5$windows[1],pattern = "_",simplify = T)[,1])-21
window_Rsite=as.numeric(str_split(filter_w10_5$windows[1],pattern = "_",simplify = T)[,2])+5
for( i in 4:12){filter_w10_5_one_exon_fsr[,i]<-filter_w10_5_one_exon_fsr[,i]/(max(na.omit(filter_w10_5_one_exon_fsr[,i])))}
one_exon_fsr_w<-filter_w10_5_one_exon_fsr%>%dplyr::filter(mouse_orthologsite%in%c(window_Fsite:window_Rsite))
one_exon_fsr_w1<-one_exon_fsr_w[,c(3:12)]
colnames(one_exon_fsr_w1)<-c("mouse_orthologsite","M_B","M_L","M_H","R_B","R_L","R_H","C_B","C_L","C_H")
long_dt <- reshape2::melt(one_exon_fsr_w1,id.vars = "mouse_orthologsite",variable.name = "sample", value.name = "FSR")

if (filter_w10_5$mouse=="liver"){
  long_dt$sample<-factor(long_dt$sample,levels=c("M_B","M_H","M_L","R_B","R_H","R_L","C_B","C_H","C_L"))
}else if (filter_w10_5$mouse=="brain"){
  long_dt$sample<-factor(long_dt$sample,levels=c("M_L","M_H","M_B","R_L","R_H","R_B","C_L","C_H","C_B"))
}

heart_color="#FF6666"
highlight_region <- data.frame(xmin = 13.5,xmax = 20.5, ymin = -Inf,ymax = Inf)
fig3f<-ggplot(long_dt,aes(x=mouse_orthologsite,y=FSR,fill=sample))+
  geom_bar(stat = "identity", width = 1)+facet_grid(sample~.)+
  scale_y_continuous(breaks = c(0,1),labels = c(0,1),limits = c(0,1), expand = c(0, 0.05)) + 
  scale_fill_manual(values = c(rep("#C599B6",2),heart_color,rep("#A5B68D",2),
                               heart_color,rep("#66B3FF",2),heart_color)) +
  labs(x="Nucleotide position (1–25) within exon: 19,291–19,869 of Kmt2b")+
  geom_vline(xintercept = 13.5, linetype = "dashed", color = "grey", linewidth = 0.5)+
  geom_vline(xintercept = 20.5, linetype = "dashed", color = "grey", linewidth = 0.5)+
  theme(plot.title = element_text(hjust = 0.5),panel.border = element_blank(),
        panel.grid = element_blank(),axis.title = element_text(size=8),
        axis.text = element_text(size=8),panel.spacing = unit(0.3, "lines"),
        strip.text = element_text(size = 8),legend.position = "none",
        strip.text.y = element_text(angle = 0,hjust = 0.5,vjust = 0.5),
        panel.background = element_rect(fill = "gray100", color = "gray80", linewidth = 0.5),
        strip.background = element_rect(fill = "gray97", color = "gray70"))+
  geom_rect(data = highlight_region,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", alpha = 0.35,inherit.aes = FALSE)
ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_0.5/fig3f_example.pdf", 
       plot = fig3f,width = 15, height = 10,units = "cm")




