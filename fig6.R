library(reshape2)
####fig6a
Hsp_d3_liverFSR<-read.table("/mnt/data6/disk/liuyao/nascentRNA_tissue/exon_intron_FSR/sample/ENSMUSG00000020048_d3ss_50nt_liverFSR.txt")
Hsp_d3_brainFSR<-read.table("/mnt/data6/disk/liuyao/nascentRNA_tissue/exon_intron_FSR/sample/ENSMUSG00000020048_d3ss_50nt_brainFSR.txt")
Hsp_d3_heartFSR<-read.table("/mnt/data6/disk/liuyao/nascentRNA_tissue/exon_intron_FSR/sample/ENSMUSG00000020048_d3ss_50nt_heartFSR.txt")

colnames(Hsp_d3_liverFSR)<-c("site","Liver")
colnames(Hsp_d3_brainFSR)<-c("site","Brain")
colnames(Hsp_d3_heartFSR)<-c("site","Heart")
Hsp_d3_lbFSR<-merge(Hsp_d3_liverFSR,Hsp_d3_brainFSR,by="site")
Hsp_d3_lbhFSR<-merge(Hsp_d3_lbFSR,Hsp_d3_heartFSR,by="site")
Hsp_long_dt <- reshape2::melt(Hsp_d3_lbhFSR,id.vars = "site",variable.name = "Organs", value.name = "FSR")
Hsp_long_dt$Organs<-factor(Hsp_long_dt$Organs,levels = c("Heart","Brain","Liver"))

fig6a<-ggplot(Hsp_long_dt, aes(x = site, y = FSR, fill = Organs)) +
  geom_bar(stat ="identity",width=1,aes(color=Organs),linewidth=0.4)+facet_grid(Organs ~ .) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5),labels = c(0, 0.5, 1, 1.5),
                     limits = c(0, 1.5),expand = c(0, 0.05)) +labs(x="")+
  theme(plot.title = element_text(hjust = 0.5),panel.border = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.title = element_text(size = 8),axis.text = element_text(size = 8),
        panel.spacing = unit(0.3, "lines"),legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),strip.text = element_text(size = 8),
        panel.background = element_rect(fill = "gray100", color = "gray80", linewidth = 0.5),
        legend.key.height = unit(0.4, "cm"),legend.key.width = unit(0.4, "cm"),
        strip.background = element_rect(fill = "gray90", color = "gray70"))

ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_new/fig6a.pdf", 
       plot = fig6a,width = 11, height = 5,units = "cm")

#####fig6b
mytheme=
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = NA),
        plot.title = element_text(hjust = 0.5, size = 8),
        panel.grid = element_blank(),axis.title = element_text(size=8),
        axis.text = element_text(size=8),legend.text = element_text(size=8),
        legend.title = element_text(size=8))

load("/mnt/data6/disk/liuyao/nascentRNA_tissue/exon_intron_FSR/oneSS_FSR/oneSS_fsr_se_rho.Rdata")
oneSS_fsr_se<-oneSS_fsr_se_rho[,c(8,2:7)]
##计算FSR的z-score
oneSS_fsr_se_z <- oneSS_fsr_se %>%
  rowwise() %>%dplyr::mutate(
    heart_d3SS_z = (heart_d3SS - mean(c(heart_d3SS, brain_d3SS, liver_d3SS))) / 
      sd(c(heart_d3SS, brain_d3SS, liver_d3SS)),
    brain_d3SS_z = (brain_d3SS - mean(c(heart_d3SS, brain_d3SS, liver_d3SS))) / 
      sd(c(heart_d3SS, brain_d3SS, liver_d3SS)),
    liver_d3SS_z = (liver_d3SS - mean(c(heart_d3SS, brain_d3SS, liver_d3SS))) / 
      sd(c(heart_d3SS, brain_d3SS, liver_d3SS)))%>%ungroup()

se_long <- oneSS_fsr_se_z %>%select(species, heart_SE, brain_SE, liver_SE) %>%
  pivot_longer(cols = c(heart_SE, brain_SE, liver_SE),names_to = "organ",values_to = "SE") %>%
  dplyr::mutate(organ = gsub("_SE$", "", organ)) 

d3ss_long <- oneSS_fsr_se_z %>%select(species, heart_d3SS_z, brain_d3SS_z, liver_d3SS_z) %>%
  pivot_longer(cols = c(heart_d3SS_z, brain_d3SS_z, liver_d3SS_z),
               names_to = "organ",values_to = "d3SS_z") %>%
  dplyr::mutate(organ = gsub("_d3SS_z$", "", organ))

oneSS_long <- se_long %>%inner_join(d3ss_long, by = c("species", "organ")) %>%
  mutate(organ = factor(organ, levels = c("heart", "brain", "liver")),
         species= factor(species, levels = c("mouse", "rat", "cavia")))

fig6b<-ggplot(oneSS_long,aes(x = d3SS_z, y = SE, group = species, color = species)) +
  geom_point(aes(shape = organ),size = 2 )+geom_line()+ mytheme +
  labs(x = "Z-score of Mean FSR\n(50nt downstream of 3'SS)", y = "Splicing Efficiency",
       color = "Species", shape = "Organ") 

ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_new/fig6/fig6b_zFSR_SE_hsp90b1_mrc.pdf",
       plot = fig6b,width = 9,height = 6.5,units = "cm")

####fig6e
####ASO_ENSMUSG00000020048_2ASOs
remove_farthest_from_median <- function(x) {
  med <- median(x);distances <- abs(x - med)
  farthest_index <- which.max(distances)
  return(x[-farthest_index])}

ASO_048_1<-read.table("/mnt/data6/disk/liuyao/nascentRNA_tissue/exon_intron_FSR/sample/ASO_048_1_result.txt",header = T,sep = "\t")
ASO_048_2<-read.table("/mnt/data6/disk/liuyao/nascentRNA_tissue/exon_intron_FSR/sample/ASO_048_2_result.txt",header = T,sep = "\t")
ASO_048_1<-ASO_048_1[,c(3,4,6,8)];ASO_048_2<-ASO_048_2[,c(3,4,6,8)]
ASO_048<-rbind(ASO_048_1,ASO_048_2)
ASO_048<-ASO_048%>%mutate(group=c(rep(1:(nrow(ASO_048)/3),each=3)))

ASO_048_cqMean=as.data.frame(matrix(ncol=3))
colnames(ASO_048_cqMean)<-c("sample","gene","cqMean")
for (i in 1:(nrow(ASO_048)/3)){
  ASO_048_1<-ASO_048%>%dplyr::filter(group==i)
  if (ASO_048_1$Cq.Error[1]>0.1){
    remove_one<-remove_farthest_from_median(ASO_048_1$Cq)
    if (diff(remove_one)<=0.1){cqMean=round(mean(remove_one),2)}else{cqMean=NA}
  }else{cqMean=round(mean(ASO_048_1$Cq),2)}
  ASO_048_cqMean[i,]<-list(ASO_048_1$Sample.Name[1],ASO_048_1$Gene.Name[1],cqMean)
}

ASO_048_sample<-ASO_048_cqMean%>%dplyr::filter(!gene=="actin")
ASO_048_ref<-ASO_048_cqMean%>%dplyr::filter(gene=="actin")
colnames(ASO_048_ref)[3]<-"ref_cq"

ASO_048_sample2<-merge(ASO_048_sample,ASO_048_ref[,c(1,3)],by = "sample")
ASO_048_sample2<-ASO_048_sample2%>%mutate(delta_cq=cqMean-ref_cq,FC=2^(ref_cq-cqMean))
ASO_data<-ASO_048_sample2%>%mutate(group=str_split(sample,pattern = "_",simplify = T)[,1])
ASO_data<-ASO_048_sample2%>%mutate(group=substr(sample, 1, nchar(sample) - 5))
ASO_data1<-ASO_data[,c(1,2,5)]
ASO_data1$gene<-rep(c("spliced", "unspliced"),times=nrow(ASO_data1)/2)
ASO_data1_reps_SE <- ASO_data1 %>%
  pivot_wider(names_from = gene, values_from = delta_cq) %>%
  mutate(delta_Cq = spliced - unspliced,splicing_ratio = 2^(-delta_Cq),
         spliced = (splicing_ratio / (1 + splicing_ratio)) * 100,unspliced = 100 - spliced,
         ratio= unspliced/spliced,group=substr(sample, 1, nchar(sample) - 5)) %>%
  mutate(Sample = case_when(group == "NC" ~ "Control",group == "048_1" ~ "ASO1",group == "048_2" ~ "ASO2"),
         Sample = factor(Sample, levels = c("Control", "ASO1", "ASO2")))

t.test(ASO_data1_reps_SE$ratio[1:3],ASO_data1_reps_SE$ratio[7:9])##p-value = 0.00404
t.test(ASO_data1_reps_SE$ratio[4:6],ASO_data1_reps_SE$ratio[7:9])##p-value = 0.001478

ASO_reps_SE_ratio<-ASO_data1_reps_SE%>%group_by(Sample)%>%
  summarise(mean_ratio=mean(ratio),sd_ratio=sd(ratio))
signif_data <- data.frame(xmin = c(1, 1),xmax = c(2, 3),y = c(0.52, 0.6)) 
fig6e<-ggplot(ASO_reps_SE_ratio, aes(x = Sample, y = mean_ratio)) +
  geom_col(aes(fill = Sample), position = position_dodge(0.8), width = 0.4, show.legend = FALSE) +
  scale_fill_manual(values = c("grey80", "#A7C7CE", "#A7C7CE")) +
  geom_errorbar(aes(ymin = mean_ratio - sd_ratio, ymax = mean_ratio + sd_ratio),
                position = position_dodge(0.8), width = 0.15, linewidth = 0.5) +
  geom_point(data = ASO_data1_reps_SE,aes(x = Sample, y = ratio),  
             position = position_jitter(width = 0.1),size = 0.8, alpha = 0.6) +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.1))+
  theme_bw() + xlab("")+ylab("Unspliced/Spliced")+ggtitle(expression(italic("Hsp90b1")))+
  theme(panel.grid= element_blank(),plot.title = element_text(hjust = 0.5, size = 8),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = NA),
        axis.title = element_text(size = 8),axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),legend.title = element_text(size = 8)) +
  geom_segment(data = signif_data, aes(x = xmin, xend = xmax, y = y, yend = y),inherit.aes = FALSE)+
  geom_segment(data = signif_data, aes(x = xmin, xend = xmin, y = y - 0.02, yend = y),inherit.aes = FALSE, linewidth = 0.5) + #添加左端竖线
  geom_segment(data = signif_data, aes(x = xmax, xend = xmax, y = y - 0.02, yend = y),  inherit.aes = FALSE, linewidth = 0.5) +# 添加右端竖线
  annotate("text", x = 1.5, y = 0.54, label = "italic(P)~'= 0.004'", parse = TRUE, size = 2) +
  annotate("text", x = 2, y = 0.62, label = "italic(P)~'= 0.001'", parse = TRUE, size = 2)

ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_0.5/fig6e.pdf", 
       plot = fig6e,width =6.4 , height =6.5 , device = cairo_pdf,units = "cm")  

