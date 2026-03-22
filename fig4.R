###fig4a
#####mouse
spe="M"
lbh=as.data.frame(matrix(nrow = 0,ncol = 4))
colnames(lbh)<-c("name","group","mean_fsr","tissue")
for (tissue in c("heart","liver","brain")){
  MeanFSR_add_Mean_Random<-read.table(sprintf("/mnt/data6/disk/liuyao/nascentRNA_tissue/splicingsite_FSR/exon_fsr_se/03.add_Random/%s_%s_exon_add_1Random_fsr.txt",spe,tissue),
                                      header = T,stringsAsFactors = F)
  colnames(MeanFSR_add_Mean_Random)[5:12]<-c("U_5SS","D_5SS","U_3SS","D_3SS","U_5SS_r","D_5SS_r","U_3SS_r","D_3SS_r")
  for (col in c(5:8)){
    r_col=col+4
    MeanFSR_add_Mean_Random_f<-na.omit(MeanFSR_add_Mean_Random[,c(1:4,col,r_col)])
    tissue_data<-pblapply(cl=2,c(5:6),function(n){
      aa<-data.frame(name=colnames(MeanFSR_add_Mean_Random_f)[n],
                     group=unlist(str_split(colnames(MeanFSR_add_Mean_Random_f)[n],pattern = "_r"))[1],
                     mean_fsr=na.omit(MeanFSR_add_Mean_Random_f[,n]))
    })%>%rbind.fill()
    tissue_data<-tissue_data%>%mutate(tissue=gsub("^(.)", "\\U\\1", tissue, perl = TRUE))
    lbh<-rbind(lbh,tissue_data)
  }
}

lbh$tissue <- factor(lbh$tissue, levels = c("Brain", "Liver", "Heart"))
lbh$group <- factor(lbh$group, levels = c("U_5SS", "D_5SS", "U_3SS", "D_3SS"))
lbh <- lbh %>%mutate(facet_var = paste(tissue, name, sep = "_"))
median_data1 <- lbh %>%dplyr::group_by(name, tissue, group) %>% 
  dplyr::summarise(median_val = median(mean_fsr))
median_data1<-median_data1%>%mutate(names=paste(tissue,name,sep = "_" ))

fill_colors <- c("Brain_D_5SS"= "#5F8B4C", "Brain_U_5SS_r" = "#313647", "Brain_U_5SS"= "#8CE4FF" ,"Brain_D_5SS_r"  = "gray",
                 "Brain_D_3SS"= "#D78FEE", "Brain_U_3SS_r" = "gray", "Brain_U_3SS"= "#F57C61" ,"Brain_D_3SS_r"  = "#313647",
                 "Liver_D_5SS"= "#5F8B4C", "Liver_U_5SS_r" = "#313647", "Liver_U_5SS"= "#8CE4FF" ,"Liver_D_5SS_r"  = "gray",
                 "Liver_D_3SS"= "#D78FEE", "Liver_U_3SS_r" = "gray", "Liver_U_3SS"= "#F57C61" ,"Liver_D_3SS_r"  = "#313647",
                 "Heart_D_5SS"= "#5F8B4C", "Heart_U_5SS_r" = "#313647", "Heart_U_5SS"= "#8CE4FF" ,"Heart_D_5SS_r"  = "gray",
                 "Heart_D_3SS"= "#D78FEE", "Heart_U_3SS_r" = "gray", "Heart_U_3SS"= "#F57C61" ,"Heart_D_3SS_r"  = "#313647")
combined_title <- paste0("Mouse\n","5' splice (donor) site",
                         strrep(" ", 40), "3' splice (acceptor) site") 
fig4a<-ggplot(lbh, aes(x = mean_fsr, color = facet_var)) +
  geom_density(linewidth = 0.8, alpha = 0.8, adjust = 2) +
  facet_grid(rows = vars(tissue),cols = vars(group),
             labeller = labeller(group = c("U_5SS" = "50nts upstream","D_5SS" = "50nts downstream",
                                           "U_3SS" = "50nts upstream","D_3SS" = "50nts downstream"))) +
  scale_color_manual(values = fill_colors) +
  geom_vline(data = median_data1,  aes(xintercept = median_val, color = names),
             linetype = "dashed", linewidth = 0.8, alpha = 0.7) +
  labs(x = "Average FSR",y = "Density",color = "Group",title = combined_title) +theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 20),size = 12 ),
        strip.text.y.right = element_text(size = 8,angle = -90,hjust = 0.5),
        strip.text = element_text(size = 8),panel.grid = element_blank(),
        strip.background = element_rect(fill = "gray95", color = "black"),
        axis.text = element_text(size = 8),axis.title = element_text(size = 10),
        legend.text = element_text(size = 8),legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))

ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_0.5/fig4a_density_1random.pdf", 
       plot = fig4a,width = 26, height = 12,units = "cm")







####fig4b 
load("/mnt/data6/disk/liuyao/nascentRNA_tissue/splicingsite_FSR/exon_fsr_se/02.exondata_fsr_se/exondata_ss_fsr_se_s.Rdata")
mrc_spearman_cor_melt<-as.data.frame(matrix(nrow = 0, ncol = 5))
colnames(mrc_spearman_cor_melt)<-c("tissue","rho","p_value","location","significance")
for (l in c(1,4,7)){
  mouse_spearman_cor<-exondata_ss_fsr_se_rho[l:(l+2),]
  mouse_spearman_cor_melt<-pbapply::pblapply(cl=4,c(2,4,6,8),function(i){
    one<-mouse_spearman_cor[,c(1,i,(i+1))]
    name=str_split(colnames(one)[2],pattern = "_rho",simplify = T)[,1]
    colnames(one)<-c("sample","rho","v");one$v<-as.numeric(one$v)
    one<-one%>%dplyr::mutate(location=name)
  })%>%rbind.fill()
  mrc_spearman_cor_melt<-rbind(mrc_spearman_cor_melt,mouse_spearman_cor_melt)
}
mrc_spearman_cor_melt$p_adjusted_bonf <- p.adjust(mrc_spearman_cor_melt$v, method = "bonferroni")
mrc_spearman_cor_melt<-mrc_spearman_cor_melt%>%dplyr::mutate(significance_adjusted=ifelse(p_adjusted_bonf<0.05,"P<0.05","P>=0.05"))
mrc_spearman_cor_melt$rho<-as.numeric(mrc_spearman_cor_melt$rho)
mrc_spearman_cor_melt1<-mrc_spearman_cor_melt%>%dplyr::mutate(Species=str_split(sample,pattern = "_",simplify = T)[,1],
                                                       Organs=str_split(sample,pattern = "_",simplify = T)[,2])

mrc_spearman_cor_melt1$Species[mrc_spearman_cor_melt1$Species == "M"] <- "Mouse"
mrc_spearman_cor_melt1$Species[mrc_spearman_cor_melt1$Species == "R"] <- "Rat"
mrc_spearman_cor_melt1$Species[mrc_spearman_cor_melt1$Species == "C"] <- "Cavia"

group_means <- mrc_spearman_cor_melt1 %>%dplyr::group_by(location) %>%
  dplyr::summarise(mean_rho = mean(rho, na.rm = TRUE),.groups = 'drop')
y_groups <- unique(mrc_spearman_cor_melt1$location)

mrc_spearman_cor_melt1$Species<-factor(mrc_spearman_cor_melt1$Species,levels = c("Mouse","Rat","Cavia"))
mrc_spearman_cor_melt1$Organs<-factor(mrc_spearman_cor_melt1$Organs,levels = c("brain","liver","heart"))
mrc_spearman_cor_melt1$location<-factor(mrc_spearman_cor_melt1$location,levels = c("down_3ss","up_3ss","down_5ss","up_5ss"))

fig4b<-ggplot(mrc_spearman_cor_melt1, aes(x = rho, y =location )) +
  geom_hline(yintercept = y_groups, linetype = "dashed", color = "gray50", linewidth = 0.3)+
  geom_vline(xintercept = 0,linetype = "dashed",color = "black",linewidth = 0.4)+
  geom_point(data = filter(mrc_spearman_cor_melt1, significance_adjusted == "P<0.05"),
             aes(shape = Species, color = Organs, fill = Organs),size = 3, stroke = 0.8) +
  geom_point(data = filter(mrc_spearman_cor_melt1, significance_adjusted == "P>=0.05"), 
             aes(shape = Species, color = Organs),size = 3, fill = NA, stroke = 0.8) +
  scale_shape_manual(name = "Species",values = c(Mouse = 21, Rat = 22, Cavia = 24)) + 
  geom_point(data = group_means,aes(x = mean_rho, y = location),shape = 23, size = 5.3,color = "#3C3D37",fill = "#3C3D37")+
  scale_color_manual(name = "Organs",values = c(heart = "#C8B8D4", brain = "#CCF3FF", liver = "#FFCCCC"),
                     aesthetics = c("color", "fill")) +
  scale_y_discrete(labels = c("up_5ss" = "50-nt region \nupstream of 5'SS",
                              "down_5ss" = "50-nt region \ndownstream of 5'SS", "up_3ss" = "50-nt region \nupstream of 3'SS",
                              "down_3ss" = "50-nt region \ndownstream of 3'SS"))+theme_bw() +
  theme(panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),axis.text.y = element_text(size = 8,hjust = 0.5,margin = margin(r = 0)),
        axis.title = element_text(size=8),axis.text = element_text(size=8),
        legend.text = element_text(size=8),legend.title = element_text(size=8)) +
  labs(x = "The spearman correlation between the average FSR across the \n50-nt flanking regions of 5'/3' splice sites and splicing efficiency",
       y = "") +guides(shape = guide_legend(title = "Species", override.aes = list(color = "black", fill = "white")),
                       color = guide_legend(title = "Organs",override.aes = list(shape = 21, size = 4))) 
ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_0.5/fig4b.pdf", 
       plot = fig4b,width = 13, height = 8,units = "cm")



####fig4c
breaks <- c(-Inf,0.15,0.25,0.35,0.45,Inf)
breaks_region<-c("(-Inf,0.15]","(0.15,0.25]","(0.25,0.35]","(0.35,0.45]","(0.45,Inf]")
spe="M";tissue="brain"
load(sprintf("/mnt/data6/disk/liuyao/nascentRNA_tissue/splicingsite_FSR/exon_fsr_se/02.exondata_fsr_se/%s_%s_ss_nearfsr_se.Rdata",spe,tissue))
se_fsr_u5<-na.omit(fsr_se[,c(7,11)]);se_fsr_u5<-se_fsr_u5%>%mutate(location=paste(spe,tissue,"U_5SS",sep = "_"));colnames(se_fsr_u5)[1]<-"mean_fsr"
se_fsr_u5$group <- cut(se_fsr_u5$mean_fsr, breaks = breaks, include.lowest = TRUE,labels = breaks_region)

spe="C";tissue="brain"
load(sprintf("/mnt/data6/disk/liuyao/nascentRNA_tissue/splicingsite_FSR/exon_fsr_se/02.exondata_fsr_se/%s_%s_ss_nearfsr_se.Rdata",spe,tissue))
se_fsr_d5<-na.omit(fsr_se[,c(9,11)]);se_fsr_d5<-se_fsr_d5%>%mutate(location=paste(spe,tissue,"D_5SS",sep = "_"));colnames(se_fsr_d5)[1]<-"mean_fsr"
se_fsr_d5$group <- cut(se_fsr_d5$mean_fsr, breaks = breaks, include.lowest = TRUE,labels = breaks_region)

alldata_se_fsr<-rbind(se_fsr_u5,se_fsr_d5)
alldata_se_fsr$location<-factor(alldata_se_fsr$location,levels = c("M_brain_U_5SS","C_brain_D_5SS"))

fig_4c<-ggplot(alldata_se_fsr, aes(x = group, y = spilcing_efficiency)) +
  geom_boxplot(width = 0.6, fatten = 1,outlier.size = 0.5) + 
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.25))+
  facet_wrap(~ location, ncol = 1,
             labeller = labeller(location = c("M_brain_U_5SS" = "Upstream of 5'SS in mouse brain",
                                              "C_brain_D_5SS" = "Downstream of 5'SS in cavia brain")))+
  labs(x = "Average FSR across 50-nt \nflanking regions of 5' splice sites", 
       y = "Spilcing efficiency") +theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 8),panel.grid = element_blank(),
        axis.title = element_text(size = 8), axis.text.x = element_text(angle = 25, hjust = 1),
        axis.text = element_text(size=8),legend.text = element_text(size=8),
        legend.title = element_text(8),strip.text = element_text(size=7),
        strip.background = element_rect(fill = NA, color = "black"),legend.position = "none")

ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_0.5/fig_4c.pdf", 
       plot = fig_4c,width = 6, height = 8.5,units = "cm")

