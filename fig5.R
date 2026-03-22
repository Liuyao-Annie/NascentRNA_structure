###fig5a
mrc_od2<-read.table("/mnt/data6/disk/liuyao/nascentRNA_tissue/splicingsite_FSR/exon_fsr_se/04.oddsratio/se_fsr_oddsratio_new.txt",
                    header = T,stringsAsFactors = F)
mrc_od2<-mrc_od2%>%mutate(species=str_split(sample,pattern = "_",simplify = T)[,1],
                          organs=str_split(sample,pattern = "_",simplify = T)[,2])
mrc_od2$species[mrc_od2$species == "M"] <- "Mouse"
mrc_od2$species[mrc_od2$species == "R"] <- "Rat"
mrc_od2$species[mrc_od2$species == "C"] <- "Cavia"
mrc_od2$species<-factor(mrc_od2$species,levels = c("Mouse","Rat","Cavia"))
mrc_od2$organs<-factor(mrc_od2$organs,levels = c("brain","liver","heart"))
mrc_od2$location<-factor(mrc_od2$location,levels = c("up_5SS","down_5SS","up_3SS","down_3SS"))

fig5_theme<-theme(plot.background = element_rect(fill = "white"),legend.position = "none",
                  panel.background = element_rect(fill = "white"),
                  plot.title = element_text(hjust = 0.5),panel.grid = element_blank(),
                  axis.title = element_text(size = 8),axis.text = element_text(size = 8),
                  legend.text = element_text(size = 8),legend.title = element_text(size = 8),
                  strip.background = element_rect(fill = "white",linewidth = 0.4),
                  strip.text = element_text(size = 8))

mrc_od2_adj <- mrc_od2 %>%
  mutate(odds_ratio_adj = odds_ratio - 1,
         bootstrap_OR_25_adj = bootstrap_OR_25 - 1,
         bootstrap_OR_975_adj = bootstrap_OR_975 - 1)

fig_5a<-ggplot(mrc_od2_adj, aes(x = organs, y = odds_ratio_adj, fill = species)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(ymin = bootstrap_OR_25_adj, ymax = bootstrap_OR_975_adj),
                width = 0.2, position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(-1,1.6),breaks = c(-1,-0.5,0,0.5, 1,1.5), 
                     labels = c(0,0.5, 1,1.5, 2,2.5))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) + 
  facet_wrap(~ location, ncol = 2, scales = "free_y",
             labeller = labeller(location = c("up_5SS" = "50nts upstream of 5'SS",
                                              "down_5SS" = "50nts downstream of 5'SS",
                                              "up_3SS" = "50nts upstream of 3'SS",
                                              "down_3SS" = "50nts downstream of 3'SS"))) +
  labs(x = NULL, y = "odds ratio") +theme_bw()+fig5_theme  

ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_0.5/fig_5a_from1_col1.pdf", 
       plot = fig_5a,width = 16, height = 12,units = "cm")


####fig5b
random_percol<-function(mr_blh_se_fsr_ss){
  for( i in c(4:ncol(mr_blh_se_fsr_ss))){
    mr_blh_se_fsr_ss[,i]<-sample(mr_blh_se_fsr_ss[,i])
  }
  return(mr_blh_se_fsr_ss)
}

spearman<-function(data1,data2){
  cor_test_result <- suppressWarnings(
    try(cor.test(data1,data2, method = "s"), silent = TRUE))
  rho <- if(inherits(cor_test_result, "try-error")) NA_real_ else cor_test_result$estimate
  return(rho)
}
pairorgans_FSR_SE_match_introns<-function(mr_blh_se_fsr_ss,x){
  SE_match<-c();FSR_SE_match<-c()
  for (line in 1:nrow(mr_blh_se_fsr_ss)){
    oneline=mr_blh_se_fsr_ss[line,]
    data<-data.frame(FSR1=c(oneline[1,4],oneline[1,6],oneline[1,8]),SE1=c(oneline[1,5],oneline[1,7],oneline[1,9]),
                     FSR2=c(oneline[1,10],oneline[1,12],oneline[1,14]),SE2=c(oneline[1,11],oneline[1,13],oneline[1,15]))
    se_rho<-spearman(data$SE1,data$SE2)
    if (!is.na(se_rho) && se_rho > 0){
      SE_match<-c(SE_match,oneline$mouse_intron)
      fsr1_rho<-spearman(data$FSR1,data$SE1);fsr2_rho<-spearman(data$FSR2,data$SE2)
      if (x==2){
        if(!is.na(fsr1_rho) && !is.na(fsr2_rho) && fsr1_rho > 0 && fsr2_rho > 0){
          FSR_SE_match<-c(FSR_SE_match,oneline$mouse_intron)
        } 
      }else{
        if(!is.na(fsr1_rho) && !is.na(fsr2_rho) && fsr1_rho < 0 && fsr2_rho < 0){
          FSR_SE_match<-c(FSR_SE_match,oneline$mouse_intron)
        } 
      }
    }
  }
  match_introns<-list(SE_match_intronnumber=length(SE_match),
                      FSR_SE_match_intronnumber=length(FSR_SE_match),
                      SE_match_introns=SE_match,FSR_SE_match_introns=FSR_SE_match)
  return(match_introns)
}

spe="M"
m_blh_fsr_se<-read.table(sprintf("/mnt/data6/disk/liuyao/nascentRNA_tissue/splicingsite_FSR/exon_fsr_se/02.exondata_fsr_se/%s_blh_fsr_se_add.txt",spe),
                         header = T,stringsAsFactors = F);m_blh_fsr_se<-m_blh_fsr_se[,c(1,2,15,16,3:5,13,14,6:8,11,12,9:10)]
spe="R"
r_blh_fsr_se<-read.table(sprintf("/mnt/data6/disk/liuyao/nascentRNA_tissue/splicingsite_FSR/exon_fsr_se/02.exondata_fsr_se/%s_blh_fsr_se_add.txt",spe),
                         header = T,stringsAsFactors = F);r_blh_fsr_se<-r_blh_fsr_se[,c(1,2,15,16,3:5,13,14,6:8,11,12,9:10)]
spe="C"
c_blh_fsr_se<-read.table(sprintf("/mnt/data6/disk/liuyao/nascentRNA_tissue/splicingsite_FSR/exon_fsr_se/02.exondata_fsr_se/%s_blh_fsr_se_add.txt",spe),
                         header = T,stringsAsFactors = F);c_blh_fsr_se<-c_blh_fsr_se[,c(1,2,15,16,3:5,13,14,6:8,11,12,9:10)]

mrc_one2one_intron_all<-read.table("/mnt/data6/disk/liuyao/nascentRNA_tissue/exon_mapping/mouse_rat_cavia_one2one_intron.txt",
                                   header = T,stringsAsFactors = F)

mrc_one2one_intron<-mrc_one2one_intron_all%>%dplyr::select("mouse_intron","rat_intron","cavia_intron")
group<-c("up_5SS","down_5SS","down_3SS","up_3SS")
set.seed(12)
FSR_SE_match_randomRate<-as.data.frame(matrix(nrow=0,ncol=5))
colnames(FSR_SE_match_randomRate)<-c("M_R","M_C","R_C","M_R_C","location")
for (x in c(1,2,3,4)){
  m_blh_se_fsr_ss<-na.omit(m_blh_fsr_se[,c(1,x+1,6,x+6,11,x+11,16)]);colnames(m_blh_se_fsr_ss)[1]<-"mouse_intron"
  r_blh_se_fsr_ss<-na.omit(r_blh_fsr_se[,c(1,x+1,6,x+6,11,x+11,16)]);colnames(r_blh_se_fsr_ss)[1]<-"rat_intron"
  c_blh_se_fsr_ss<-na.omit(c_blh_fsr_se[,c(1,x+1,6,x+6,11,x+11,16)]);colnames(c_blh_se_fsr_ss)[1]<-"cavia_intron"
  
  m_blh_se_fsr_ss<-merge(mrc_one2one_intron,m_blh_se_fsr_ss,by="mouse_intron")
  mr_blh_se_fsr_ss<-merge(m_blh_se_fsr_ss,r_blh_se_fsr_ss,by="rat_intron")
  mc_blh_se_fsr_ss<-merge(m_blh_se_fsr_ss,c_blh_se_fsr_ss,by="cavia_intron")
  
  rc_blh_se_fsr_ss<-merge(mrc_one2one_intron,r_blh_se_fsr_ss,by="rat_intron")
  rc_blh_se_fsr_ss<-merge(rc_blh_se_fsr_ss,c_blh_se_fsr_ss,by="cavia_intron")
  
  random_100<-pbapply::pblapply(cl=50,1:100,function(r){
    mr_blh_se_fsr_ss_r<-random_percol(mr_blh_se_fsr_ss)
    mc_blh_se_fsr_ss_r<-random_percol(mc_blh_se_fsr_ss)
    rc_blh_se_fsr_ss_r<-random_percol(rc_blh_se_fsr_ss)
    
    mr_FSR_SE_match_introns<-pairorgans_FSR_SE_match_introns(mr_blh_se_fsr_ss_r,x)
    mc_FSR_SE_match_introns<-pairorgans_FSR_SE_match_introns(mc_blh_se_fsr_ss_r,x)
    rc_FSR_SE_match_introns<-pairorgans_FSR_SE_match_introns(rc_blh_se_fsr_ss_r,x)
    
    counts<-data.frame(SE_match=c(mr_FSR_SE_match_introns[[1]],mc_FSR_SE_match_introns[[1]],
                                  rc_FSR_SE_match_introns[[1]]),
                       FSR_SE_match=c(mr_FSR_SE_match_introns[[2]],mc_FSR_SE_match_introns[[2]],
                                      rc_FSR_SE_match_introns[[2]]))
    
    overlap_SE_match_introns<-intersect(mr_FSR_SE_match_introns[[3]],mc_FSR_SE_match_introns[[3]])
    overlap_SE_match_introns<-intersect(overlap_SE_match_introns,rc_FSR_SE_match_introns[[3]])
    overlap_FSR_SE_match_introns<-intersect(mr_FSR_SE_match_introns[[4]],mc_FSR_SE_match_introns[[4]])
    overlap_FSR_SE_match_introns<-intersect(overlap_FSR_SE_match_introns,rc_FSR_SE_match_introns[[4]])
    counts[4,]<-c(length(overlap_SE_match_introns),length(overlap_FSR_SE_match_introns))
    
    counts$ratio=counts$FSR_SE_match/counts$SE_match
    
    data<-data.frame(M_R=counts$ratio[1],M_C=counts$ratio[2],
                     R_C=counts$ratio[3],M_R_C=counts$ratio[4])
  })%>%rbind.fill()
  random_100<-random_100%>%mutate(location=group[x])
  FSR_SE_match_randomRate<-rbind(FSR_SE_match_randomRate,random_100)
}

save(FSR_SE_match_randomRate,file = "/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_0.5/Rdata/FSR_SE_match_random100_seed12.Rdata")

all_data_long <- FSR_SE_match_randomRate %>%
  pivot_longer(cols = c("M_R", "M_C", "R_C", "M_R_C"),
               names_to = "Group", values_to = "Ratio") %>%select(location, Group, Ratio)
all_data_long$location <- factor(all_data_long$location, levels = c("up_5SS", "down_5SS", "up_3SS", "down_3SS"))
all_data_long$Group <- factor(all_data_long$Group, levels = c("M_R", "M_C", "R_C", "M_R_C"))

fig_5b <- ggplot(all_data_long, aes(x = Group, y = Ratio, fill = Group)) +
  geom_violin(alpha = 0.5, trim = FALSE,linewidth = 0.5) +
  geom_boxplot(alpha = 0.3, width = 0.1, outlier.shape =NA,linewidth = 0.25) +
  geom_segment(data = real_data,aes(x = as.numeric(factor(Group)) - 0.08, 
               xend = as.numeric(factor(Group)) + 0.08, y = ratio, yend = ratio),
               color = "red", alpha = 0.8, size = 1, lineend = "round") +
  scale_fill_manual(values = c("#B1DDF4", "#E3C0DA", "#FFCC99", "#FF9999")) +
  coord_cartesian(ylim = c(0,0.6)) +
  facet_wrap(~ location, ncol = 2, scales = "free_x",
             labeller = labeller(location = c("up_5SS" = "50nts upstream of 5'SS",
                                              "down_5SS" = "50nts downstream of 5'SS",
                                              "up_3SS" = "50nts upstream of 3'SS",
                                              "down_3SS" = "50nts downstream of 3'SS")))+theme_bw() +
  labs(x = NULL, y = "Proportion of introns with concordant FSR trends \nin SE-conserved introns among species") +
  theme(strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 8),legend.position = "none",
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        axis.title = element_text(size = 9),axis.text = element_text(size = 8) )
save(all_data_long,real_data,fig_5b,file = "/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_0.5/Rdata/FSR_SE_match_random100_fig_seed12.Rdata" )
ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_0.5/fig_5b_random100_col1.pdf",
       plot = fig_5b,width = 16, height = 14,units = "cm")

