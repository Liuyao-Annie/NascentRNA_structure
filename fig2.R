
####fig2a
setwd("/mnt/data6/disk/liuyao/nascentRNA_tissue/correlation/")
result=as.data.frame(matrix(nrow = 0, ncol = 5))
colnames(result)<-c("order","pearson_r","Significance","group","species")
for (spe in c("M","R","C")){
  if (spe=="M"){spe1="Mouse"}else if (spe=="R"){spe1="Rat"}else{spe1="Cavia"}
  load(sprintf("%s_tissus_overlapgene.Rdata",spe))
  specie_lh<-read.table(sprintf("%s_lh_gene_pearson_o.txt",spe),header = T,stringsAsFactors = F)
  specie_bh<-read.table(sprintf("%s_bh_gene_pearson_o.txt",spe),header = T,stringsAsFactors = F)
  specie_bl<-read.table(sprintf("%s_bl_gene_pearson_o.txt",spe),header = T,stringsAsFactors = F)
  specie_lh_bh<-merge(specie_lh[,1:3],specie_bh[,1:3],by="genename")
  specie_lh_bh_bl<-merge(specie_lh_bh,specie_bl[,1:3],by="genename")
  specie_lh_bh_bl<-specie_lh_bh_bl%>%
    mutate(flag_lh=ifelse(pvalue_l_h<0.05,"P < 0.05","P >= 0.05"),
           flag_bh=ifelse(pvalue_b_h<0.05,"P < 0.05","P >= 0.05"),
           flag_bl=ifelse(pvalue_b_l<0.05,"P < 0.05","P >= 0.05"),
           mean_r=(pearson_l_h+pearson_b_h+pearson_b_l)/3)
  
  specie_lh_bh_bl<-specie_lh_bh_bl%>%
    mutate(flag_blh=ifelse(flag_lh=="P < 0.05" & flag_bh=="P < 0.05" & flag_bl=="P < 0.05" ,"P < 0.05","P >= 0.05"))
  specie_lh_bh_bl<-specie_lh_bh_bl[order(specie_lh_bh_bl$mean_r),]
  specie_lh_bh_bl$order<-1:nrow(specie_lh_bh_bl)
  
  data1<-specie_lh_bh_bl%>%dplyr::select(c("order","mean_r","flag_blh"))%>%mutate(group="Averaged similarity")
  colnames(data1)<-c("order","pearson_r","Significance","group")
  data2<-specie_lh_bh_bl%>%dplyr::select(c("order","pearson_l_h","flag_lh"))%>%mutate(group="Liver-Heart")
  colnames(data2)<-c("order","pearson_r","Significance","group")
  data3<-specie_lh_bh_bl%>%dplyr::select(c("order","pearson_b_h","flag_bh"))%>%mutate(group="Brain-Heart")
  colnames(data3)<-c("order","pearson_r","Significance","group")
  data4<-specie_lh_bh_bl%>%dplyr::select(c("order","pearson_b_l","flag_bl"))%>%mutate(group="Brain-Liver")
  colnames(data4)<-c("order","pearson_r","Significance","group")
  data<-rbind(data1,data2,data3,data4)
  data<-data%>%mutate(species=spe1)
  result<-rbind(result,data)
}
result$species<-factor(result$species,levels = c("Mouse","Rat","Cavia"))
fig2a<-ggplot(result, aes(order,pearson_r,fill = Significance)) +geom_bar(stat = "identity")+
  facet_grid(rows = vars(group), cols = vars(species), scales = "free_x") +
  scale_fill_manual(values = c("red","grey"))+
  ylab("Inter-organ structural similarity (Pearson's correlation coefficient of FSR)" )+
  xlab("Genes")+theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=8),panel.grid = element_blank(),
        axis.title = element_text(size=8),axis.text = element_text(size=8),
        legend.text = element_text(size=8),legend.title = element_text(size=8),
        strip.text = element_text(size = 8),legend.position = "none")

ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_0.5/fig2a.pdf", 
       plot = fig2a,width =19 , height =12 , device = cairo_pdf,units = "cm")


###fig2b
###mouse_cds_phscore
chrid=1
allgene_cds_phscore=as.data.frame(matrix(ncol = 4,nrow = 0))
colnames(allgene_cds_phscore)<-c("genename","average_score","row_count","CDS_len")
for (chrid in c(1:19,"X","Y")){
  da<-read.table(sprintf("/mnt/data6/disk/huangyanying/comparison/CDS/phscore_length/chr%s_phscore_length.txt",chrid),
                 header = T,stringsAsFactors = F)
  allgene_cds_phscore<-rbind(allgene_cds_phscore,da)
}

allgene_cds_phscore<-na.omit(allgene_cds_phscore)
allgene_cds_phscore<-allgene_cds_phscore%>%mutate(coverage=row_count/CDS_len)
allgene_cds_phscore_f<-allgene_cds_phscore%>%dplyr::filter(coverage>=0.5)
allgene_cds_phscore_f<-allgene_cds_phscore_f[,c(1,2)]
colnames(allgene_cds_phscore_f)<-c("genename","phscore")

###mouse_cds_fsr_pearson  averaged inter-organ structural similarities fig2b-1
setwd("/mnt/data6/disk/liuyao/nascentRNA_tissue/correlation/")
spe="M"
mouse_lh<-read.table(sprintf("./cds_tissues_p/%s_lh_cds_pearson_o.txt",spe),header = T,stringsAsFactors = F)
mouse_bh<-read.table(sprintf("./cds_tissues_p/%s_bh_cds_pearson_o.txt",spe),header = T,stringsAsFactors = F)
mouse_bl<-read.table(sprintf("./cds_tissues_p/%s_bl_cds_pearson_o.txt",spe),header = T,stringsAsFactors = F)

mouse_lh_bh<-merge(mouse_lh[,1:3],mouse_bh[,1:3],by="genename")
mouse_lh_bh_bl<-merge(mouse_lh_bh,mouse_bl[,1:3],by="genename")
mouse_lh_bh_bl<-mouse_lh_bh_bl%>%
  mutate(flag_lh=ifelse(pvalue_l_h<0.05,"P < 0.05","P >= 0.05"),
         flag_bh=ifelse(pvalue_b_h<0.05,"P < 0.05","P >= 0.05"),
         flag_bl=ifelse(pvalue_b_l<0.05,"P < 0.05","P >= 0.05"),
         pearson_r=(pearson_l_h+pearson_b_h+pearson_b_l)/3)

mouse_mean_pearson<-mouse_lh_bh_bl%>%dplyr::select(c("genename","pearson_r"))
mouse_cdsphastcos_r<-merge(mouse_mean_pearson,allgene_cds_phscore_f,by="genename")

cor.test(mouse_cdsphastcos_r$pearson_r,mouse_cdsphastcos_r$phscore,method="s")
##r=0.141702  p.value=6.136e-04
mouse_phastcos_r_p<-ggplot(mouse_cdsphastcos_r, aes(x =phscore, y = pearson_r)) + geom_point(size=0.5) +  
  labs(title = "Mouse",x = "Evolutionary conservation (phastCons)", 
       y = "Inter-organ structural similarity")+
  annotate("text",x=0.25,y=0.5,size=2.5,label=expression(italic(rho)* " = 0.14  " * italic(P) * " < " * 10^{-3}))+theme_bw()+
  theme(plot.background = element_rect(fill = "white"),panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5,size=8),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.title = element_text(size=8),
        axis.text = element_text(size=8),legend.text = element_text(size=8),
        legend.title = element_text(size=8))
save(mouse_cdsphastcos_r,mouse_phastcos_r_p,
     file = "/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/mouse_phastCons.Rdata")
ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_new/fig2b_mouse.pdf", 
       plot = mouse_phastcos_r_p,width =6.5 ,height =7 ,device = cairo_pdf , units = "cm")

#####rat_cds_phscore
chrid=1
allgene_cds_phscore=as.data.frame(matrix(ncol = 4,nrow = 0))
colnames(allgene_cds_phscore)<-c("genename","average_score","row_count","CDS_len")
for (chrid in c(1:19,"X","Y")){
  da<-read.table(sprintf("/mnt/data6/disk/huangyanying/comparison/Rat_CDS/CDS_phscore_length/chr%s_phscore_length.txt",chrid),
                 header = T,stringsAsFactors = F)
  allgene_cds_phscore<-rbind(allgene_cds_phscore,da)
}
allgene_cds_phscore<-na.omit(allgene_cds_phscore)
allgene_cds_phscore<-allgene_cds_phscore%>%mutate(coverage=row_count/CDS_len)
allgene_cds_phscore_f<-allgene_cds_phscore%>%dplyr::filter(coverage>=0.5)
allgene_cds_phscore_f<-allgene_cds_phscore_f[,c(1,2)]
colnames(allgene_cds_phscore_f)<-c("genename","phscore")


###rat_cds_fsr_pearson  averaged inter-organ structural similarities fig2b-2
setwd("/mnt/data6/disk/liuyao/nascentRNA_tissue/correlation/")
spe="R"
rat_lh<-read.table(sprintf("./cds_tissues_p/%s_lh_cds_pearson_o.txt",spe),header = T,stringsAsFactors = F)
rat_bh<-read.table(sprintf("./cds_tissues_p/%s_bh_cds_pearson_o.txt",spe),header = T,stringsAsFactors = F)
rat_bl<-read.table(sprintf("./cds_tissues_p/%s_bl_cds_pearson_o.txt",spe),header = T,stringsAsFactors = F)

rat_lh_bh<-merge(rat_lh[,1:3],rat_bh[,1:3],by="genename")
rat_lh_bh_bl<-merge(rat_lh_bh,rat_bl[,1:3],by="genename")
rat_lh_bh_bl<-rat_lh_bh_bl%>%
  mutate(flag_lh=ifelse(pvalue_l_h<0.05,"P < 0.05","P >= 0.05"),
         flag_bh=ifelse(pvalue_b_h<0.05,"P < 0.05","P >= 0.05"),
         flag_bl=ifelse(pvalue_b_l<0.05,"P < 0.05","P >= 0.05"),
         pearson_r=(pearson_l_h+pearson_b_h+pearson_b_l)/3)

rat_mean_pearson<-rat_lh_bh_bl%>%dplyr::select(c("genename","pearson_r"))
rat_cdsphastcos_r<-merge(rat_mean_pearson,allgene_cds_phscore_f,by="genename")
cor.test(rat_cdsphastcos_r$pearson_r,rat_cdsphastcos_r$phscore,method="s")$p.value
##r=0.2155594  p.value=1.983e-07

rat_phastcos_r_p<-ggplot(rat_cdsphastcos_r, aes(x =phscore, y = pearson_r)) + geom_point(size=0.5) +  
  labs(title = "Rat",x = "Evolutionary conservation (phastCons)", 
       y = "Inter-organ structural similarity")+
  scale_y_continuous(limits = c(0.1,0.7),breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7))+
  annotate("text",x=0.22,y=0.65,size=2.5,label=expression(italic(rho)* " = 0.22  " * italic(P) * " < " * 10^{-6}))+theme_bw()+
  theme(plot.background = element_rect(fill = "white"),panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5,size=8),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.title = element_text(size=8),
        axis.text = element_text(size=8),legend.text = element_text(size=8),
        legend.title = element_text(size=8))
save(rat_cdsphastcos_r,rat_phastcos_r_p,
     file = "/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/rat_phastCons.Rdata")

ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_new/fig2b_rat.pdf", 
       plot = rat_phastcos_r_p,width =6.5 ,height =7 ,device = cairo_pdf , units = "cm")

####fig2b cavia_cds_phscore
cavia_phastcos_score<-read.table("/mnt/data5/disk/liuyao/cavia_data/score/cavia_allexon_phastcos_17way_ratmod.txt",header = F)
colnames(cavia_phastcos_score)<-c("chrM","genename","CDS_startsite","CDS_stopsite","exon_len",
                                  "scaffold_name","mean_phastcos","phastcos_siteCounts")
cavia_phastcos_score_f<-na.omit(cavia_phastcos_score)

allgene_phastcos_score<-cavia_phastcos_score_f%>%dplyr::group_by(genename)%>%
  dplyr::summarise(phastcos_score=mean(mean_phastcos),
            CDS_len=sum(exon_len),all_phastcos_Counts=sum(phastcos_siteCounts))

allgene_phastcos_score<-allgene_phastcos_score%>%mutate(coverage=all_phastcos_Counts/CDS_len)
allgene_phastcos_score_f<-allgene_phastcos_score%>%dplyr::filter(coverage>=0.5)

setwd("/mnt/data6/disk/liuyao/nascentRNA_tissue/correlation/")
spe="C"
cavia_lh<-read.table(sprintf("./cds_tissues_p/%s_lh_cds_pearson_o.txt",spe),header = T,stringsAsFactors = F)
cavia_bh<-read.table(sprintf("./cds_tissues_p/%s_bh_cds_pearson_o.txt",spe),header = T,stringsAsFactors = F)
cavia_bl<-read.table(sprintf("./cds_tissues_p/%s_bl_cds_pearson_o.txt",spe),header = T,stringsAsFactors = F)

cavia_lh_bh<-merge(cavia_lh[,1:3],cavia_bh[,1:3],by="genename")
cavia_lh_bh_bl<-merge(cavia_lh_bh,cavia_bl[,1:3],by="genename")
cavia_lh_bh_bl<-cavia_lh_bh_bl%>%
  mutate(flag_lh=ifelse(pvalue_l_h<0.05,"P < 0.05","P >= 0.05"),
         flag_bh=ifelse(pvalue_b_h<0.05,"P < 0.05","P >= 0.05"),
         flag_bl=ifelse(pvalue_b_l<0.05,"P < 0.05","P >= 0.05"),
         pearson_r=(pearson_l_h+pearson_b_h+pearson_b_l)/3)
cavia_mean_pearson<-cavia_lh_bh_bl%>%dplyr::select(c("genename","pearson_r"))


cavia_cdsphastcos_r<-merge(cavia_mean_pearson,allgene_phastcos_score_f,by="genename")
cor.test(cavia_cdsphastcos_r$pearson_r,cavia_cdsphastcos_r$phastcos_score,method="s",exact = F)
##r=0.136048  p.value=0.000172

cavia_phastcos_r_p<-ggplot(cavia_cdsphastcos_r, aes(x = phastcos_score, y = pearson_r)) +
  geom_point(size = 0.5) +
  labs(title = "Cavia",x = "Evolutionary conservation (phastCons)",
       y = "Inter-organ structural similarity") +
  annotate("text",x = 0.24,y = 0.53,size = 2.5,
           label=expression(italic(rho)* " = 0.14  " * italic(P) * " < " * 10^{-3}))+theme_bw()+
  theme(plot.background = element_rect(fill = "white"),panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5,size=8),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.title = element_text(size=8),
        axis.text = element_text(size=8),legend.text = element_text(size=8),
        legend.title = element_text(size=8))
save(cavia_cdsphastcos_r,cavia_phastcos_r_p,
     file = "/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/cavia_phastCons.Rdata")
ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_0.5/fig2b_cavia.pdf", 
       plot = cavia_phastcos_r_p,width =6.5 ,height =7 ,device = cairo_pdf , units = "cm")

fig2b<-ggarrange(mouse_phastcos_r_p,rat_phastcos_r_p,cavia_phastcos_r_p,ncol=3)
ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_0.5/fig2b_3spe.pdf", 
       plot = fig2b,width =19 ,height =6.5 ,device = cairo_pdf , units = "cm")


###fig2c
###数据来源于"/mnt/data6/disk/liuyao/nascentRNA_tissue/conservation/mr_conservation_inter_organ_cor.R"
partial_correlation_forfig<-data.frame(rho=c(0.146,0.125,0.144,0.222,0.269,0.211,0.136,0.147,0.131),
                                       type=(rep(c("Conservation level~\nstructural similarity",
                                                   "Partial correlation controlling \nthe nascent RNA expression",
                                                   "Partial correlation controlling \nthe total RNA expression"),times=3)),
                                       Organ=rep(c("Mouse","Rat","Cavia"),each=3))

partial_correlation_forfig$Organ<-factor(partial_correlation_forfig$Organ,levels = c("Cavia","Rat","Mouse"))
partial_correlation_forfig$type<-factor(partial_correlation_forfig$type,
                                        levels = c("Partial correlation controlling \nthe total RNA expression",
                                                   "Partial correlation controlling \nthe nascent RNA expression",
                                                   "Conservation level~\nstructural similarity"))

fig2c<-ggplot(partial_correlation_forfig,aes(x = Organ, y = rho, fill = type)) +
  geom_col(width = 0.8, position= position_dodge(width=0.9))+
  scale_y_continuous(limits = c(0,0.37), breaks = c(0,0.1,0.2,0.3))+
  coord_flip() +labs(x = "", y = "Correlation coefficient(rho)",fill = "Group") +theme_bw()+
  scale_fill_manual(values = c("#bfd0e1","#Fcd7af","#efbcb8"))+
  theme(plot.background = element_rect(fill = "white"),panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5),panel.grid = element_blank(),
        axis.title = element_text(size=8),axis.text = element_text(size=8),
        legend.text = element_text(size=8),legend.title = element_text(size=8)) 

ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_0.5/fig2c_new2.pdf", 
       plot = fig2c,width = 10.5, height =5.2,units = "cm")

##fig2d
load("/mnt/data6/disk/liuyao/nascentRNA_tissue/correlation/heatmap_fsr/gene_fsr_forheatmap4.Rdata")
pearson_rms<-data.frame(matrix(nrow = 9,ncol = 9))
i=1
for (i in 1:9){
  aa<-round(as.numeric(unlist(str_split(forheatmap1[,i],pattern = "_"))),4)
  pearson_rms[,i]<-aa[1:36 %% 4 == 1]
}
sample<-c("M_B","M_H","M_L","R_B","R_H","R_L","C_B","C_H","C_L")
colnames(pearson_rms)<-sample;rownames(pearson_rms)<-sample

####聚类
hc<-hclust(dist(pearson_rms),method = "average") #对行进行聚类
rowInd<-hc$order #将聚类后行的顺序存为rowInd
hc<-hclust(dist(t(pearson_rms)),method = "average")  #对矩阵进行转置，对原本的列进行聚类
colInd<-hc$order  #将聚类后列的顺序存为colInd
pearson_rms<-pearson_rms[rowInd,colInd]
pearson_rms$Row<-rownames(pearson_rms) ####加一行方便melt
# 使用melt函数，将数据框转换为长格式##library(reshape2)
melted_data <- reshape2::melt(pearson_rms, id.vars = "Row", variable.name = "Column", value.name = "Correlation")
melted_data$Correlation[melted_data$Row == melted_data$Column] <- NA  ###把对角线值变成NA
sample_order <- pearson_rms$Row  
gene_order <- pearson_rms$Row 

# 绘制热图
p<-ggplot(melted_data, aes(x = Column, y = Row, fill = Correlation)) +geom_tile() +
  scale_fill_gradient2(low = "#FFFFB8", high = "red",mid = "#FF6C00", 
                       midpoint = 0.2, limit = c(0, 0.4), na.value = "white") +
  labs(x="",y="") +theme_minimal()+ 
  scale_x_discrete(limits = sample_order) +scale_y_discrete(limits = gene_order) +
  ggtitle("Structural similarities between various \ntypes of samples(FSR similarity)")+
  theme(plot.title = element_text(hjust = 0.5,size = 8),
        axis.title = element_text(size=8),axis.text = element_text(size=8),
        legend.text = element_text(size=8),legend.title = element_text(size=8)) 
# 绘制行聚类树并加到热图左侧
library(ggtree)
library(aplot)
row_tree <- ggtree(ape::as.phylo(hc), layout = "rectangular", branch.length = "none")+theme(legend.position = "none")
fig2d<-p%>%insert_left(row_tree,width = 0.1)# 将聚类树插入热图的左侧

ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_new/fig2d_heatmap.pdf",
       plot = fig2d,width = 10, height = 7.5,units = "cm")


##fig2e
load("/mnt/data6/disk/liuyao/nascentRNA_tissue/ViennaRNA/fig/Dot_Bracket_structure_similarity_forheatmap.Rdata")
structure_rms<-forheatmap_identical_structure_rms
rownames(structure_rms) <- gsub("heart", "H", gsub("liver", "L", gsub("brain", "B", rownames(structure_rms))))
colnames(structure_rms) <- gsub("heart", "H", gsub("liver", "L", gsub("brain", "B", colnames(structure_rms))))
####聚类
hc<-hclust(dist(structure_rms),method = "average") 
rowInd<-hc$order 
hc<-hclust(dist(t(structure_rms)),method = "average") 
colInd<-hc$order 
structure_rms<-structure_rms[rowInd,colInd]
structure_rms$Row<-rownames(structure_rms) 
melted_data <- reshape2::melt(structure_rms, id.vars = "Row", variable.name = "Column", value.name = "Proportion")
melted_data$Proportion[melted_data$Row == melted_data$Column] <- NA 
sample_order <- structure_rms$Row
gene_order <- structure_rms$Row 
# 绘制热图
gene_structure_p<-ggplot(melted_data, aes(x =Column , y = Row, fill = Proportion)) +
  geom_tile() +
  scale_fill_gradient2(low = "#FFFFB8", high = "red",mid = "#FF6C00",
                       midpoint = 0.45, limit = c(0.35,0.6), na.value = "white") + 
  labs(x="",y="") +theme_minimal()+ 
  scale_x_discrete(limits = sample_order) +scale_y_discrete(limits = gene_order) +##加聚类信息
  ggtitle("Structural similarities between various types \nof samples (Dot-Bracket_consistency)")+
  theme(plot.title = element_text(hjust = 0.5,size=8),
        axis.title = element_text(size=8),axis.text = element_text(size=8),
        legend.text = element_text(size=8),legend.title = element_text(size=8)) 
# 绘制行聚类树并加到热图左侧
row_tree <- ggtree(ape::as.phylo(hc), layout = "rectangular", branch.length = "none") 
fig2e<-gene_structure_p%>%insert_left(row_tree,width = 0.1)
ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_new/fig2e_heatmap.pdf", 
       plot = fig2e,width = 10, height = 7.5,units = "cm")

