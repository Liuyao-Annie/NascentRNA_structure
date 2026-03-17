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
####fig1b-d
load("/mnt/data6/disk/huangyanying/mouse_plot/04_pea/02conbine_pea/all_C.Rdata")
data<-data_C;pearson<-pearson_C
##换行名列名格式
result <- gsub("brain", "B", colnames(data))            # 替换 "brain" → "B"
result <- gsub("heart", "H", result) 
result <- gsub("liver", "L", result) 
result <- gsub("(D)(\\d)", "\\1-\\2", result)  # D1 → D-1（正则捕获组）
result <- gsub("(N)(\\d)", "\\1-\\2", result) 
#colnames(data)<-result
colnames(pearson)<-result[-1]
rownames(pearson)<-result[-1]
####聚类
hc<-hclust(dist(pearson),method = "average") #对行进行聚类
rowInd<-hc$order #将聚类后行的顺序存为rowInd
hc<-hclust(dist(t(pearson)),method = "average")  #对矩阵进行转置，对原本的列进行聚类
colInd<-hc$order  #将聚类后列的顺序存为colInd
pearson<-pearson[rowInd,colInd]
row=colnames(pearson)
col=row

pearson<-pearson[row,col]
# 绘制热图  
# 假设你有一个矩阵mat，表示热图的数据  
# 你可能需要先将这个矩阵转换为长格式的数据框，以便ggplot2可以使用  
# 例如，使用reshape2包中的melt函数
#library(reshape2)
df_long <- reshape2::melt(pearson)  
names(df_long) <- c("Samples", "Col", "value")  

hotmap <- function(df_long,start,mid,mid_value) {
  hotmap_fig<-ggplot(df_long, aes(x = Col, y = Samples, fill = value)) +  
    geom_tile() + 
    scale_fill_gradientn(colors = c("blue","white", "red"), # 定义颜色从白色到红色  
                                        breaks = c(start,mid, 1),          # 定义颜色映射的断点  
                                        labels = c(start,mid, "1"),      # 为断点提供标签（在图例中显示）  
                                        limits = c(start, 1),          # 设置颜色映射的范围  
                                        values = c(0,mid_value,1),
                                        na.value = "grey80") +  theme_minimal() +    # 如果数据中有NA值，设置一个默认颜色  
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5,size=8),
          legend.direction = "horizontal",legend.position = "top",legend.margin = margin(t = 10, r = 5, b = 0.01, l = 5),
          legend.key.width = unit(0.6, "cm"),legend.key.height = unit(0.25, "cm"),
          axis.title = element_text(size=8),axis.text = element_text(size=5),
          legend.text = element_text(size=8),legend.title = element_text(size=8)) + 
    coord_fixed() +  labs(x = "Samples", y = "Samples")+
    ggtitle("Pearson's correlation coefficient for the\n number of reads mapped to each genes")
  return(hotmap_fig)
}

start=0.5;mid=0.85;mid_value=0.7
df_long_C<-df_long
start=0.3;mid=0.65;mid_value=0.5
# 使用ggplot2绘制热图  
hotmap_C_new<-hotmap(df_long_C,0.3,0.65,0.5)

load("/mnt/data6/disk/huangyanying/update_M_heart/MR_counts_heatmap.Rdata")
hotmap_M_new<-hotmap(df_long_M,0.4,0.75,0.6)
hotmap_R_new<-hotmap(df_long_R,0.5,0.85,0.7)
save(df_long_M,df_long_C,df_long_R,hotmap_M_new,hotmap_R_new,hotmap_C_new,
     file ="/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig1b_d_MCR_counts_heatmap.Rdata" )
load("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig1b_d_MCR_counts_heatmap.Rdata")
hotmap_MCR_counts<-ggarrange(hotmap_M_new,hotmap_R_new,hotmap_C_new,ncol=3)
# ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/MCR_counts.pdf", plot = hotmap_MCR_counts,
#        width = 24, height = 8)
ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_new/fig1b_d_MCR_counts.pdf", plot = hotmap_MCR_counts,
       width = 20, height = 8,units = "cm")
       




###fig1b inside 
load("/mnt/data6/disk/huangyanying/mouse_plot/04_pea/02conbine_pea/new_all_M.Rdata")
new_df <-data_M %>% dplyr::select(`M-heart-D1`, `M-heart-D2`) # 或者使用 df[, c("A", "B")]
new_df$`M-heart-D1` <- log10(new_df$`M-heart-D1` + 1)  
new_df$`M-heart-D2` <- log10(new_df$`M-heart-D2` + 1)
new_df <- na.omit(new_df)  
# 绘制散点图  
p_value=0.999
# 计算一个合适的位置来放置 p 值文本  
# 这里我们简单地使用数据范围的一个百分比作为位置，但您可能需要根据实际情况调整  
x_pos <- max(new_df$`M-heart-D1`) * 0.7  # x 位置稍微超出数据的最大值  
y_pos <- max(new_df$`M-heart-D2`) * 0.95  # y 位置在数据的最大值下方一点  
pp=ggplot(new_df, aes(x = `M-heart-D1`, y = `M-heart-D2`)) +  geom_point(size=0.01) +  
  scale_x_continuous(breaks = c(0, 1, 3, 5, 7),  # 设置刻度位置10¹, "10²", "10⁴", "10⁶", "10⁸"
                     labels = c("0", "10¹", "10³", "10⁵", "10⁷"),limits = c(0, 7))+    # 设置刻度标签
  scale_y_continuous(breaks = c(0, 1, 3, 5, 7),  # 设置刻度位置10¹, "10²", "10⁴", "10⁶", "10⁸"
                     labels = c("0", "10¹", "10³", "10⁵", "10⁷"),limits = c(0, 7))+    # 设置刻度标签
  theme(plot.background = element_rect(fill = NA),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA))+  # 添加边框，填充色设置为 NA（因为我们只想要边框）)+
  labs(title = "",x = "Number of reads", y = "Number of reads")+
  theme(plot.title = element_text(hjust = 0.5),#panel.border = element_blank(),
        axis.ticks = element_line(linewidth = 0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks.length = unit(0.04, "cm"),
        axis.title = element_text(size=5.5),axis.text = element_text(size=5.5),
        legend.text = element_text(size=5.5),legend.title = element_text(size=5.5))
ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_new/fig1b_inside.pdf", 
       plot = pp,width = 2.5, height = 3, device = cairo_pdf,units = "cm")


###fig 1e
###gene
##横坐标为average numbers of reads of the gene，纵坐标为两个样本之间每个基因的FSR的P值
pearson1=read.table("/mnt/data6/disk/huangyanying/mouse_plot/01_rawdata/filtersite/C/C-brain-N_union.txt",header = T)
density1=read.table("/mnt/data6/disk/huangyanying/mouse_plot/FSR02/filter_density/C-brain-N_freq1.txt",header = T)
colnames(density1)=c("genename","freq","genelength","exon","intron","density")
pearson2=merge(pearson1,density1,by="genename",all = T)
reps_pearson2<-pearson2[order(pearson2$density),]
reps_pearson2_1<-na.omit(reps_pearson2)
fig1e<-ggplot(reps_pearson2_1,aes(density,rep1_rep2))+#geom_point(size=2,alpha=0.4)+
  scale_x_log10(limits=c(1,1000),breaks = c(1,10,100,1000),
                labels =as.character(c(1,10,100,1000)))+
  scale_y_continuous(limits = c(0,1),breaks = c(0,0.25,0.50,0.75,1),
                     labels =as.character(c(0,0.25,0.50,0.75,1)))+
  labs(x="Reads coverage of the genes",
       y="Pearson’s correlation of RT stop counts between \nbiological repeats of NAI-N3 treated Cavia brain",parse=TRUE)+
  theme_bw()+geom_pointdensity(size=0.8,alpha=0.6,adjust = 0.5)+
  geom_point(data = reps_pearson2_1[4000, ],color = "red", size = 1,alpha=1) + #例子标红（U5 snRNA）
  theme(plot.title = element_text(hjust = 0.5),##panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=8),axis.text = element_text(size=8),
        legend.text = element_text(size=8),legend.title = element_text(size=8),
        legend.position="none")
ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_new/fig1e.pdf", 
       plot = fig1e,width = 7, height = 6.5,units = "cm") 

###fig 1f
#######    C_brain_U5_RT_stopsites_reps                                                                                                                                  

C_brain_U5_rep1<-read.table("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig/C_brain_U5/C_brain_U5_rep1.txt",
                            header = F,stringsAsFactors = F)
C_brain_U5_rep2<-read.table("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig/C_brain_U5/C_brain_U5_rep2.txt",
                            header = F,stringsAsFactors = F)
C_brain_U5_rep3<-read.table("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig/C_brain_U5/C_brain_U5_rep3.txt",
                            header = F,stringsAsFactors = F)

U5_C_brain_rep1<-as.data.frame(table(C_brain_U5_rep1$V2))
U5_C_brain_rep2<-as.data.frame(table(C_brain_U5_rep2$V2))
U5_C_brain_rep3<-as.data.frame(table(C_brain_U5_rep3$V2))
U5_C_brain_rep1$Var1<-as.numeric(as.character(U5_C_brain_rep1$Var1))
U5_C_brain_rep2$Var1<-as.numeric(as.character(U5_C_brain_rep2$Var1))
U5_C_brain_rep3$Var1<-as.numeric(as.character(U5_C_brain_rep3$Var1))

add1<-data.frame(Var1=setdiff(1:116,U5_C_brain_rep1$Var1),Freq=0)
add2<-data.frame(Var1=setdiff(1:116,U5_C_brain_rep2$Var1),Freq=0)
add3<-data.frame(Var1=setdiff(1:116,U5_C_brain_rep3$Var1),Freq=0)

U5_C_brain_rep1_add<-rbind(U5_C_brain_rep1,add1)
U5_C_brain_rep2_add<-rbind(U5_C_brain_rep2,add2)
U5_C_brain_rep3_add<-rbind(U5_C_brain_rep3,add3)

U5_C_brain_rep1_add<-U5_C_brain_rep1_add%>%dplyr::mutate(type="rep1")
U5_C_brain_rep2_add<-U5_C_brain_rep2_add%>%dplyr::mutate(type="rep2")
U5_C_brain_rep3_add<-U5_C_brain_rep3_add%>%dplyr::mutate(type="rep3")

rep1_rep2_rep3<-rbind(U5_C_brain_rep1_add,U5_C_brain_rep2_add,U5_C_brain_rep3_add)

rep1_rep2_rep3$Var1 <- as.numeric(as.character(rep1_rep2_rep3$Var1))
save(rep1_rep2_rep3,file = "/mnt/data6/disk/liuyao/nascentRNA_tissue/fig/C_brain_U5/C_brain_U5_RTStops_reps.Rdata")

##add D1
C_brain_U5_D1<-read.table("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig/C_brain_U5/C_brain_D1.txt",
                          header = F,stringsAsFactors = F)
U5_C_brain_D1<-as.data.frame(table(C_brain_U5_D1$V2))
U5_C_brain_D1$Var1<-as.numeric(as.character(U5_C_brain_D1$Var1))
add<-data.frame(Var1=setdiff(1:116,U5_C_brain_D1$Var1),Freq=0)
U5_C_brain_D1_add<-rbind(U5_C_brain_D1,add)
U5_C_brain_D1_add<-U5_C_brain_D1_add%>%dplyr::mutate(type="D1")
rep1_rep2_rep3_D1<-rbind(rep1_rep2_rep3,U5_C_brain_D1_add)

plot_data<-rep1_rep2_rep3_D1%>%dplyr::filter(type%in%c("rep1","rep2"))
C_brain_U5_rep1_rep2<-ggplot(rep1_rep2_rep3_D1%>%dplyr::filter(type%in%c("rep1","rep2")),aes(x=Var1,y=Freq,fill=type))+
  geom_bar(stat = "identity")+
  facet_grid(type~.) +ylab("RT stops")+xlab("Nucleotide positions(nt) of U5 snRNA")+
  scale_x_continuous(breaks = c(1,20,40,60,80,100,116))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),  # 白色背景 + 黑色边框
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),          # 面板边框（可选）
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.title = element_text(size=8),axis.text = element_text(size=8),
        legend.text = element_text(size=8),legend.title = element_text(size=8),
        strip.background = element_rect(fill = "white", color = "black"),#,       # 背景颜色 # 边框颜色
        strip.text = element_text(size=8),legend.position = "none")
ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_new/fig1f.pdf", 
       plot = C_brain_U5_rep1_rep2,width = 7, height = 5,units = "cm")


library(magick)
image_read("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_new/1g.PNG") %>% 
  image_write("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_new/fig1g.pdf", format = "pdf")
library(png)
img <- readPNG("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_new/1g.PNG")
plot<-ggplot() +
  annotation_raster(
    img,
    xmin = 0, xmax = 1,  # 根据需求调整
    ymin = 0, ymax = 1 * (nrow(img) / ncol(img))  # 按图片比例计算高度
  ) +
  coord_fixed() +  # 锁定x/y轴比例
  theme_void()  
# 导出 PDF
ggsave("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/fig_new/fig1g.pdf", 
       plot, device = "pdf", width = 10, height = 10, units = "cm")

install.packages("rsvg")
rsvg::rsvg_pdf("output.svg", "output.pdf")
# ####fig1b-d
# load("/mnt/data6/disk/huangyanying/mouse_plot/04_pea/02conbine_pea/new_all_M.Rdata")
# data<-data_M;pearson<-pearson_M
# load("/mnt/data6/disk/huangyanying/mouse_plot/04_pea/02conbine_pea/all_R.Rdata")
# data<-data_R;pearson<-pearson_R
# load("/mnt/data6/disk/huangyanying/mouse_plot/04_pea/02conbine_pea/all_C.Rdata")
# data<-data_C;pearson<-pearson_C
# ##换行名列名格式
# result <- gsub("brain", "B", colnames(data))            # 替换 "brain" → "B"
# result <- gsub("heart", "H", result) 
# result <- gsub("liver", "L", result) 
# result <- gsub("(D)(\\d)", "\\1-\\2", result)  # D1 → D-1（正则捕获组）
# result <- gsub("(N)(\\d)", "\\1-\\2", result) 
# #colnames(data)<-result
# colnames(pearson)<-result[-1]
# rownames(pearson)<-result[-1]
# ####聚类
# hc<-hclust(dist(pearson),method = "average") #对行进行聚类
# rowInd<-hc$order #将聚类后行的顺序存为rowInd
# hc<-hclust(dist(t(pearson)),method = "average")  #对矩阵进行转置，对原本的列进行聚类
# colInd<-hc$order  #将聚类后列的顺序存为colInd
# pearson<-pearson[rowInd,colInd]
# row=colnames(pearson)
# col=row
# 
# pearson<-pearson[row,col]
# # 绘制热图  
# # 假设你有一个矩阵mat，表示热图的数据  
# # 你可能需要先将这个矩阵转换为长格式的数据框，以便ggplot2可以使用  
# # 例如，使用reshape2包中的melt函数
# #library(reshape2)
# df_long <- reshape2::melt(pearson)  
# names(df_long) <- c("Samples", "Col", "value")  
# 
# df_long_M<-df_long
# start=0.4;mid=0.75;mid_value=0.6
# df_long_R<-df_long
# start=0.5;mid=0.8;mid_value=0.6
# df_long_C<-df_long
# start=0.3;mid=0.65;mid_value=0.5
# # 使用ggplot2绘制热图  
# hotmap_C_new<-ggplot(df_long_C, aes(x = Col, y = Samples, fill = value)) +  
#   geom_tile() +  scale_fill_gradientn(colors = c("blue","white", "red"), # 定义颜色从白色到红色  
#                                       breaks = c(start,mid, 1),          # 定义颜色映射的断点  
#                                       labels = c(start,mid, "1"),      # 为断点提供标签（在图例中显示）  
#                                       limits = c(start, 1),          # 设置颜色映射的范围  
#                                       values = c(0,mid_value,1),
#                                       na.value = "grey80") +     # 如果数据中有NA值，设置一个默认颜色  
#   theme_minimal() +  
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         plot.title = element_text(hjust = 0.5,size=8),
#         legend.direction = "horizontal",
#         legend.position = "top",
#         legend.margin = margin(t = 10, r = 5, b = 0.01, l = 5),
#         legend.key.width = unit(0.6, "cm"),
#         #legend.justification = c(0.2, 0),  # 水平居中，垂直顶部对齐
#         legend.key.height = unit(0.25, "cm"),
#         axis.title = element_text(size=8),axis.text = element_text(size=8),
#         legend.text = element_text(size=8),legend.title = element_text(size=8)) + 
#   coord_fixed() +  labs(x = "Samples", y = "Samples")+
#   ggtitle("Pearson's correlation coefficient for the\n number of reads mapped to each genes")
# save(df_long_M,df_long_C,df_long_R,hotmap_M_new,hotmap_R_new,hotmap_C_new,
#      file ="/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/MCR_counts_heatmap.Rdata" )
# load("/mnt/data6/disk/liuyao/nascentRNA_tissue/fig2/MCR_counts_heatmap.Rdata")
# hotmap_MCR_counts<-ggarrange(hotmap_M_new,hotmap_R_new,hotmap_C_new,ncol=3)