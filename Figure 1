##################################Figure 1. Human lung cell atlas in the Chinese LungSCeQTL cohort######################################
R
library(data.table)
library(tidyverse)
library(Seurat) 
library(ggpubr)
library(RColorBrewer)
library(Matrix)
library(Rmagic)
library(qs) 
library(ggbreak)

lineageorder=c("Myeloid","Mast","TNKcell","Bcell","Epithelial","Stroma","Endothelial")
lineage_color <- c("#DC050C", "#F5BA3E", "#F6C141", "#F7EF5D", "#CAE0AB", "#1965B0", "#882E72")
cellorder21 <- c("MP", "CM", "Mono", "DC", "Neu",
                 "Mast",
                 "NKT", "CD8T", "CD4T",
                 "NaiveB", "Plasma",
                 "AT2", "AT1", "Club", "Ciliated", "Basal", 
                 "Meso","Fib", "MyoF",
                 "VECs", "LECs")
revcellorder21=rev(cellorder21)  
cellorder21_color <- c("#DC050C", "#E33B15", "#E96A1F", "#EF8829", "#F2A233", 
                       "#F5BA3E", 
                       "#F6C141", "#F6D047","#F6DF4E", 
                       "#F7EF5D", "#FFFF99", 
                       "#CAE0AB", "#A7D295", "#82C480", "#5BB66B", "#43AA4E",
                       "#0F3D82", "#1965B0", "#A6CEE3", 
                       "#882E72", "#D6C1DE")
revcellorder21_color=rev(cellorder21_color) 

sceqtl=qread('sceqtl_562255_anno.qs') #读取
Metadata <- fread("Metadata_sceqtl.xls")
head(Metadata)
dim(Metadata) # 562255

# Figure 1A Study Design###############################################################

# Figure 1B Number of cells per donor############################################
a <- data.frame(table(Metadata$ID))
pdf(file = "cells_perindividual_hist.pdf", width = 7.5, height = 7.5)
hist(a$Freq,
     xlab = "Individuals",
     ylab = "Cells per individual",
     main = "",
     col = "#1965B0"
)
dev.off()

b=subset(a,Freq>6000)
a$Freq1=ifelse(a$Freq>6000,6001,a$Freq)
summary(a$Freq1)
pdf(file = "cells_perindividual_hist1.pdf", width = 8, height = 8)
hist(a$Freq1,
     xlab = "Individuals",
     ylab = "Cells per individual",
     main = "",,border = "white",
     col = "#1965B0",
)
dev.off()

# Figure 1C Atlas of cell annotation############################################
# 21 celltype 
p=DimPlot(object=sceqtl,
        pt.size=0.01,
        group.by='subcelltype',
        raster=FALSE,
        cols=cellorder21_color, 
        label = T,
        repel = T)+  
  NoLegend()+  
  labs(title = "") 
pdf(file='UMAPplot_subcelltype.pdf',width=6,heigh=6)
print(p)
dev.off()
png(file='UMAPplot_subcelltype.png',width=15,heigh=15,units="cm",res=350)
print(p)
dev.off()

# Figure 1D Violin plot for maker gene expression profiles########################
library(MySeuratWrappers) 
# The previous marker list highly expressed genes of each cell type  
markers <- c('LYZ','MARCO','MKI67','FCN1','CD1C','FCGR3B','KIT',
'NKG7','CD3D','IL7R',
'MS4A1','MZB1',
'EPCAM','SFTPA1','AGER','SCGB1A1','TPPP3','KRT15',
'ITLN1','COL1A2','ACTA2',
'VWF','CCL21')
p=VlnPlot(sceqtl, features = markers,  group.by='subcelltype',
        stacked=T,pt.size=0,  
        cols = cellorder21_color,#颜色  
        direction = "horizontal",   #水平作图  
        x.lab = '', y.lab = '')+    #横纵轴不标记任何东西  
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())   #不显示坐标刻度 
pdf(file='Violinplot_subcelltype.pdf',width=10,heigh=6)
print(p)
dev.off()

# Table S3 List of significant differentially expressed genes from each cell type####
head(sceqtl@meta.data)
a=data.frame(table(sceqtl$orig.ident))
randomorig=as.character(sample(a$Var1,4)) #"NX28" "NX9"  "NX18" "NX29"
sceqtlsub=subset(sceqtl,orig.ident %in% randomorig)
rm(sceqtl)
logFCfilter=0.5
adjPvalFilter=0.05
sceqtlsub <- sceqtlsub %>% dplyr::glimpse()
Idents(sceqtlsub) <-sceqtlsub@meta.data$subcelltype
levels(sceqtlsub) 
markers <- FindAllMarkers(object=sceqtlsub,only.pos=F,min.pct=0.25,logfc.threshold=logFCfilter)
# write.table(markers,file='allcelltype_markers.xls',sep='\t',row.names=F,quote=F)

markers_sub<-subset(markers,((as.numeric(as.vector(markers$avg_log2FC)))>logFCfilter)&(as.numeric(as.vector(markers$p_val_adj))<adjPvalFilter))
write.table(markers_sub,file='allcelltype_markers_logFCf5_adjPval05.xls',sep='\t',row.names=F,quote=F)

# Table S4 The median and interquartile range of cell numbers and proportion per donor for each cell subtypes#########
# 7 Major Cell types(Lineage)
Metadata$celltype <- factor(Metadata$celltype, levels = lineageorder)
table(Metadata$celltype)
# 21 subCell types
table(Metadata$subcelltype)
Metadata$subcelltype <- factor(Metadata$subcelltype,levels = cellorder21)

table_cell <- table(Metadata$ID, Metadata$subcelltype)
head(table_cell)
table_cell <- as.data.frame.array(table_cell)
table_cell <- table_cell[, cellorder21]
head(table_cell)
numsummary=apply(table_cell,2,summary)
numsummary=data.frame(t(numsummary))
sd=data.frame(apply(table_cell,2,sd))
names(sd)="SD"
numsummary=cbind(numsummary,sd)
numsummary$IQR=paste0(round(numsummary$Median,2)," (",round(numsummary$X1st.Qu.,2),",",round(numsummary$X3rd.Qu.,2),")")
numsummary$celltype=rownames(numsummary)
numsummary=numsummary[,c('celltype','IQR')]
names(numsummary)=c('celltype','numIQR')

table_cell <- table(Metadata$ID, Metadata$subcelltype) 
head(table_cell)
datp_cell <- prop.table(table_cell, 1)
datp_cell <- datp_cell[, cellorder21] 
datp_cell <- round(datp_cell * 100, 2)
head(datp_cell)
head(datp_cell)
prosummary=apply(datp_cell,2,summary)
prosummary=data.frame(t(prosummary))
sd=data.frame(apply(datp_cell,2,sd))
names(sd)="SD"
prosummary=cbind(prosummary,sd)
prosummary$IQR=paste0(round(prosummary$Median,2)," (",round(prosummary$X1st.Qu.,2),",",round(prosummary$X3rd.Qu.,2),")")
prosummary$celltype=rownames(prosummary)
prosummary=prosummary[,c('celltype','IQR')]
names(prosummary)=c('celltype','proIQR')

cellnum <- table(Metadata$subcelltype) 
cellnum <- as.data.frame.array(cellnum)
cellnum$celltype=rownames(cellnum)
cellnum=cellnum[,c('celltype','cellnum')]

cellpro <- matrix(table(Metadata$subcelltype))
cellpro <- prop.table(cellpro,2)
cellpro <- as.data.frame.array(cellpro)
rownames(cellpro)=rownames(cellnum)
cellpro$celltype=rownames(cellpro)
names(cellpro)=c('cellpro','celltype')

res <- cellnum %>%
  left_join(cellpro) %>%
  left_join(numsummary) %>%
  left_join(prosummary)
fwrite(res,"subcelltype_summary.xls",sep="\t",quote=F)

# Figure 1E Number of each cell type per donor############################################
table_cell <- data.frame(table(Metadata$subcelltype))
head(table_cell)
names(table_cell) <- c("celltype", "numbersum")
table_cell$numbersumlog <- log10(table_cell$numbersum)
rownames(table_cell) <- table_cell$celltype
table_cell <- table_cell[cellorder21, ]
table_cell$celltype <- factor(table_cell$celltype, cellorder21, ordered = TRUE)
cellorder21_factor <- factor(cellorder21, level = cellorder21, ordered = TRUE)

p <- ggplot(table_cell, aes(x = celltype, y = numbersumlog, fill = cellorder21_factor)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cellorder21_color) +
  labs(x = "", y = "Cells (Log10)") +
  theme(panel.grid = element_blank(), 
    panel.background = element_rect(color = "black", fill = "transparent")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("numbersumloghistplot_subcelltype.png", p, width = 15, height = 6)

pdf(file="numbersumloghistplot_subcelltype.pdf",width = 6, height = 4)
print(p)  
dev.off()

# Figure 1F Proportion of each cell type per donor############################################
table_cell <- table(Metadata$ID, Metadata$subcelltype)
head(table_cell)
datp_cell <- prop.table(table_cell, 1)
datp_cell <- as.data.frame.array(datp_cell)
datp_cell <- datp_cell[, cellorder21]
datp_cell <- round(datp_cell * 100, 2)
head(datp_cell)
head(datp_cell)
datp_cell$sampleid <- rownames(datp_cell)
head(datp_cell) # sum(datp_cell[1,])
datp_cell=data.table(datp_cell)
data <- melt(datp_cell,
  id.vars = "sampleid", 
  measure.vars = cellorder21, 
  variable.name = "celltype", 
  value.name = "pro"
) 
head(data)
class(data$pro)
table(is.na(data$pro))
data$pro=as.numeric(data$pro)
data$celltype <- factor(data$celltype,levels=cellorder21)

pdf(file = "proboxplotbreak_subcelltype.pdf", width = 4 , heigh = 6, onefile=FALSE)
ggplot(data, aes(y = celltype, x = pro, fill = pro)) +
  stat_boxplot(geom = "errorbar", width = 0.5) + 
  geom_boxplot(
    position = position_dodge(0.1),
    size = 0.1,
    width = 0.85, 
    fill = revcellorder21_color,
    outlier.shape = NULL, outlier.colour = "NA" 
  ) +
  theme_bw() + 
  theme(
    # legend.position = "top",
    panel.background = element_blank(), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    axis.text.x = element_text(size = 10, colour = "black", face = "plain"),
    axis.text.y = element_text(size = 10, colour = "black", face = "plain"),
  ) +
  ylab("") + 
  xlab("") +
  
  scale_x_break(c(35,55),space = 0.02) + 
  scale_x_break(c(60,85),space = 0.02) +
  theme_classic()
dev.off()
