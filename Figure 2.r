###########################Figure 2. Mapping eQTL across cell types in the Chinese LungSCeQTL cohort###############################

# Figure 2A The number of eSNP1 and eSNP2 identified per cell type#######################
R
library(data.table)
library(tidyverse)
library(ggupset) 
library(RColorBrewer)

resdir="alleqtl"
dir="clump"
setwd(dir)
lineage=c("Myeloid","Mast","TNKcell","Bcell","Epithelial","Stroma","Endothelial")
lineage17=c(rep("Myeloid",4),'Mast' ,rep("TNKcell",3), rep("Bcell",2), rep("Epithelial",4), rep("Stroma",2), "Endothelial")
lineage_color <- c("#DC050C", "#F5BA3E", "#F6D047", "#F7EF5D", "#CAE0AB", "#1965B0", "#882E72")
lineage17_color = c(rep("#DC050C",4),"#F5BA3E",rep("#F6D047",3),rep("#F7EF5D",2),rep("#CAE0AB",4),rep("#1965B0",2),rep("#882E72"))

cellorder17=c('MP','CM','Mono','DC',"Mast",'NKT','CD8T','CD4T','NaiveB','Plasma',
'AT2','AT1','Club','Ciliated','Fib','MyoF','VECs') 
color17=c("#DC050C","#E33B15","#E96A1F","#EF8829","#F5BA3E",
          "#F6C141","#F6D047","#F6DF4E","#F7EF5D","#FFFF99",
          "#CAE0AB","#A7D295","#82C480","#5BB66B",
          "#1965B0","#A6CEE3","#882E72")
cellcolor_df=data.frame(cbind(lineage17,lineage17_color,cellorder17,color17))


all_esnp=fread(sprintf("%s/alleQTL_match_17.xls",dir))
table(all_esnp$celltype)
nrow(all_esnp)
nrow(unique(all_esnp[,c("celltype","gene_id","variant_id")]))
nrow(unique(all_esnp[,c("gene_id","variant_id")]))
nrow(unique(all_esnp[,c("celltype","gene_id")]))
nrow(unique(all_esnp[,"variant_id"]))
nrow(unique(all_esnp[,"gene_id"]))

paste0(nrow(all_esnp),"(",length(unique(all_esnp$gene_variant)),")")
eSNP1=subset(all_esnp,esnp=="eSNP1")
paste0(nrow(eSNP1),"(",length(unique(eSNP1$gene_variant)),")")
eSNP2=subset(all_esnp,esnp=="eSNP2")
paste0(nrow(eSNP2),"(",length(unique(eSNP2$gene_variant)),")")

cellesnp=table(all_esnp$celltype,all_esnp$esnp)
cellesnp <- cellesnp[cellorder17,]  

data<-as.data.frame(table(all_esnp$celltype,all_esnp$esnp))
names(data)=c('celltype',"eQTL_rank","Freq")
data$celltype = factor(data$celltype,levels = cellorder17 ,ordered = T )
data$eQTL_rank = factor(data$eQTL_rank,levels = c('eSNP2','eSNP1'), label=c("Round2", "Round1"), ordered = T ) 

pdf(file='Stackedbar_eSNP1andeSNP2.pdf',width = 9 ,height = 6, onefile=FALSE)
ggplot(data=data,aes(celltype,Freq,fill=eQTL_rank))+
  geom_bar(stat="identity",position="stack", width=0.7,size=0.25)+ 
  scale_fill_manual(values=c("#A6CEE3","#1965B0"))+
  labs(x = "",y = "Number of independent eQTLs")+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(panel.background=element_rect(), 
        axis.line=element_line(colour="black",size=0.25), 
        axis.title=element_text(size=13,color="black"),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1.1,size=12,color="black"),
        axis.text.y = element_text(size=12,color="black"), 
        legend.position= c(0.9,0.85)
  )
dev.off()

# Figure 2B Correlation between the number of independent eQTLs and eGenes per cell type with the number of cells 
# cellnum
Metadata <- fread("Metadata_sceqtl.xls")
head(Metadata)
Metadata <- Metadata[,-1]
dim(Metadata)
table_cell<- data.frame(table(Metadata$subcelltype)) 
head(table_cell)
names(table_cell)=c('celltype','numbersum')
table_cell$numbersumlog=log10(table_cell$numbersum)
rownames(table_cell)=table_cell$celltype
table_cell=table_cell[cellorder17,]
table_cell$celltype <- factor(table_cell$celltype, cellorder17,ordered=TRUE)

all_esnp=fread(sprintf("%s/alleQTL_match_17.xls",dir))
table(all_esnp$celltype)

# eQTL
freq_eqtl<- data.frame(table(all_esnp$celltype)) 
head(freq_eqtl)
names(freq_eqtl)=c('celltype','freq_eqtl')

# eGene
eSNP1=subset(all_esnp,esnp=="eSNP1")
paste0(nrow(eSNP1),"(",length(unique(eSNP1$gene_variant)),")")
table_egene<- data.frame(table(eSNP1$celltype))
head(table_egene)
names(table_egene)=c('celltype','freq_egene')

dat=merge(table_cell,freq_eqtl,by="celltype",sort=F)
dat=merge(dat,table_egene,by="celltype",sort=F)
head(dat)

dat=merge(dat,cellcolor_df,by.x="celltype",by.y="cellorder17",sort=F)
head(dat)
write.table(dat, sprintf("%s/eQTLandeGenenum_cellnum_plot.xls",dir),sep="\t",quote=F,col.names=T,row.names=F)

# celltype Number and independent eQTL Number 
cor(dat$numbersumlog,dat$freq_eqtl,method="pearson",use="complete.obs")# 0.78199
cor.test(dat$numbersumlog,dat$freq_eqtl,method="pearson",use="complete.obs")

# celltype Number and independent eGene Number 
cor(dat$numbersumlog,dat$freq_egene,method="pearson",use="complete.obs")#  0.7847465
cor.test(dat$numbersumlog,dat$freq_egene,method="pearson",use="complete.obs")

lm_model1 <- lm(freq_eqtl ~ numbersumlog, data = dat)
slope1 <- coef(lm_model1)[2]
intercept1 <- coef(lm_model1)[1]

lm_model2 <- lm(freq_egene ~ numbersumlog, data = dat)
slope2 <- coef(lm_model2)[2]
intercept2 <- coef(lm_model2)[1]

p <- ggplot() +
  theme_bw() +  
  geom_point(data = dat, aes(x = numbersumlog, y = freq_eqtl, color = celltype, shape = "eQTL"), size = 3) +
  geom_point(data = dat, aes(x = numbersumlog, y = freq_egene, color = celltype, shape = "eGene"), size = 3) +
  geom_abline(aes(slope = slope1, intercept = intercept1), linetype = "solid") +  
  geom_abline(aes(slope = slope2, intercept = intercept2), linetype = "dashed") + 
  scale_color_manual(values = dat$color17, name = "Cell type") + 
  scale_shape_manual(values = c(16, 17)) +  
  labs(x = "Number of Cells (Log10)", y = "Number of eQTL and eGene") +
  theme(panel.grid = element_blank(),
        legend.position = "right") + 
  guides(color = guide_legend(position = "bottom"), 
  shape = guide_legend(position = "top"))
pdf(file=sprintf("%s/eQTLandeGenenum_cellnum_plot.pdf",dir), width = 6.5, height = 8)
print(p)
dev.off()

# Figure 2C Pairwise sharing by the magnitude of eQTLs among cell types###############################################
# For each pair of cell types, we considered the eQTLs that were significant (lfsr < 0.05) in at least one of the two cell types and plotted the proportion of these shared in magnitude (i.e., effect estimates with the same direction and within a factor of 2 in size). Brackets around cell type labels highlight cell types that show higher sharing within their respective lineage. Mean of pairwise percentage sharing per lineage is shown in black. 
# Examples include the allelic effects of eSNPs in AT2 cells versus NKT cells AT1 cells .

library(cowplot)
library(ggcorrplot)
library(colorRamp2)
library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)
library(ggpubr)
library(ggthemes)

mydir="/data1/sceqtl/eQTL/share_imput/Mashr" #39服务器
datadir=sprintf("%s/data",mydir)
plotdir=sprintf("%s/plot",mydir)

# magnitude
result=read.table(sprintf("%s/sharing_res.xls",plotdir),sep="\t",h=T,stringsAsFactors = F)
result1=result[,c('mycell.i','mycell.j','sharedpro_magnitude')]
result2=result[,c('mycell.i','mycell.j','sharedpro_magnitude')]
names(result2)=c('mycell.j','mycell.i','sharedpro_magnitude')
result3=data.frame(cellorder17,cellorder17,rep(1,17))
names(result3)=c('mycell.i','mycell.j','sharedpro_magnitude')
plotres=rbind(result1,result2,result3)

plotres_wide=spread(plotres,mycell.i,sharedpro_magnitude)
rownames(plotres_wide)=plotres_wide$mycell.j
plotres_wide=plotres_wide[,-1]
plotres_wide=plotres_wide[cellorder17,]
plotres_wide=plotres_wide[,cellorder17]
write.table(plotres_wide, file=sprintf("%s/sharing_magnitude_plotfile.xls",plotdir),sep="\t",quote=F,col.names=T,row.names=T)

shared_magnitude<-read.table(sprintf("%s/sharing_magnitude_plotfile.xls",plotdir),sep="\t",h=T,stringsAsFactors = F)
min(shared_magnitude) 
max(shared_magnitude) 
max(shared_magnitude[shared_magnitude != 1], na.rm = TRUE) 

cols <- colorRampPalette(rev(c("#D73027","#EF6E48","#FC8D59","#FEE090","#FFFFBF","#E0F3F8","#91BFDB","#4575B4")))(64)
pdf(file=sprintf("%s/corrplot_mash_sharing_magnitude_full.pdf",plotdir),width=12,height=12)
corrplot(corr = as.matrix(shared_magnitude),  is.corr = FALSE,method = "shade",type = "full", win.asp = 1, tl.pos = "lt", col.lim=c(0,1),tl.col="black",col = cols)
dev.off()

# example AT2,AT1 beta
data.strong=readRDS(sprintf("%s/data_strong/%s.%s.data.strong.rds",datadir,"AT1","AT2"))
m2=readRDS(sprintf("%s/res/%s.%s.m.res.rds",datadir,"AT1","AT2"))
dim(m2$result$PosteriorMean) 
dim(data.strong$Shat) 
beta <- data.strong$Bhat
standard.error <- data.strong$Shat
names(m2$result) 
pm.mash        <- m2$result$PosteriorMean
pm.mash.beta   <- pm.mash * standard.error
lfsr           <- m2$result$lfsr
pm.mash.beta=pm.mash.beta[rowSums(lfsr<0.05)>0,]
dim(pm.mash.beta) 
max(pm.mash.beta)
min(pm.mash.beta)
cor(pm.mash.beta)

mydf <- data.frame(pm.mash.beta) %>%
  mutate(ratios = AT2 / AT1) %>%
  mutate(if_shared = case_when(
			ratios >= 0.5 & ratios <= 2 ~ 'shared',
      TRUE ~ 'no_shared'))
mydf$if_shared <- factor(mydf$if_shared, levels = c("shared", "no_shared"))
table(mydf$if_shared)
mydf <- mydf[order(factor(mydf$if_shared, levels = c("no_shared","shared"))), ]

plot_example1 <- ggplot(mydf, aes(AT2, AT1, colour = if_shared)) +
  geom_point(size=0.1) + 
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  theme_minimal(base_size = 14) +  
  theme(legend.position = "none") +
  scale_y_continuous(name = 'Posterior beta of AT1') +
  scale_x_continuous(name = "Posterior beta of AT2") +
  scale_color_grey(start = 0.2, end = 0.8)  
pdf(file=sprintf("%s/Scatterplot_pmbeta_%s_%s.pdf",plotdir,"AT2","AT1"),width=6,height=6)
plot_example1 + theme(panel.border = element_rect(fill=NA,color="gray", size=2, linetype="solid"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

# example AT2,NKT beta
data.strong=readRDS(sprintf("%s/data_strong/%s.%s.data.strong.rds",datadir,"AT2","NKT"))
m2=readRDS(sprintf("%s/res/%s.%s.m.res.rds",datadir,"AT2","NKT"))
dim(m2$result$PosteriorMean) 
dim(data.strong$Shat) 
beta <- data.strong$Bhat
standard.error <- data.strong$Shat
names(m2$result) 
pm.mash        <- m2$result$PosteriorMean
pm.mash.beta   <- pm.mash * standard.error
lfsr           <- m2$result$lfsr
pm.mash.beta=pm.mash.beta[rowSums(lfsr<0.05)>0,]
dim(pm.mash.beta)
max(pm.mash.beta)
min(pm.mash.beta)
cor(pm.mash.beta)

mydf <- data.frame(pm.mash.beta) %>%
  mutate(ratios = AT2 / NKT) %>%
  mutate(if_shared = case_when(
			ratios >= 0.5 & ratios <= 2 ~ 'shared',
      TRUE ~ 'no_shared'))
mydf$if_shared <- factor(mydf$if_shared, levels = c("shared", "no_shared"))
table(mydf$if_shared)

mydf <- mydf[order(factor(mydf$if_shared, levels = c("no_shared","shared"))), ]

plot_example1 <- ggplot(mydf, aes(AT2, NKT, colour = if_shared)) +
  geom_point(size=0.1) + 
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  theme_minimal(base_size = 14) +  
  theme(legend.position = "none") +
  scale_y_continuous(name = 'Posterior beta of NKT') +
  scale_x_continuous(name = "Posterior beta of AT2") +
  scale_color_grey(start = 0.2, end = 0.8) 
pdf(file=sprintf("%s/Scatterplot_pmbeta_%s_%s.pdf",plotdir,"AT2","NKT"),width=6,height=6)
plot_example1 + theme(panel.border = element_rect(fill=NA,color="gray", size=2, linetype="solid"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

# Figure 2D Heatmap for pairwise correlations in allelic effects across cell types#####################################################################################
# The size of each square reflects the number of shared eQTL genes, and the color represents the proportion of non-significant eSNPs after conditional regression, with higher values indicating greater genetic control sharing. 
# Examples include the allelic effects of eSNPs in AT2 cells before and after conditioning on lead eSNPs from NKT cells versus AT1 cells. 
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(reshape2)

workdir = 'ConditionBetweenCelltypes'
outdir = 'plot'
setwd(outdir)

df_correlation2 <- fread("df_correlation2.xls")
df_correlation2$significance_status <- factor(df_correlation2$significance_status, levels = c("significant", "non-significant"))

df_for_corrplot <- df_correlation2 %>% 
  group_by(firstC, secondC, significance_status) %>%  
  summarise (n=n()) %>%
  group_by(firstC,secondC) %>% 
  mutate(freq=n/sum(n)) %>%   
  filter(significance_status=="non-significant")

dat <- table(df_correlation2$firstC,df_correlation2$secondC)
dat1 <-dat %>% 
       reshape2::melt(id.vars="x",variable.name="y")  
df_for_corrplot$cell <- paste(df_for_corrplot$firstC,df_for_corrplot$secondC,sep='_')
dat1$cell <- paste(dat1$Var1,dat1$Var2,sep='_')
nosig_ggplot <- merge(df_for_corrplot,dat1,by='cell',all = TRUE)

head(nosig_ggplot)
head(nosig_ggplot)
nosig_ggplot$firstC=nosig_ggplot$Var1
nosig_ggplot$secondC=nosig_ggplot$Var2
nosig_ggplot <- select(nosig_ggplot,-c(Var1,Var2))
nosig_ggplot$value=ifelse(is.na(nosig_ggplot$n)==TRUE,NA,nosig_ggplot$value)
colnames(nosig_ggplot)[c(5:7)] <-c('Number_raw','Proportion','Total_number_raw')
nosig_ggplot$Number <- nosig_ggplot$Number_raw
summary(nosig_ggplot$Number)
nosig_ggplot$Total_number <- nosig_ggplot$Total_number_raw
nosig_ggplot$Number[nosig_ggplot$Number >= 100]<- 100
nosig_ggplot$Total_number[nosig_ggplot$Total_number >= 100]<- 100
length(unique(nosig_ggplot$firstC))
length(unique(nosig_ggplot$secondC))
nosig_ggplot$firstC<-factor(nosig_ggplot$firstC,levels = cellorder17)
nosig_ggplot$secondC<-factor(nosig_ggplot$secondC,levels = rev(cellorder17)) 

# Number
clrs <- colorRampPalette(rev(c("#ED133C","#ed1941","#FC7E43","#E0F3F8","#91BFDB","#346BB4")))(64)
pdf("share_nonsigprop_num1.pdf",  width = 6.5, height = 6)
ggplot()+
  geom_tile(data=nosig_ggplot,aes(firstC,secondC),fill="white",color="white")+
  geom_point(data=nosig_ggplot,aes(firstC,secondC,size=Number,color=Proportion),shape=15)+
  scale_colour_gradientn(colours = clrs)+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))+
  scale_y_discrete(position = "left")+
  labs(x=NULL,y=NULL)
dev.off()

# example Additional Plots for Figure 2C
# example 1 AT2 vs AT1 
df_example1 <- df_correlation2 %>% filter(firstC=="AT2") %>% filter(secondC=="AT1")
a=as.data.frame(table(df_example1$significance_status))
a[1,2] # sig
a[2,2] # non-sig
sum(a[,2]) #shared eGene-eSNP pairs
a[1,2]/sum(a[,2]) #sigpro
a[2,2]/sum(a[,2]) #nonsigpro

df_example1$significance_status <- factor(df_example1$significance_status, levels = c("significant", "non-significant"))
plot_example1 <- ggplot(df_example1, aes(celltype1_slope, c1_new_beta, colour = significance_status)) +
      geom_point(size=2) + 
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0) +
      theme_pander(base_size=14, nomargin=FALSE) +
      theme(legend.position="none") +
      scale_y_continuous(paste('Regression Coefficient of AT2 given AT1'), limits=c(-0.2,0.2)) +
      scale_x_continuous(paste("Regression Coefficient of AT2"),limits=c(-0.2,0.2)) +
      scale_color_grey(end = 0.8, start = 0.2)
pdf("AT2_con_AT1.pdf",  width=6, height=6)
plot_example1 + theme(panel.border = element_rect(fill=NA,color="gray", size=2, linetype="solid"), 
panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

# example 2 AT2 vs NKT 
df_example2 <- df_correlation2 %>% filter(firstC=="AT2") %>% filter(secondC=="NKT")
df_example2$significance_status <- factor(df_example2$significance_status, levels = c("significant", "non-significant"))
a=as.data.frame(table(df_example2$significance_status))
a[1,2] # sig
a[2,2] # non-sig
sum(a[,2]) #shared eGene-eSNP pairs
a[1,2]/sum(a[,2]) #sigpro
a[2,2]/sum(a[,2]) #nonsigpro
df_example2$significance_status <- factor(df_example2$significance_status, levels = c("significant", "non-significant"))
plot_example2 <- ggplot(df_example2, aes(celltype1_slope, c1_new_beta, colour = significance_status)) +
      geom_point(size=2) + 
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0) +
      theme_pander(base_size=14, nomargin=FALSE) +
      theme(legend.position="none") +
      scale_y_continuous(paste('Regression Coefficient of AT2 given NKT'), limits=c(-0.2,0.2)) +
      scale_x_continuous(paste("Regression Coefficient of AT2"),limits=c(-0.2,0.2)) +
      scale_color_grey(end = 0.8, start = 0.2)
pdf("AT2_con_NKT.pdf",  width=6, height=6)
plot_example1 + theme(panel.border = element_rect(fill=NA,color="gray", size=2, linetype="solid"), 
panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


