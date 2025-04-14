###########################Figure 3. Single-cell eQTLs replicate in Bulk219 eQTL database and other studies###############################

# Figure 3A Concordance of allelic direction of effect between independent single-cell eQTLs and the Bulk219 eQTL dataset#######################
R
rm(list = ls())
library(data.table)
library(tidyverse)
library(readxl)
library(ggsci)
library(knitr)

outdir="bulk219_share"
resdir="plot"
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
eqtl.dir <- 'fasteqtl_allcelltype'
cellcolor_df$file=sprintf("%s/%s/%s/result/%s_eqtl.allpairs.txt.gz", eqtl.dir, cellcolor_df$lineage17, cellcolor_df$cellorder17, cellcolor_df$cellorder17)

qtl = fread(sprintf("%s/alleQTL_match_17.xls",dir))
bulk=fread('bulk219_eqtl.allpairs.txt.gz')
qtl = as.data.frame(qtl)
qtl = qtl %>%
	mutate(
		uni_id = paste0(gene_id, '_', variant_id),
		Z_stat = beta/ beta_se
		) %>%
	rename(
		beta_sc = beta, se_sc = beta_se, p_sc = pvalue, 
		pfdr_sc = pfdr, Z_sc = Z_stat
		) %>%
	select(
		uni_id, celltype, beta_sc, se_sc, p_sc, pfdr_sc, Z_sc
		) 

bulk = bulk %>%
	mutate(
		uni_id = paste0(gene_id, '_', variant_id),
		Z_stat = slope/ slope_se,
		pfdr = p.adjust(pval_nominal, method = "fdr")
		) %>%
	rename(
		beta_bulk = slope, se_bulk = slope_se, p_bulk = pval_nominal, 
		pfdr_bulk = pfdr, Z_bulk = Z_stat
		) %>%
	select(
		uni_id, beta_bulk, se_bulk, p_bulk, pfdr_bulk, Z_bulk
		) 

table(qtl$uni_id %in% bulk$uni_id)
qtl_match = qtl[qtl$uni_id %in% bulk$uni_id,]
bulk_match = bulk[bulk$uni_id %in% qtl$uni_id,]

kable(qtl_match %>% group_by(celltype) %>% reframe(qtl_n = n()))

df = merge(qtl_match, bulk_match, all.x = T)
df$Z_ratio = df$Z_sc / df$Z_bulk
df$celltype = factor(df$celltype, levels = cellorder17)
fwrite(df, file = sprintf("%s/df_sc_bulk.csv",outdir))

df = df %>%
  mutate(Z_coordinate = ifelse(Z_ratio>0, 'Consistent', 'Reverse'))
table(df$Z_coordinate)
11211/nrow(df) #0.7439283

df <- table(df$celltype,df$Z_coordinate)
head(df)
datp_cell <- prop.table(df, 1) 
datp_cell <- as.data.frame.array(datp_cell)
datp_cell <- datp_cell[cellorder17, ] 
datp_cell <- round(datp_cell * 100, 2)
datp_cell
df = df %>%
  mutate(
    Z_coordinate = ifelse(Z_ratio>0, 'Consistent', 'Reverse')
  ) %>%
  freq_table(celltype, Z_coordinate) %>%
  select(
    row_cat, col_cat, n, n_row, n_total, percent_row, se_row, lcl_row, ucl_row
  )
df$col_cat = factor(df$col_cat, levels = c('Reverse','Consistent'))
df$row_cat = factor(df$row_cat, levels = c('MP','CM','Mono','DC',"Mast",'NKT','CD8T','CD4T','NaiveB','Plasma','AT2','AT1','Club','Ciliated','Fib','MyoF','VECs'))
df$row_cat = factor(df$row_cat, levels = rev(levels(df$row_cat)))

p = ggplot(df, aes(x = row_cat, y = percent_row, fill = col_cat)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_y_reverse(position = 'right') +
  scale_fill_manual(values = c("#E5E6E7","#5DADE2"))+
  theme_classic()+
  labs(x = NULL, y = "Percentage of eQTLs", fill = 'Effect direction')+
  theme(
    legend.position = 'left',
    legend.text = element_text(size = 15, family = "Helvetica"),
    legend.title = element_text(size = 15, family = "Helvetica"),
    title = element_text(size = 10, family = "Helvetica"),
    axis.title.x = element_text(size = 15, family = "Helvetica"),
    axis.title.y = element_text(size = 15, family = "Helvetica"),
    axis.text.x = element_text(colour = "black", size = 12, family = "Helvetica"),
    axis.text.y = element_text(colour = "black", size = 12, family = "Helvetica")
  )+  
  scale_x_discrete(position = "top")
ggsave(p, file = sprintf("%s/sc_Bulk_effect_direction.pdf",outdir), width = 6, height = 8.5)
