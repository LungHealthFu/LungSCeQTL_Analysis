
# Figure 5. Dissection of identifying target susceptibility genes in novel lung cancer susceptibility loci.###############################

# sc-eqtl and NLCLC GWAS colocalisation analysis 
qtl_coloc <- function(gene_id, coloc_file, sample_size) {
  eqtl_coloc <- list(
    snp = coloc_file$SNP,
    position = coloc_file$bp,
    beta = coloc_file$beta_eqtl,
    varbeta = coloc_file$varbeta_eqtl,
    type = "quant",
    N = as.numeric(sample_size),
    MAF = coloc_file$maf
  )
  gwas_coloc <- list(
    snp = coloc_file$SNP,
    position = coloc_file$bp,
    beta = coloc_file$beta_gwas,
    varbeta = coloc_file$varbeta_gwas,
    type = "cc"
  )
  result <- coloc.abf(dataset1 = eqtl_coloc, dataset2 = gwas_coloc)
  res <- t(as.data.frame(result$summary))
  rownames(res) <- gene_id
  res <- as.data.frame(res)
  return(res)
}

for (i in 1:nrow(list)) {
  print(list$celltype[i])
  dat <- parse_fst(sprintf("%s/data/%s.fst",dir,list$celltype[i])) %>%
    as.data.frame() %>%
    filter(gene_id %in% genelist$gene_id) %>%
    arrange(gene_id)
  
  sample_size <- list$N[i]
  
  # NSCLC
  gene_list <- genelist %>% pull(gene_id)
  res_all <- data.frame()
  for (j in 1:length(gene_list)) {
    sub <- dat %>% 
      filter(gene_id == gene_list[j])
    print(nrow(sub))
    if (nrow(sub) > 30) {
      res <- qtl_coloc(gene_list[j], sub, sample_size)
      res_all <- rbind(res_all, res)
    }
  }
  res_all <- res_all %>%
    rownames_to_column(var="gene_id")
  fwrite(res_all, sprintf("%s/result_LC_%s_coloc_result.csv",outdir,list$celltype[i]))
}

##twas
for j in $(seq 1 22);do
     Rscript FUSION.assoc_test.R \
     --sumstats ～/result.sumstats \
     --weights ./WEIGHTS/scTWAS.pos \
     --weights_dir ./WEIGHTS/ \
     --ref_ld_chr ./LDREF/1000G.EUR. \
     --chr $j \
     --out /output.txt
done

# Figure 5A Pie chart presents the details of annotating reported lung cancer GWAS loci by combining our sc-eQTL database with a lung snATAC-seq database by Long, E., et al.lung snATAC-seq database.#######################
R
library(data.table)
library(tidyverse)
library(ggforce)

LCnew1=LCnew %>%
  select(Locus_ID,Cytoband,rsid,colocpeak_gene,coloc_gene) %>%
  distinct() %>%
  mutate(
    colocpeaks= case_when(is.na(colocpeak_gene) ~ 0 ,TRUE ~1),
    coloc= case_when(coloc_gene %in% "" ~ 0 ,TRUE ~1)) %>%
  select(Locus_ID,Cytoband,rsid,colocpeaks,coloc) %>%
  distinct()
LCnew1=LCnew1 %>%
  mutate(gourp= case_when(colocpeaks == 1 ~ 1 ,
                          coloc == 1 ~ 2 ,
                          TRUE ~0))
table(LCnew1$gourp)

df <- table(LCnew1$gourp) 
ratio <- prop.table(df) * 100
label <- c("Level 1", "Level 2")
df <- data.frame(ratio, label)

pdf(file=sprintf("%s/pieplot.pdf",funcdir),width=7,heigh=6)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = c('#D293B8', '#E4E4EF'))+
  geom_arc_bar(data=df,stat = "pie",aes(x0=0,y0=0,r0=1,r=2,amount=ratio,fill=label)
  )
dev.off()


# Figure 3B Heatmap shows the PP4 of co-localization for target susceptibility genes in the lung cancer GWAS loci that colocalized with GWAS associations for NSCLC.#######################
# get gene order 
lctype=c("LC","AD",'SC')
colocre=fread(sprintf("%s/eqtlsiggene_colocre.csv",plotdir))
colocre$lctype_gene=paste(colocre$lctype,colocre$gene_name,sep=":")
egenecoloc_sig=subset(colocre,PPH4>=0.7) 
gene=egenecoloc_sig[,c('lctype_gene','lctype','gene_name','seqnames','start')]
gene <- gene[!duplicated(lctype_gene),] 
dim(gene) 
gene$lctype=factor(gene$lctype,level=c('LC','AD','SC'))
gene=gene[order(gene$lctype,gene$seqnames,gene$start),]
gene <- gene[order(gene$lctype), ]
gene <- gene[!duplicated(gene$gene_name),]
gene_labels <- gene$lctype_gene

egenecoloc=colocre[,c('celltype','lctype_gene','PPH4')]
head(colocre$celltype)
length(unique(egenecoloc$celltype))==18

m_pph4 <- pivot_wider(egenecoloc, names_from = "lctype_gene",values_from = "PPH4") %>% tibble::column_to_rownames("celltype") %>% as.data.frame()
m_pph4[is.na(m_pph4)] = 0
m_pph4 <- m_pph4[cell_labels,]
m_pph4 <- m_pph4[,gene_labels]

df=m_pph4
colnames(df) <- gsub("LC:", "", colnames(df))
colnames(df) <- gsub("AD:", "", colnames(df))
colnames(df) <- gsub("SC:", "", colnames(df))

new_colnames <- colnames(df)
dups <- duplicated(new_colnames)

pdf(sprintf("%s/Heatmap_PPH4_alllctype_colocsig_short_wide.pdf",plotdir),  width=10, height=5)
pheatmap(df,
         cluster_rows=FALSE,
         cluster_cols = FALSE,
         border = "white" ,
         cellwidth = 10, cellheight = 10,
         display_numbers = matrix(ifelse(df >= 0.95, "**",ifelse(df >= 0.7, "*","")), nrow = nrow(df)),
         fontsize_number = 7, 
         fontsize_row = 8, 
         fontsize_col = 7,
         angle_col = 45,  
         show_rownames=TRUE,
         show_colnames=TRUE
)
dev.off()

# Figure 5CD #######################
library(data.table)
library(tidyverse)
library(locuscomparer)

# locuszoom regional plot
p <- locuscompare(in_fn1=gwas_fn, in_fn2=eqtl_fn, 
                    population = "EAS", 
                    title1="GWAS", title2="eQTL", 
                    legend_position='topleft', 
                    marker_col1= marker_col, pval_col1=pval_col, marker_col2=marker_col, pval_col2=pval_col,
                    snp=leadsnp,
                    genome = "hg19")+ ggtitle("title")
pdf(sprintf("%s/plot/%s_%s_%s_coloc.pdf",Locusdir,mygene_name, mycell, mylctype), width = 9, height = 5)
print(p)
dev.off()

# eqtl boxplot
p <- ggplot(res1, aes(x = genotype, y = expresidual)) +
  xlab("") + ylab("") +
  ggtitle(paste0(mycell,"\nP: ", myPvalue, "\nP_FDR: ", myFDRP)) +
  geom_violin(color = mycolor, size = 1, trim = FALSE) +
  geom_boxplot(color = mycolor, fill = "white", size = 1, width = 0.1) +
  geom_line(aes(y = median_expresidual, group = 1), color = "black", size = 1) +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5,family = "ArialMT", color = "black"),
    axis.text = element_text(size = 12,family = "ArialMT", color = "black"),
    axis.line = element_line(size = 0.6, lineend = "square", color = "black"),
    axis.ticks = element_line(size = 0.6),
    panel.background = element_rect(fill = 'white'),
    plot.background = element_rect(fill = "white"))
ggsave(sprintf("%s/%s_boxplot.png", outdir, mycell), p, width = 10, height = 10, units = "cm")

# ATAC Sequencing tracks plot
mergelung_cluster=readRDS(sprintf("%s/lung_integrate_macs2peaks.rds",ncdir))
DefaultAssay(mergelung_cluster) <- "peaks"

colocsigre2=fread(sprintf("%s/alllctype_LCgwas_colocsigres_geneuni_peaks.xls",funcdir),data.table=FALSE)
myrsidlist=unique(colocsigre2$rsid) 
length(myrsidlist)
unique(colocsigre2$gene_name1)

regions=unique(colocpeak1$snp)
ranges.show <- StringToGRanges(regions = regions)
ranges.show$color <- "red"
p <- CoveragePlot(
  object = mergelung_cluster,
  region = myregions1,
  features = mygene,
  assay = 'peaks', 
  expression.assay = "SCT", 
  expression.slot = "data",
  peaks = T,
  group.by="subcelltype",
  links = FALSE,
  extend.upstream = 0,
  extend.downstream = 0,
  region.highlight = ranges.show)& scale_fill_manual(values = colororderNC)
ggsave(sprintf("%s/%s_%s.pdf", atacplotdir, mygene,myrsid), plot = p, width = 7, height = 7)
