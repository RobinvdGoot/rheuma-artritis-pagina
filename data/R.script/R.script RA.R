setwd("C:/Users/robin/OneDrive/Documents/J2P4/Github/")
install.packages("BiocManager")
BiocManager::install('Rsubread')
library(Rsubread)
#indexeren
buildindex(
  basename = 'human',
  reference = 'GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 4000,
  indexSplit = TRUE)
#mapping
align.Con1<- align(index = "human", readfile1 = "Data_RA_raw/SRR4785819_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785819_2_subset40k.fastq", output_file = "Con1.BAM")
align.Con2<- align(index = "human", readfile1 = "Data_RA_raw/SRR4785820_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785820_2_subset40k.fastq" ,output_file = "Con2.BAM")
align.Con3<- align(index = "human", readfile1 = "Data_RA_raw/SRR4785828_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785828_2_subset40k.fastq", output_file = "Con3.BAM")
align.Con4<- align(index = "human", readfile1 = "Data_RA_raw/SRR4785831_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785831_2_subset40k.fastq" ,output_file = "Con4.BAM")
align.RA1<- align(index = "human", readfile1 = "Data_RA_raw/SRR4785979_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785979_2_subset40k.fastq" ,output_file = "RA1.BAM")
align.RA2<- align(index = "human", readfile1 = "Data_RA_raw/SRR4785980_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785980_2_subset40k.fastq", output_file = "RA2.BAM")
align.RA3<- align(index = "human", readfile1 = "Data_RA_raw/SRR4785986_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785986_2_subset40k.fastq" ,output_file = "RA3.BAM")
align.RA4<- align(index = "human", readfile1 = "Data_RA_raw/SRR4785988_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785988_2_subset40k.fastq", output_file = "RA4.BAM")
BiocManager::install('Rsamtools')
library(Rsamtools)
samples <- c('Con1', 'Con2', 'Con3', 'Con4', 'RA1', 'RA2','RA3', 'RA4')
#deel 2 count matrix
count_matrix <- featureCounts(
  files = "Con1.BAM",
  annot.ext = "genomic.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE, 
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE)
# alle samples
allcasussamples<- c("Con1.BAM","Con2.BAM","Con3.BAM","Con4.BAM", "RA1.BAM", "RA2.BAM", "RA3.BAM", "RA4.BAM")
count_matrix <- featureCounts(
  files = allcasussamples,
  annot.ext = "genomic.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE)
str(count_matrix)
RA_casus <- count_matrix$counts
head(RA_casus)
colnames(RA_casus) <- c('Con1', 'Con2', 'Con3', 'Con4', 'RA1', 'RA2','RA3', 'RA4')
head(RA_casus)
write.csv(RA_casus, "count_matrix_RA.txt")
#deel 3 statistiek en analyse
RA_casus <- read.table("count matrix goeie.txt", row.names = 1, header=TRUE)
head(RA_casus)
BiocManager::install("DESeq2")
BiocManager::install("KEGGREST")
BiocManager::install("EnhancedVolcano")
BiocManager::install("pathview")

library(DESeq2)
library(KEGGREST)
library(EnhancedVolcano)
library(pathview)
# treatment data set
treatment <- c("control", "control", "control", "control", "Reuma", "Reuma", "Reuma", "Reuma")
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- c("Con1", "Con2", "Con3", "Con4", "RA1", "RA2", "RA3", "RA4")
head(treatment_table)
dds_casus <- DESeqDataSetFromMatrix(countData = RA_casus,
                              colData = treatment_table,
                              design = ~ treatment)
dds_casus <- DESeq(dds_casus)
resultaten <- results(dds_casus)
write.csv(RA_casus,
          file = ,
          row.names = TRUE)
write.table(resultaten, file = 'cassus_RA123', row.names = TRUE, col.names = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)
hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]
#visualisatie: volcano plot
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj')
dev.copy(png, 'VolcanoplotCasus.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()
resultaten[1] <- NULL
resultaten[2:5] <- NULL
pathview(
  gene.data = resultaten,
  pathway.id = "eco2026",  
  species = "eco",          
  gene.idtype = "KEGG",     
  limit = list(gene = 5)    
)
keggLink("pathway", "ANKRD30BL")
#go analyse
BiocManager::install("goseq")
BiocManager::install("geneLenDataBase")
BiocManager::install("org.Dm.eg.db")
library("goseq")
library("geneLenDataBase")

ALL=rownames(resultaten)
res<-as.data.frame(resultaten)
DEG=res %>% 
  filter(padj<0.05)
DEG=rownames(DEG)
#gene factor maken
gene.vector=as.integer(ALL%in%DEG)
names(gene.vector)=ALL
#lets explore this new vector a bit
head(gene.vector)
tail(gene.vector)
#pwf
pwf=nullp(gene.vector,"hg19","geneSymbol")
  pathview(
    gene.data = resultaten,
    pathway.id = "hsa04010",  
    species = "hsa",          
    gene.idtype = "SYMBOL",     
    limit = list(gene = 5)    )
  resultaten[1] <- NULL
  resultaten[2:5] <- NULL 
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  
  # resultaten dataframe maken
  res_df <- as.data.frame(resultaten)
  
  # gene symbols toevoegen
  res_df$SYMBOL <- rownames(res_df)
  
  # conversie SYMBOL -> ENTREZID
  gene_conversion <- bitr(
    res_df$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  # merge
  res_annotated <- merge(
    res_df,
    gene_conversion,
    by = "SYMBOL"
  )
  
  # significante genen selecteren
  sig_genes <- res_annotated %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1)
  
  # KEGG enrichment
  kegg <- enrichKEGG(
    gene = sig_genes$ENTREZID,
    organism = "hsa"
  )
  
  # controleren of resultaten bestaan
  if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
    dotplot(kegg)
  } else {
    print("Geen KEGG pathways gevonden")
  }
  colnames(res_df)
  