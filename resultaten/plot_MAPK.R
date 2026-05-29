GO.wall=goseq(pwf,"hg19","geneSymbol")

#How many enriched GO terms do we have
class(GO.wall)
head(GO.wall)
nrow(GO.wall)
enriched.GO=GO.wall$category[GO.wall$over_represented_pvalue<.05]
#NOTE: They recommend using a more stringent multiple testing corrected p value here

#How many GO terms do we have now?
class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
library(GO.db)
capture.output(for(go in enriched.GO[1:258]) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
, file="SigGo.txt")
top10<- GO.wall %>% arrange(over_represented_pvalue) %>%
  slice(1:10)
# ggplot 1.0
ggplot(top10, aes(x = reorder(category, -over_represented_pvalue),
                  y = -log10(over_represented_pvalue))) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 10 meest verrijkte GO‑categorieën",
    x = "GO‑categorie",
    y = "-log10(p‑waarde)"
  ) +
  theme_minimal()
#ggplot 2.0
ggplot(top10, aes(x = -log10(over_represented_pvalue),
                  y = reorder(category, over_represented_pvalue))) +
  geom_point(size = 4, color = "darkred") +
  labs(
    title = "Top 10 GO‑categorieën (dotplot)",
    x = "-log10(p‑waarde)",
    y = "GO‑categorie"
  ) +
  theme_minimal()
#ggplot 3.0
ggplot(top10,
       aes(x = reorder(category, over_represented_pvalue),
           y = -log10(over_represented_pvalue),
           fill = -log10(over_represented_pvalue))) +
  
  geom_col(width = 0.8) +
  
  coord_flip() +
  
  scale_fill_gradient(
    low = "skyblue",
    high = "darkblue"
  ) +
  
  labs(
    title = "Top 10 meest verrijkte GO-categorieën",
    subtitle = "Over-representation analyse",
    x = "GO-categorie",
    y = expression(-log[10](pwaarde)),
    fill = expression(-log[10](pwaarde))
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )
# KEGG maken
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)


library(org.Hs.eg.db)
library(clusterProfiler)

gene_conversion <- bitr(
  res_df$ENSEMBL,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db )
res_df <- as.data.frame(resultaten)

res_df$SYMBOL <- rownames(res_df)

res_annotated <- merge(
  res_df,
  gene_conversion,
  by = "SYMBOL"
)
gene_list <- res_annotated$log2FoldChange

names(gene_list) <- res_annotated$ENTREZID 
pathview(
  gene.data = resultaten,
  pathway.id = "hsa05323",
  species = "hsa",
  gene.idtype = "SYMBOL",     
)
sig_genes <- res_annotated %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)
kegg <- enrichKEGG(
  gene = sig_genes$SYMBOL,
  organism = "hsa"
)

dotplot(kegg)