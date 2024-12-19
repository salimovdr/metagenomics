library(dada2)
library(ggplot2)
library(vegan)
library(dplyr)
library(tidyr)


path <- "~/U_Sirius/metagenomics/data/fastqs/trimmed_quality"
fns <- list.files(path, pattern = "\\.fastq\\.gz$", full.names = TRUE)

#Фильтрация
filtpath <- file.path(path, "filtered") 
if (!dir.exists(filtpath)) dir.create(filtpath)

filtered <- file.path(filtpath, basename(fns))
out <- filterAndTrim(fns, filtered, maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)


err <- learnErrors(filtered, multithread=TRUE)
plotErrors(err) 

dada_seq <- dada(filtered, err=err, multithread=TRUE)

seqtab <- makeSequenceTable(dada_seq)
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

taxa <- assignTaxonomy(seqtab.nochim, "~/U_Sirius/metagenomics/data/silva_nr99_v138.2_toSpecies_trainset.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "~/U_Sirius/metagenomics/data/silva_v138.2_assignSpecies.fa.gz")

#Saving
write.table(seqtab.nochim, "~/U_Sirius/metagenomics/data/dada2_out/ASV_abundance_table.txt", sep="\t", quote=FALSE, col.names=NA)
write.table(taxa, "~/U_Sirius/metagenomics/data/dada2_out/Taxonomy_table.txt", sep="\t", quote=FALSE, col.names=NA)

asv_seqs <- getSequences(seqtab.nochim)
asv_headers <- paste0(">ASV", seq_along(asv_seqs))
writeLines(c(rbind(asv_headers, asv_seqs)), "~/U_Sirius/metagenomics/data/dada2_out/ASV_sequences.fasta")
#End saving


#Ординация (можно начинать с этого места, без предобработки)
asv_table <- read.table("~/U_Sirius/metagenomics/data/dada2_out/ASV_abundance_table.txt", header=TRUE, sep="\t", row.names=1)

asv_relative <- sweep(asv_table, 1, rowSums(asv_table), "/")

#PCA
pca <- prcomp(asv_relative, scale.=TRUE)

pca_df <- data.frame(Sample=rownames(asv_relative), PC1=pca$x[,1], PC2=pca$x[,2])

ggplot(pca_df, aes(x=PC1, y=PC2)) +
  geom_point(size=3, alpha=0.7) +
  theme_minimal() +
  labs(title="PCA of ASV Data", x="PC1", y="PC2")


bray_dist <- vegdist(asv_relative, method="bray")

pcoa <- cmdscale(bray_dist, eig=TRUE, k=2)

pcoa_df <- data.frame(
  Sample = rownames(asv_relative),
  PC1 = pcoa$points[,1],
  PC2 = pcoa$points[,2]
)

ggplot(pcoa_df, aes(x = PC1, y = PC2)) +
  geom_point(size=3, alpha=0.7) +
  theme_minimal() +
  labs(title = "PCoA with Bray-Curtis", x = "PC 1", y = "PC 2")

#Преобразование
asv_ids <- paste0("ASV", seq_len(ncol(asv_table)))
colnames(asv_table) <- asv_ids
#asv_table <- t(asv_table) 
asv_table <- as.data.frame(asv_table)
asv_table <- cbind("#NAME" = rownames(asv_table), asv_table)  
rownames(asv_table) <- sapply(strsplit(rownames(asv_table), "_"), `[`, 2)

write.table(asv_table, file = "~/U_Sirius/metagenomics/data/dada2_out/ASV_abundance_table_2.txt", sep = "\t", row.names = FALSE, quote = FALSE)

taxonomy_table <- as.data.frame(taxa)
taxonomy_table <- cbind("#TAXONOMY" = asv_ids, taxonomy_table)  
write.table(taxonomy_table, file = "~/U_Sirius/metagenomics/data/dada2_out/Taxonomy_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)

asv_sequences <- data.frame(ASV_ID = asv_ids, Sequence = colnames(seqtab.nochim))
fasta_file <- file("~/U_Sirius/metagenomics/data/dada2_out/ASV_sequences.fasta", open = "w")
for (i in seq_len(nrow(asv_sequences))) {
  writeLines(paste0(">", asv_sequences$ASV_ID[i]), fasta_file)
  writeLines(asv_sequences$Sequence[i], fasta_file)
}
close(fasta_file)

#Добавление метаинформации
metadata <- read.table("~/U_Sirius/metagenomics/data/sample-metadata.txt", header=FALSE, sep="\t", row.names=1)
colnames(metadata) <- c("subject", "antibiotic", "days")

metadata$Sample <- rownames(metadata)
asv_table$Sample <- rownames(asv_table)

asv_table <- merge(asv_table, metadata, by="Sample")

rownames(asv_table) <- asv_table$Sample
asv_table$Sample <- NULL

asv_table <- asv_table[, -1]

asv_table <- asv_table[, c("subject", "antibiotic", "days", setdiff(colnames(asv_table), c("subject", "antibiotic", "days")))]

asv_data <- asv_table[, !(colnames(asv_table) %in% c("subject", "antibiotic", "days"))]

asv_relative <- sweep(asv_data, 1, rowSums(asv_data), "/")

bray_dist <- vegdist(asv_relative, method="bray")

pcoa <- cmdscale(bray_dist, eig=TRUE, k=2)

pcoa_df <- data.frame(
  Sample = rownames(asv_table),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2],
  subject = asv_table$subject,
  antibiotic = asv_table$antibiotic,
  days = asv_table$days
)

ggplot(pcoa_df, aes(x = PC1, y = PC2, color = subject, shape = antibiotic, size = days)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "PCA",
    x = "PС 1",
    y = "PС 2",
    color = "Subject",
    shape = "Antibiotic Usage",
    size = "Days Since Start"
  )

# Кластеризация
set.seed(218)  
kmeans_result <- kmeans(pcoa_df[, c("PC1", "PC2")], centers=3)

pcoa_df$cluster <- kmeans_result$cluster

head(pcoa_df)

ggplot(pcoa_df, aes(x = PC1, y = PC2, color = factor(cluster), shape = antibiotic, size = days)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "PCA",
    x = "PC 1",
    y = "PC 2",
    color = "Cluster",
    shape = "Antibiotic Usage",
    size = "Days Since Start"
  )

# Исследвоание taxa по clusters
colnames(taxonomy_table)[1] <- "taxonomy"

taxonomy_table$Full_Taxonomy <- apply(taxonomy_table[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")], 1, function(x) {
  paste(na.omit(x), collapse = "; ")
})

asv_table$Sample <- rownames(asv_table)

asv_long <- asv_table %>%
  pivot_longer(
    cols = -c(cluster, Sample, subject, antibiotic, days),  
    names_to = "ASV",                                      
    values_to = "Abundance"                                
  )

asv_with_taxonomy <- merge(asv_long, taxonomy_table[, c("ASV", "Full_Taxonomy")], by = "ASV")

# Группировка по кластеру и Full_Taxonomy
diversity_by_cluster <- asv_with_taxonomy %>%
  group_by(cluster, Full_Taxonomy) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(cluster, desc(TotalAbundance))  

# Выбираем топ-10 таксонов для каждого кластера
top_taxa_by_cluster <- diversity_by_cluster %>%
  group_by(cluster) %>%
  slice_max(TotalAbundance, n = 10) %>%
  ungroup()


