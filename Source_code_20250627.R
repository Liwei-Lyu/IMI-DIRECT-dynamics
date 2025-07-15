library(vegan)
library(tidyverse)

meta <- read.delim("metadata.txt", sep="\t", row.names=1)

baseline_counts <- read.delim("baseline_species_counts.txt")
endline_counts <- read.delim("endline_species_counts.txt")

feature_table <- merge(baseline_counts, endline_counts, by="Taxon") 
rownames(feature_table) <- feature_table$Taxon
feature_table$Taxon <- NULL
feature_table <- t(feature_table)

min_reads <- min(rowSums(feature_table))
rarefied_data <- rrarefy(feature_table, sample=min_reads)

taxa_meta <- data.frame(
  feature_ID = paste0("feature_", 1:ncol(rarefied_data)),
  feature_name = colnames(rarefied_data)
)

extract_rank <- function(string, rank) {
  str_extract(string, paste0(rank, "_[^|]*")) %>%
    str_replace(paste0(rank, "_"), "")
}

taxa_meta <- taxa_meta %>%
  mutate(
    kingdom = ifelse(grepl("Bacteria", feature_name), "Bacteria",
                     ifelse(grepl("Viruses", feature_name), "Bacteriophage", "Fungi")),
    phylum = extract_rank(feature_name, "phylum"),
    genus = extract_rank(feature_name, "genus"),
    species = extract_rank(feature_name, "species")
  )

calculate_diversity <- function(count_matrix) {
  data.frame(
    SampleID = rownames(count_matrix),
    Richness = specnumber(count_matrix),
    Shannon = diversity(count_matrix, index="shannon"),
    Simpson = diversity(count_matrix, index="simpson")
  )
}

diversity_df <- calculate_diversity(rarefied_data)

diversity_meta <- merge(diversity_df, meta, by="SampleID")

plot_diversity <- function(data, metric) {
  ggplot(data, aes(x=Timepoint, y=.data[[metric]], color=Group)) +
    geom_boxplot() +
    geom_point(aes(group=SubjectID)) +
    labs(y=metric, title=paste(metric, "Diversity"))
}

plot_diversity(diversity_meta, "Shannon")
plot_diversity(diversity_meta, "Simpson")

calculate_prevalence <- function(count_matrix, timepoint) {
  presence <- count_matrix > 0
  colMeans(presence) * 100
}

baseline_prevalence <- calculate_prevalence(rarefied_data[grep("baseline", rownames(rarefied_data)),])
endline_prevalence <- calculate_prevalence(rarefied_data[grep("endline", rownames(rarefied_data)),])

prevalence_df <- data.frame(
  feature_ID = names(baseline_prevalence),
  Baseline = baseline_prevalence,
  Endline = endline_prevalence,
  Change = endline_prevalence - baseline_prevalence
) %>%
  merge(taxa_meta, by="feature_ID")

calculate_relative_abundance <- function(count_matrix) {
  total_counts <- rowSums(count_matrix)
  sweep(count_matrix, 1, total_counts, FUN="/")
}

relative_abundance <- calculate_relative_abundance(feature_table)

analyze_feature_subsets <- function(rel_abundance_matrix, feature_sets) {
  results <- list()
  
  for(set_name in names(feature_sets)) {
    subset_matrix <- rel_abundance_matrix[, feature_sets[[set_name]]]
    
    results[[set_name]] <- list(
      summary_stats = summary(rowSums(subset_matrix, na.rm=TRUE)),
      taxa_names = taxa_metadata[colnames(subset_matrix), "taxa_name"]
    )
  }
  
  return(results)
}

feature_sets <- list(
  core_union = union_features,
  core_shared = shared_features,
  prevalent = prevalent_features
)

subset_results <- analyze_feature_subsets(relative_abundance, feature_sets)

subset_results$core_union$summary_stats
subset_results$core_union$taxa_names



distance_matrix <- vegdist(relative_abundance_data, method = "bray") %>% as.matrix()

dist_all <- function(distance) {
  diag(distance) <- NA
  sample_ids <- unique(gsub("_(bas|end)$", "", rownames(distance)))
  bas_index <- paste(sample_ids, "_bas", sep = "")
  end_index <- paste(sample_ids, "_end", sep = "")
  intra_distances <- diag(distance[bas_index, end_index])
  inter_distances_bas <- sapply(seq_along(bas_index), function(i) {
    median(distance[bas_index[i], bas_index[-i]], na.rm = TRUE)
  })
  inter_distances_end <- sapply(seq_along(end_index), function(i) {
    median(distance[end_index[i], end_index[-i]], na.rm = TRUE)
  })
  data.frame(SubjectID = sample_ids, intra_distances, inter_distances_bas, inter_distances_end)
}

all_distance <- dist_all(distance_matrix)

intra_group_dissimilarity <- sapply(unique(metadata$SubjectID), function(individual) {
  individual_samples <- which(metadata$SubjectID == individual)
  individual_group <- metadata$Group[individual_samples][1]
  other_samples <- which(metadata$Group == individual_group & !metadata$SubjectID %in% individual)
  inter_dissimilarities <- c(as.matrix(distance_matrix)[individual_samples, other_samples])
  data.frame(SubjectID = individual, intra_group_dissimilarity = median(inter_dissimilarities, na.rm = TRUE))
}, simplify = FALSE) %>% do.call(rbind, .)

all_distance <- merge(all_distance, intra_group_dissimilarity, by = "SubjectID")
all_distance$DMI <- (all_distance$inter_distances_bas + all_distance$inter_distances_end) / 2 - all_distance$intra_distances

merged_data <- merge(all_distance, metadata, by = "SubjectID")

long_data <- pivot_longer(merged_data,
                          cols = c("intra_distances", "inter_distances_bas", "inter_distances_end", "intra_group_dissimilarity"),
                          names_to = "Type", values_to = "Value") %>%
  mutate(Type = factor(Type, levels = c("intra_distances", "inter_distances_bas", "inter_distances_end", "intra_group_dissimilarity")))

p_val_table <- long_data %>%
  wilcox_test(Value ~ Type, paired = TRUE, exact = TRUE) %>%
  mutate(FDR = p.adjust(p, method = "BH"),
         signif_label = case_when(FDR < 0.001 ~ "***", FDR < 0.01 ~ "**", FDR < 0.05 ~ "*", FDR < 0.1 ~ "ˆ", TRUE ~ "ns"))
