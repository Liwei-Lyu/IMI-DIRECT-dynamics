# Load required packages (consolidated list)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(vegan)
library(phyloseq)
library(Maaslin2)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(ggrepel)
library(openxlsx)
library(lme4)
library(emmeans)
library(performance)

# Data loading ----
input.name <- "phanta"
obser.name <- "_bac"
obser.fullname <- "_bacterial species"
out_dir <- "/Users/nhx573/Library/CloudStorage/OneDrive-UniversityofCopenhagen/work_in_2020-2023/DIRECT_and_stability_studyfrom20230206/results/community_analysis/"
working_dir <- "/Users/nhx573/Library/CloudStorage/OneDrive-UniversityofCopenhagen/work_in_2020-2023/DIRECT_and_stability_studyfrom20230206/processed_data/phanta1.1.0/DIRECT_phanta"
setwd(working_dir)

# Load metadata and marker data
marker_df <- read.delim("biomarkers/marker_df.txt", sep = "\t", row.names = 1)
meta <- read.delim("metadata/meta_20240719.txt", sep = ",", row.names = 1)
meta_adv <- meta[, !colnames(meta) %in% c("Anthropometric", "OGTT", "MGS")]

# Get paired samples
paired_subject.ID <- unique(meta[meta$Timepoint == "Endline" & meta$SequencingID != "unknown", "SubjectID"])
paired_sample.ID <- meta[meta$SubjectID %in% paired_subject.ID, "SampleID"]

# Load additional metadata
meta_all <- readRDS("metadata/meta_all20240803.R")
meta_all_wide <- readRDS("metadata/meta_all_wide20240803.R")
meta_adv_short <- readRDS("results/community_analysis/meta_adv_short.rds")
meta_adv_short_short <- unique(meta_adv_short[, c("SubjectID", "Group")])

# Process lifestyle data
lifestyle <- read.delim2("metadata/Lifestyle_score.txt")
lifestyle_bas <- lifestyle %>% 
  mutate(SubjectID = SampleID,
         SampleID = paste0(SampleID, "_bas")) %>% 
  select("SampleID", "Smoking.Status.BL.n", "Value.vm.hpf.mean.BL.Med.n", "HDI.Med.n", "LifestyleScore", "SubjectID")

lifestyle_end <- lifestyle_bas %>% 
  mutate(SampleID = paste0(SubjectID, "_end"))
lifestyle_all <- rbind(lifestyle_bas, lifestyle_end)
rownames(lifestyle_all) <- lifestyle_all$SampleID

# Merge lifestyle data with metadata
lifestyle_all <- lifestyle_all[meta_adv_short$SampleID,]
meta_adv_short <- merge(meta_adv_short, lifestyle_all[, c("SampleID", "Smoking.Status.BL.n", "Value.vm.hpf.mean.BL.Med.n", "HDI.Med.n", "LifestyleScore")], by = "SampleID")

# Define sample groups
T2D_ids <- meta_adv_short[meta_adv_short$Group == "T2D", "SampleID"]
IGR_ids <- meta_adv_short[meta_adv_short$Group == "IGR", "SampleID"]
NGR_ids <- meta_adv_short[meta_adv_short$Group == "NGR", "SampleID"]
nk_ids <- meta_adv_short[meta_adv_short$Group == "nk", "SampleID"]

T2D_subjectids <- unique(meta_adv_short[meta_adv_short$Group == "T2D", "SubjectID"])
IGR_subjectids <- unique(meta_adv_short[meta_adv_short$Group == "IGR", "SubjectID"])
NGR_subjectids <- unique(meta_adv_short[meta_adv_short$Group == "NGR", "SubjectID"])

# Load feature tables
baseline_feature.table <- read.delim("Phanta_baseline_outputs/counts_species.txt", check.names = FALSE)
endline_feature.table <- read.delim("Phanta_followup/counts_species.txt", check.names = FALSE)
baseline.totalreads <- read.delim("Phanta_baseline_outputs/total_reads.tsv", check.names = FALSE)
endline.totalreads <- read.delim("Phanta_followup/total_reads.tsv", check.names = FALSE)
totalreads <- rbind(baseline.totalreads, endline.totalreads)

# Process feature tables
baseline_feature.table$Taxon_Lineage_with_IDs <- NULL
endline_feature.table$Taxon_Lineage_with_IDs <- NULL

feature.table <- merge(baseline_feature.table, endline_feature.table, by = "Taxon_Lineage_with_Names")
rownames(feature.table) <- feature.table$Taxon_Lineage_with_Names
feature.table$Taxon_Lineage_with_Names <- NULL
feature.table <- as.data.frame(t(feature.table))

# Match sample IDs with metadata
feature.table <- feature.table[intersect(rownames(feature.table), meta$SequencingID),]
lookup <- setNames(meta$SampleID, meta$SequencingID)
rownames(feature.table) <- lookup[rownames(feature.table)]
feature.table <- feature.table[paired_sample.ID,]

# Define color schemes
color.scheme <- c("nk" = "#BF9E91", "IGR" = "#BF826B", "NGR" = "#73564C")
color.scheme.tp <- c("Endline" = "#CB8D8B", "Baseline" = "#E1BD83")
color.scheme.tpgp <- c("nk_Endline" = "black", "IGR_Endline" = "blue", "NGR_Endline" = "#038C3E",
                       "T2D_Endline" = "red", "nk_Baseline" = "grey", "IGR_Baseline" = "dodgerblue",
                       "NGR_Baseline" = "green", "T2D_Baseline" = "orange")