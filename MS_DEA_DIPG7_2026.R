#### Mass Spectrometry - Metabolomics of DIPG007 cells treated with TEPA ####
## Project: Orbitrap Exploris 480 coupled with 3000 nanoLC
## AIM: Data classification and filtering to sort only metabolites in the analysis

# author: Antonietta Salerno
# date: 22/07/2026

BiocManager::install(c("clusterProfiler", "enrichplot", "org.Hs.eg.db", "msigdbr"))


# 1. Load Libraries
library(rlang)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(readxl)
library(tidyverse)
library(stringr)
library(reshape2)
library(ggrepel)

setwd("~/OneDrive - UNSW/FIlipDIPG")


# ==============================================================================
#### STEP 1: Prepare Your Data ####


# 1. Load and Clean Raw Input Data

# Read raw metabolomics data and treat custom placeholder strings as NAs
raw_data <- read_csv(
  "Input_for_MA_DIPG007_p38_24h.csv", 
  na = c("na#", "NA", "NaN"), 
  show_col_types = FALSE
)

# Clean string artifacts (whitespace, quotation marks) and filter out empty records
cleaned_data <- raw_data %>% 
  filter(!is.na(sample)) %>% 
  mutate(
    sample = str_remove_all(sample, "'"),
    sample = str_trim(sample),
    label  = str_remove_all(label, "'"),
    label  = str_trim(label)
  ) %>% 
  filter(sample != "" & label != "")


# 2. Build Numeric Matrix (Metabolites as Rows, Samples as Columns)
data_filtered <- cleaned_data %>% 
  select(-label)

# Reshape into a numeric matrix with metabolites along rows and samples along columns
mat_raw <- data_filtered %>% 
  pivot_longer(cols = -sample, names_to = "Metabolite", values_to = "Intensity") %>% 
  pivot_wider(names_from = sample, values_from = Intensity) %>% 
  column_to_rownames("Metabolite") %>% 
  as.matrix()

# 3. Log2 Transformation & Sample-Wise Median Normalization
min_nonzero <- min(mat_raw[mat_raw > 0], na.rm = TRUE)
mat_log2 <- log2(mat_raw + (min_nonzero / 2))

sample_medians <- apply(mat_log2, 2, median, na.rm = TRUE)
global_median  <- median(sample_medians)
norm_offsets   <- sample_medians - global_median

mat_normalized <- sweep(mat_log2, 2, norm_offsets, FUN = "-")

# 4. Export Normalized Data Table
data_normalised_df <- as_tibble(mat_normalized, rownames = "Metabolite")
write_tsv(data_normalised_df, "metabolites_DIPG007_normalized_log2_tData.tsv")


# 5. Quality Control: Pre- vs. Post-Normalization Diagnostic Boxplots

# Open PDF device for diagnostic visualization
pdf("QC_Metabolomics_Normalization_Boxplots.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(7, 4, 3, 1))

# Boxplot before normalization (raw log2-transformed intensities)
boxplot(
  mat_log2, 
  las = 2, 
  outline = FALSE, 
  main = "Pre-Normalization (log2 Raw)", 
  ylab = expression(Log[2]~Intensity),
  col = "grey85"
)

# Boxplot after median centering (medians aligned across all samples)
boxplot(
  mat_normalized, 
  las = 2, 
  outline = FALSE, 
  main = "Post-Normalization (Median-Centered)", 
  ylab = expression(Normalized~Log[2]~Intensity),
  col = "#2B5C8F"
)

dev.off()

# ==============================================================================

# ==============================================================================
#### STEP 2. Generate Volcano plot for differentiallty expressed metabolites in TEPA ####

dea_data <- read_csv("DEA_FDRMetaboanalyst_040826_Anto_FULL_limma.csv", na = c("na#", "NA", "NaN"))

dea_data <- dea_data %>%
  dplyr::rename(
    metabolite = 1, 
    FC          = 2, 
    log2FC      = 3, 
    p_adjusted  = 4, 
    neg_log10p  = 5  
  )

glimpse(dea_data)

# 1. Target pattern matching
target_patterns <- c(
  "ATP", "Adenosine triphosphate", "Adenosine 5'-triphosphate", "Adenosine diphosphate",
  "GTP", "Guanosine triphosphate", "Guanosine 5'-triphosphate", "Guanosine diphosphate",
  "PEP", "Phosphoenolpyruvate", "Phosphoenolpyruvic acid", "Glutamic acid",
  "GSH", "Glutathione", "L-Glutathione", "Reduced glutathione", "gamma-Glutamylleucine",
  "Ascorbate 6-phosphate", "alpha-Tocopherol","Carnosine", "NADPH"
)

# 2. Define Significance Status 
fc_cutoff <- 0.5
p_cutoff  <- 0.1

# 3. Mutate: Calculate adjusted p-value AND its -log10 transformation
dea_data <- dea_data %>% 
  mutate(
    p.adj       = p_adjusted,
    neg_log10_p = -log10(p.adj),
    
    Significance = case_when(
      log2FC >= fc_cutoff  & p.adj < p_cutoff ~ "Up",
      log2FC <= -fc_cutoff & p.adj < p_cutoff ~ "Down",
      TRUE                                     ~ "Not Significant"
    )
  )

# 4. Calculate summary counts
num_up   <- sum(dea_data$Significance == "Up")
num_down <- sum(dea_data$Significance == "Down")

# 5. Extract Top 5 Up/Down + Target Metabolites
top_up   <- dea_data %>% filter(Significance == "Up") %>% slice_max(order_by = log2FC, n = 5)
top_down <- dea_data %>% filter(Significance == "Down") %>% slice_min(order_by = log2FC, n = 5)

target_hits <- dea_data %>% 
  filter(str_detect(metabolite, regex(paste(target_patterns, collapse = "|"), ignore_case = TRUE)))

top_hits <- bind_rows(top_up, top_down, target_hits) %>% 
  distinct(metabolite, .keep_all = TRUE)

# 6. Volcano Plot

plot_volcano <- function(df, title_text, filename, pval_cut = 0.005, logfc_cut = 1.0) {
  
  df_plot <- df %>%
    mutate(
      regulation = case_when(
        log2FC >= logfc_cut & p.adj <= pval_cut ~ "Up",
        log2FC <= -logfc_cut & p.adj <= pval_cut ~ "Down",
        TRUE ~ "Not Significant"
      )
    )
  
  # Deduplicate significant genes 
  sig_labels <- df_plot %>%
    filter(regulation != "Not Significant" & !is.na(metabolite) & metabolite != "") %>%
    group_by(metabolite) %>%
    slice_min(order_by = p.adj, n = 1) %>%
    ungroup()
  
  # Generate plot
  p <- ggplot(df_plot, aes(x = log2FC, y = -log10(p.adj))) +
    geom_point(
      data = subset(df_plot, regulation == "Not Significant"), 
      color = "grey80", 
      alpha = 0.5, 
      size = 1.8
    ) +
    geom_point(
      data = subset(df_plot, regulation != "Not Significant"), 
      aes(color = regulation), 
      alpha = 0.85, 
      size = 2.4
    ) +
    scale_color_manual(values = c("Down" = "#2B5C8F", "Up" = "#D95F02")) +
    geom_vline(xintercept = c(-logfc_cut, logfc_cut), linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_hline(yintercept = -log10(pval_cut), linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_text_repel(
      data = sig_labels,
      aes(label = metabolite),
      size = 2.8,
      fontface = "bold",
      box.padding = 0.25,
      point.padding = 0.2,
      max.overlaps = Inf,
      segment.color = "grey60",
      segment.size = 0.3,
      min.segment.length = 0
    ) +
    labs(
      title = title_text,
      subtitle = paste0("|log2FC| > ", logfc_cut, " & P-value < ", pval_cut),
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(italic(P)-value)),
      color = "Regulation"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
  
  # Save plot in pdf
  ggsave(filename, plot = p, width = 9.5, height = 7.5)
  return(p)
}

volcano_MS_DS_TEPA <- plot_volcano(
  df         = dea_data,
  title_text = "Differential metabolite expression in DIPG6 TEPA vs Control",
  filename   = "MS_Volcano_DIPG6_TEPA.pdf"
)


# ==============================================================================
#### STEP 3 . Plot selected differentially expressed metabolites with boxplot ####

norm_data <- read_tsv("metabolites_DIPG007_normalized_log2_tData.tsv", show_col_types = FALSE)


# 1. Define selected metabolites
metaSel <- c("Ascorbate 6-phosphate", "alpha-Tocopherol", "NADPH", "Glutamic acid")

# 3. Prepara la tabella con log2FC e FDR da dea_data

dea_labels <- dea_data %>%
  filter(.data[[meta_col_dea]] %in% metaSel) %>%
  mutate(
    metabolite = factor(.data[[meta_col_dea]], levels = metaSel),
    stars = case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**",
      p_adjusted < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    stat_label = paste0(
      "log2FC = ", sprintf("%.2f", log2FC), "\n",
      "FDR = ", formatC(p_adjusted, format = "e", digits = 2), " (", stars, ")"
    )
  )

# 4. Intensity matrix in long format
meta_col_dea <- colnames(dea_data)[1] # first column

plot_data <- norm_data %>%
  filter(.data[[meta_col_norm]] %in% metaSel) %>%
  pivot_longer(
    cols = -all_of(meta_col_norm),
    names_to = "Sample",
    values_to = "Intensity"
  ) %>%
  mutate(
    metabolite = factor(.data[[meta_col_norm]], levels = metaSel),
    # Recognise both Cntl and TEPA from column names 
    Group = case_when(
      grepl("Cntl", Sample, ignore.case = TRUE) ~ "Cntrl",
      grepl("TEPA", Sample, ignore.case = TRUE) ~ "TEPA",
      TRUE ~ NA_character_
    ),
    Group = factor(Group, levels = c("Cntrl", "TEPA")),
    log2_intensity = log2(Intensity)
  )

# 5. Boxplot: 2 boxes (Cntrl vs TEPA) for each of the 4 metabolites 
p_box_metabolites <- ggplot(plot_data, aes(x = Group, y = log2_intensity, fill = Group)) +
  geom_boxplot(
    width = 0.5, 
    outlier.shape = NA, 
    alpha = 0.75, 
    color = "black", 
    linewidth = 0.5
  ) +
  geom_point(
    position = position_jitter(width = 0.15, seed = 42),
    shape = 21,
    size = 2.5,
    fill = "white",
    color = "black",
    stroke = 0.8
  ) +
  # Box label con log2FC e FDR
  geom_text(
    data = dea_labels,
    aes(x = 1.5, y = Inf, label = stat_label),
    inherit.aes = FALSE,
    vjust = 1.25,
    size = 3.2,
    fontface = "bold.italic",
    color = "grey20"
  ) +
  facet_wrap(~ metabolite, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Cntrl" = "#2B5C8F", "TEPA" = "#D73027")) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.35))) +
  labs(
    title = "Selected Antioxidant / Redox Metabolites",
    subtitle = "Normalized Peak Intensity distribution (Cntrl vs TEPA)",
    x = NULL,
    y = expression(Normalized~Peak~Intensity~(Log[2])),
    fill = "Condition"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey95", color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# 6. Visualise and save
print(p_box_metabolites)
ggsave("MS_Boxplot_Selected_Metabolites_Cntrl_vs_TEPA.pdf", plot = p_box_metabolites, width = 7.5, height = 6.5)



# ==============================================================================
#### STEP 4. GSEA - TEPA vs CONTROL####

dea_df <- read_csv("~/OneDrive - UNSW/FIlipDIPG/Figure3_AS/Metabolomics_AS_250826/DEA_FDRMetaboanalyst_040826_Anto_FULL_limma.csv")
colnames(dea_df)[1] <- "metabolite"
colnames(dea_df)[3] <- "log2FC"

library(KEGGREST)
library(fgsea)
library(tidyverse)


# A - Download e Build List of KEGG pathways for mouse ####

# Download all mice pathways
mmu_pathways_list <- keggList("pathway", "mmu")
mmu_pathway_ids   <- names(mmu_pathways_list)

# Extract compounds (C-numbers es. C00051) for each pathway
mouse_metabolite_pathways <- map(mmu_pathway_ids, function(path_id) {
  tryCatch({
    entry <- keggGet(path_id)
    compounds <- entry[[1]]$COMPOUND
    if (!is.null(compounds)) {
      return(names(compounds))
    } else {
      return(NULL)
    }
  }, error = function(e) NULL)
})

# Clean all pathways by removing species suffix
clean_pathway_names <- gsub(" - Mus musculus \\(house mouse\\)", "", unname(mmu_pathways_list))
names(mouse_metabolite_pathways) <- clean_pathway_names
mouse_metabolite_pathways <- compact(mouse_metabolite_pathways)

# Download KEGG database for mapping
kegg_compounds <- keggList("compound")
kegg_compounds_clean <- tolower(gsub("[^a-zA-Z0-9]", "", kegg_compounds))

# ==============================================================================

# B - Map metabolite names in KEGG Compound ID ####


find_kegg_id_fast <- function(names_vec) {
  clean_queries <- tolower(gsub("[^a-zA-Z0-9]", "", names_vec))
  map_chr(clean_queries, function(query) {
    if (query == "" || is.na(query)) return(NA_character_)
    match_idx <- which(grepl(paste0("^", query, "$|", query), kegg_compounds_clean))
    if (length(match_idx) > 0) {
      return(names(kegg_compounds)[match_idx[1]])
    } else {
      return(NA_character_)
    }
  })
}
# ==============================================================================

# C- Run function and produce barplot ####


run_kegg_msea_mouse <- function(dea_df, comparison_name, pathway_list = mouse_metabolite_pathways, 
                                padj_cutoff = 0.1, pval_fallback = 0.05, top_n = 20) {
  
  # A. Identify metabolite column 
  name_col <- if ("metabolite" %in% colnames(dea_df)) {
    "metabolite"
  } else {
    colnames(dea_df)[1]
  }
  
  # B. Map on KEGG IDs
  dea_mapped <- dea_df %>%
    mutate(
      clean_metabolite = str_remove(.data[[name_col]], "\\.\\d+$"),
      clean_metabolite = stringr::str_trim(clean_metabolite)
    ) %>%
    mutate(kegg_id = gsub("cpd:", "", find_kegg_id_fast(clean_metabolite))) %>%
    filter(!is.na(kegg_id) & kegg_id != "")
  
  if (nrow(dea_mapped) == 0) {
    message(paste0("No mapped metabolite for: ", comparison_name))
    return(NULL)
  }
  
  # C. Ranked vector and lookup table
  id_to_name <- setNames(dea_mapped$clean_metabolite, dea_mapped$kegg_id)
  
  dea_mapped <- dea_mapped %>%
    mutate(
      rank_stat = if (all(c("log2FC", "neg_log10p") %in% colnames(dea_mapped))) {
        sign(log2FC) * neg_log10p
      } else if ("log2_fc" %in% colnames(dea_mapped)) {
        as.numeric(log2_fc)
      } else if ("log2FC" %in% colnames(dea_mapped)) {
        as.numeric(log2FC)
      } else if ("t" %in% colnames(dea_mapped)) {
        as.numeric(t)
      } else {
        as.numeric(logFC)
      }
    ) %>%
    filter(!is.na(rank_stat) & is.finite(rank_stat)) %>%
    arrange(desc(abs(rank_stat))) %>%
    distinct(kegg_id, .keep_all = TRUE) %>%
    arrange(desc(rank_stat))
  
  ranked_stats <- setNames(dea_mapped$rank_stat, dea_mapped$kegg_id)
  
  # D. Execute fgsea
  fgsea_res <- fgsea(
    pathways = pathway_list,
    stats    = ranked_stats,
    minSize  = 3,
    maxSize  = 500,
    eps      = 1e-10
  )
  
  # E. Clean and translate compound into leading edges
  fgsea_clean <- fgsea_res %>%
    as_tibble() %>%
    arrange(pval) %>%
    mutate(
      leadingEdge_names = map(leadingEdge, ~ unname(id_to_name[.x])),
      leadingEdge       = map_chr(leadingEdge_names, ~ paste(.x[!is.na(.x)], collapse = ", ")),
      direction         = if_else(NES > 0, "Up", "Down")
    ) %>%
    dplyr::select(pathway, pval, padj, NES, size, leadingEdge, direction)
  
  clean_filename <- gsub("[^A-Za-z0-9_]", "_", comparison_name)
  write_tsv(fgsea_clean, paste0("MSEA_Mouse_KEGG_results_", clean_filename, ".tsv"))
  
  # F. Select for plotting
  plot_data <- fgsea_clean %>%
    filter(padj < padj_cutoff) %>%
    arrange(desc(abs(NES))) %>%
    slice_head(n = top_n)
  
  subtitle_text <- paste0("Mouse KEGG Pathways (FDR < ", padj_cutoff, ")")
  
  if (nrow(plot_data) == 0) {
    plot_data <- fgsea_clean %>%
      filter(pval < pval_fallback) %>%
      arrange(desc(abs(NES))) %>%
      slice_head(n = top_n)
    subtitle_text <- paste0("Mouse KEGG Pathways (p < ", pval_fallback, ")")
  }
  
  if (nrow(plot_data) == 0) {
    message(paste0("No significant pathway for: ", comparison_name))
    return(fgsea_clean)
  }
  
  # G. Barplot
  p_msea <- ggplot(plot_data, aes(x = reorder(pathway, NES), y = NES, fill = direction, alpha = size)) +
    geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
    scale_fill_manual(
      values = c("Down" = "#2B5C8F", "Up" = "#D73027"),
      labels = c("Downregulated", "Upregulated"),
      name   = "Regulation"
    ) +
    scale_alpha_continuous(range = c(0.45, 1.0), name = "Metabolite Count") +
    coord_flip() +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(
      title    = paste0("Mouse KEGG Enrichment: ", gsub("_", " ", comparison_name)),
      subtitle = subtitle_text,
      x        = "Metabolic Pathway",
      y        = "Normalized Enrichment Score (NES)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle   = element_text(hjust = 0.5, size = 10, color = "grey30"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  ggsave(paste0("MS_MSEA_Mouse_KEGG_Barplot_DIPG6_TEPA", clean_filename, ".pdf"), plot = p_msea, width = 8.5, height = 6.5)
  return(fgsea_clean)
}

# Execute function
kegg_mouse_TEPA <- run_kegg_msea_mouse(
  dea_df          = dea_data, 
  comparison_name = "Mouse_TEPA_vs_Control"
)

head(kegg_mouse_TEPA)

# Save results table
readr::write_tsv(kegg_mouse_TEPA, "MS_KEGG_MSEA_Mouse_DIPG6_TEPA_vs_Control_results.tsv")

