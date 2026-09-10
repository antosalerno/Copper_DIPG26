#### Phoshoproteomics analysis of DIPGVI and DIPG7 cells treated with TEPA ####

## Project: KingFisher platform (Thermo Fisher, USA) with MagReSyn® Zr-IMAC HP beads (Resyn939Biosciences, South Africa)
## AIM: Analyse differential site phosphorylation between TEPA and control
# author: Antonietta Salerno
# date: 18/08/2026

BiocManager::install(c("clusterProfiler", "enrichplot", "org.Hs.eg.db", "msigdbr"))

# 1. Load Libraries
library(rlang)
library(stats)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(msigdbr) # Contains the gene sets (Hallmarks, KEGG, etc.)
library(readxl)
library(tidyverse)
library(stringr)
library(reshape2)
library(ggrepel)

setwd("~/OneDrive - UNSW/FIlipDIPG")


# ==============================================================================
# STEP 1: Data preparation and normalisation ####


### A- Normalise total proteome for both cell lines ####
process_total_proteome <- function(file_path, line_name) {
  read_tsv(file_path, na = c("na#", "NA", "NaN", "")) %>%
    # Rimuovi contaminanti e reverse hits
    filter(
      is.na(Reverse) | Reverse != "+",
      is.na(`Potential contaminant`) | `Potential contaminant` != "+",
      is.na(`Only identified by site`) | `Only identified by site` != "+"
    ) %>%
    # Mantieni identificatori e colonne di intensità (es. LFQ intensity o Intensity)
    select(`Protein IDs`, `Majority protein IDs`, starts_with("LFQ intensity ")) %>%
    pivot_longer(
      cols = starts_with("LFQ intensity "),
      names_to = "raw_sample_name",
      values_to = "raw_prot_intensity"
    ) %>%
    mutate(
      cell_line = line_name,
      raw_prot_intensity = na_if(raw_prot_intensity, 0),
      log2_prot_intensity = log2(raw_prot_intensity)
    ) %>%
    # Normalizzazione per mediana del campione
    group_by(raw_sample_name) %>%
    mutate(
      norm_prot_intensity = log2_prot_intensity - median(log2_prot_intensity, na.rm = TRUE)
    ) %>%
    ungroup()
}

prot_dipg6 <- process_total_proteome("proteinGroups_VI.txt", "DIPG6")
prot_dipg7 <- process_total_proteome("proteinGroups_007.txt", "DIPG7")

### B - Prepare phosphoproteome datasets ####
process_phospho <- function(file_path, line_name) {
  read_tsv(file_path, na = c("na#", "NA", "NaN", "")) %>%
    filter(
      is.na(Reverse) | Reverse != "+",
      is.na(`Potential contaminant`) | `Potential contaminant` != "+",
      `Localization prob` >= 0.75
    ) %>%
    select(
      `Protein group IDs`,
      Proteins,
      `Positions within proteins`,
      `Amino acid`,
      starts_with("Intensity ") # O 'LFQ intensity' / 'Reporter intensity' a seconda del setup
    ) %>%
    pivot_longer(
      cols = starts_with("Intensity "),
      names_to = "raw_sample_name",
      values_to = "raw_phospho_intensity"
    ) %>%
    mutate(
      cell_line = line_name,
      raw_phospho_intensity = na_if(raw_phospho_intensity, 0),
      log2_phospho_intensity = log2(raw_phospho_intensity),
      site_id = paste(Proteins, paste0(`Amino acid`, `Positions within proteins`), sep = "_")
    ) %>%
    # Normalizzazione per mediana del campione
    group_by(raw_sample_name) %>%
    mutate(
      norm_phospho_intensity = log2_phospho_intensity - median(log2_phospho_intensity, na.rm = TRUE)
    ) %>%
    ungroup()
}

phospho_dipg6 <- process_phospho("Phospho (STY)Sites_VI.txt", "DIPG6")
phospho_dipg7 <- process_phospho("Phospho (STY)Sites_007.txt", "DIPG7")
# ==============================================================================


# ==============================================================================
#### STEP 2. Merge datasets ####

phospho_combined <- bind_rows(phospho_dipg6, phospho_dipg7)

# 2. (Opzionale) Matrice larga per PCA, heatmap o differential analysis (limma)
phospho_matrix <- phospho_combined %>%
  # Crea un identificatore univoco per ciascun campione/condizione
  unite("sample_id", cell_line, raw_sample_name, sep = "_") %>%
  select(site_id, sample_id, norm_phospho_intensity) %>%
  pivot_wider(
    names_from = sample_id,
    values_from = norm_phospho_intensity
  )
# ==============================================================================

# ==============================================================================
#### STEP 3. Valid values filters ####

# 1. Convert to matrix with sites as rownames
mat <- phospho_matrix %>%
  column_to_rownames("site_id") %>%
  as.matrix()

# 2. Identify columns for each experimental group

cols_d6_ctrl <- grep("DIPG6.*Control|DIPGVI.*Control", colnames(mat), ignore.case = TRUE)
#cols_d6_tepa <- grep("DIPG6.*TEPA|DIPGVI.*TEPA", colnames(mat), ignore.case = TRUE)
cols_d7_ctrl <- grep("DIPG7.*Control|DIPG0007.*Control", colnames(mat), ignore.case = TRUE)
#cols_d7_tepa <- grep("DIPG7.*TEPA|DIPG0007.*TEPA", colnames(mat), ignore.case = TRUE)

# 3. Keep the site if it's quantified in at least 2 replicates in the same condition
min_reps <- 2

keep_sites_group <- (
  rowSums(!is.na(mat[, cols_d6_ctrl, drop = FALSE])) >= min_reps |
    rowSums(!is.na(mat[, cols_d6_tepa, drop = FALSE])) >= min_reps |
    rowSums(!is.na(mat[, cols_d7_ctrl, drop = FALSE])) >= min_reps |
    rowSums(!is.na(mat[, cols_d7_tepa, drop = FALSE])) >= min_reps
)

mat_filtered <- mat[keep_sites_group, ]


# Check how many sites left
nrow(mat_filtered) # From 8812 >> 8643 sites

# ==============================================================================

# ==============================================================================
#### STEP 4. Removing NAs ####

### A - Imputation by Perseus####

# Phosphoproteomics datasets typically contain missing values (NAs) arising both from 
# stochastic instrumental sampling (MCAR: Missing Completely At Random) and low-abundance 
# peptides falling below the mass spectrometer's limit of detection (MNAR: Missing Not At Random).
#
# Imputation is performed to retain low-abundance sites and allow parametric statistical 
# testing across complete matrices.
#
# Perseus-like standard approach: replace missing values by drawing random values 
# from a downshifted normal distribution mimicking low-abundance background noise 
# (parameters: downshift = 1.8 SD, width = 0.3 SD relative to the observed sample distribution).

impute_perseus <- function(data_matrix, width = 0.3, downshift = 1.8) {
  set.seed(42)
  apply(data_matrix, 2, function(x) {
    na_idx <- which(is.na(x))
    if (length(na_idx) == 0) return(x)
    
    mu <- mean(x, na.rm = TRUE)
    sigma <- sd(x, na.rm = TRUE)
    
    imputed_mu <- mu - downshift * sigma
    imputed_sigma <- sigma * width
    
    x[na_idx] <- rnorm(length(na_idx), mean = imputed_mu, sd = imputed_sigma)
    return(x)
  })
}

mat_imputed <- impute_perseus(mat_filtered)

# Save imputed matrix
mat_imputed %>%
  as.data.frame() %>%
  rownames_to_column("site_id") %>%
  write_tsv("phospho_matrix_filtered_imputed.txt")


### B - High variance correction ####

# 1. Calculate the variance row by row
row_vars <- apply(mat_imputed, 1, var, na.rm = TRUE)

# 2. Calculate the quantile setting the NAs removal
variance_cutoff <- quantile(row_vars, probs = 0.30, na.rm = TRUE)

# 3. Filter the matrix excluding potential NAs generates
keep_var <- !is.na(row_vars) & (row_vars > variance_cutoff)
mat_imputed <- mat_imputed[keep_var, ]

# Check dimensions of the new matrix
dim(mat_imputed) # nrows: 6050
# ==============================================================================

# ==============================================================================
#### STEP 5. Map protein accession codes to human readable names ####


# 1. Extract distinct site IDs and isolate the primary clean UniProt ID
site_mapping <- tibble(site_id = rownames(mat_imputed)) %>%
  mutate(
    # Extract the protein part before the underscore (e.g. "P04637" from "P04637_S15")
    raw_protein = sub("_.*", "", site_id),
    # Take only the first ID if semicolon-separated
    uniprot_primary = sub(";.*", "", raw_protein),
    # Strip isoform suffix if present (e.g. "P04637-2" -> "P04637")
    uniprot_clean = sub("-.*", "", uniprot_primary),
    # Extract modification residue and position (e.g. "S15")
    residue_pos = sub(".*_", "", site_id)
  )

# Query org.Hs.eg.db for UniProt IDs
mapped_annotations <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(site_mapping$uniprot_clean),
  columns = c("SYMBOL", "GENENAME", "ENTREZID"),
  keytype = "UNIPROT"
) %>%
  # Deduplicate in case of 1:many mappings by taking the first match
  distinct(UNIPROT, .keep_all = TRUE)

# Join back to create the final annotation table
site_annotation <- site_mapping %>%
  left_join(mapped_annotations, by = c("uniprot_clean" = "UNIPROT")) %>%
  mutate(
    # Fallback to UniProt ID if SYMBOL is missing/NA
    gene_symbol = if_else(is.na(SYMBOL) | SYMBOL == "", uniprot_clean, SYMBOL),
    protein_desc = if_else(is.na(GENENAME), "Unknown", GENENAME),
    # Create clean, publication-ready labels (e.g. TP53_S15)
    label_id = paste0(gene_symbol, "_", residue_pos)
  )
# ==============================================================================


# ==============================================================================
#### STEP 6. Differential Analysis (Limma) DIPG6 vs DIPG7 Cntrl ####
library(limma)
library(tidyverse)
library(ggrepel)

# 1. Prepare metadata
sample_info <- tibble(sample_id = colnames(mat_imputed)) %>%
  mutate(
    cell_line = if_else(grepl("DIPG6|DIPGVI", sample_id), "DIPG6", "DIPG7"),
    treatment = if_else(grepl("TEPA", sample_id, ignore.case = TRUE), "TEPA", "Control"),
    group     = factor(paste(cell_line, treatment, sep = "_"))
  )

# 2. Create a design matrix without intercept
design <- model.matrix(~ 0 + group, data = sample_info)
colnames(design) <- levels(sample_info$group)

# 3. Fit global linear model
fit <- lmFit(mat_imputed, design)

# 4. Create contrast matrix
cont_matrix <- makeContrasts(
  # A. Basal cell line differences
  D6_vs_D7_Ctrl = DIPG6_Control - DIPG7_Control,
  
  # B. Response to individual cell line
  TEPA_in_D6    = DIPG6_TEPA - DIPG6_Control,
  TEPA_in_D7    = DIPG7_TEPA - DIPG7_Control,
  
  # C. Main shared effect
  Main_TEPA     = 0.5 * (DIPG6_TEPA - DIPG6_Control) + 0.5 * (DIPG7_TEPA - DIPG7_Control),
  
  # D. Differential response (interaction)
  Interaction   = (DIPG6_TEPA - DIPG6_Control) - (DIPG7_TEPA - DIPG7_Control),
  
  levels = design
)

# 5. Fit contrasts with Bayesian model
fit_cont <- contrasts.fit(fit, cont_matrix)
fit_eb   <- eBayes(fit_cont)

# 6. Extract annotate and map table of differentially expressed sites
site_annotation_df <- if (!is.data.frame(site_annotation)) {
  as.data.frame(site_annotation) %>% rownames_to_column("site_id")
} else {
  site_annotation
}

extract_dea <- function(fit_obj, coef_name, annot_df, filename) {
  topTable(fit_obj, coef = coef_name, number = Inf) %>%
    rownames_to_column("site_id") %>%
    left_join(annot_df, by = "site_id") %>%
    arrange(P.Value) %>%
    {
      write_tsv(., filename)
      .
    }
}

dea_D6_vs_D7_Ctrl <- extract_dea(fit_eb, "D6_vs_D7_Ctrl", site_annotation_df, "DEA_DIPG6_Ctrl_vs_DIPG7_Ctrl.tsv")
dea_TEPA_in_D6    <- extract_dea(fit_eb, "TEPA_in_D6", site_annotation_df, "DEA_TEPA_in_DIPG6.tsv")
dea_TEPA_in_D7    <- extract_dea(fit_eb, "TEPA_in_D7", site_annotation_df, "DEA_TEPA_in_DIPG7.tsv")
dea_Interaction   <- extract_dea(fit_eb, "Interaction", site_annotation_df, "DEA_Interaction_D6_vs_D7.tsv")
dea_Main_TEPA   <- extract_dea(fit_eb, "Main_TEPA", site_annotation_df, "DEA_Main_TEPA.tsv")

# 7. Plot volcano

plot_volcano <- function(df, title_text, filename, pval_cut = 0.005, logfc_cut = 1.0, top_n_each = 15) {
  
  # A. Fallback / check for label_id column
  if (!"label_id" %in% colnames(df)) {
    df <- df %>%
      mutate(
        clean_residue = toupper(stringr::str_extract(site_id, "[STYsty][0-9]+")),
        label_id = dplyr::case_when(
          !is.na(gene_symbol) & gene_symbol != "" & !is.na(clean_residue) ~ paste0(gene_symbol, "_", clean_residue),
          !is.na(gene_symbol) & gene_symbol != "" ~ as.character(gene_symbol),
          TRUE ~ as.character(label_id)
        )
      )
  }
  
  # B. Assign regulation status
  df_plot <- df %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      regulation = dplyr::case_when(
        logFC >= logfc_cut & P.Value <= pval_cut ~ "Up",
        logFC <= -logfc_cut & P.Value <= pval_cut ~ "Down",
        TRUE ~ "Not Significant"
      )
    )
  
  # C. Explicitly extract top N Up and top N Down separately
  top_up <- df_plot %>%
    dplyr::filter(regulation == "Up" & !is.na(label_id) & label_id != "") %>%
    dplyr::slice_min(order_by = P.Value, n = top_n_each, with_ties = FALSE)
  
  top_down <- df_plot %>%
    dplyr::filter(regulation == "Down" & !is.na(label_id) & label_id != "") %>%
    dplyr::slice_min(order_by = P.Value, n = top_n_each, with_ties = FALSE)
  
  sig_labels <- dplyr::bind_rows(top_up, top_down)
  
  # Diagnostic check in console
  message(paste0("[", title_text, "] Labelled: ", nrow(top_up), " Up, ", nrow(top_down), " Down"))
  
  # D. Generate Volcano Plot
  p <- ggplot(df_plot, aes(x = logFC, y = -log10(P.Value))) +
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
    ggrepel::geom_text_repel(
      data = sig_labels,
      aes(label = label_id),
      size = 2.8,
      fontface = "bold",
      box.padding = 0.35,
      point.padding = 0.25,
      max.overlaps = Inf,
      segment.color = "grey50",
      segment.size = 0.4,
      min.segment.length = 0
    ) +
    labs(
      title = title_text,
      subtitle = paste0("|log2FC| > ", logfc_cut, " & P-value < ", pval_cut, 
                        " (Top ", top_n_each, " Up & ", top_n_each, " Down labelled)"),
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
  
  # 5. Save and return
  ggsave(filename, plot = p, width = 9.5, height = 7.5)
  return(p)
}
# Execute for all contrasts 

# Volcano 1: Basal differences between DIPG6 vs DIPG7
volcano_D6_vs_D7_Ctrl <- plot_volcano(
  df         = dea_D6_vs_D7_Ctrl,
  title_text = "Basal Phosphorylation: DIPG6 Control vs DIPG7 Control",
  filename   = "Volcano_DIPG6_Ctrl_vs_DIPG7_Ctrl_SITES.pdf"
)

# Volcano 2: TEPA effect in DIPG6
volcano_TEPA_in_D6 <- plot_volcano(
  df         = dea_TEPA_in_D6,
  title_text = "Differential Phosphorylation: TEPA vs Control in DIPG6",
  filename   = "Volcano_TEPA_in_DIPG6_SITES.pdf"
)

# Volcano 3: TEPA effect in DIPG7
volcano_TEPA_in_D7 <- plot_volcano(
  df         = dea_TEPA_in_D7,
  title_text = "Differential Phosphorylation: TEPA vs Control in DIPG7",
  filename   = "Volcano_TEPA_in_DIPG7_SITES.pdf"
)

# Volcano 4: Interaction (Difference in TEPA response between D6 and D7)
volcano_Interaction <- plot_volcano(
  df         = dea_Interaction,
  title_text = "Differential Response to TEPA: Interaction (DIPG6 vs DIPG7)",
  filename   = "Volcano_Interaction_D6_vs_D7_SITES.pdf"
)

volcano_Main_TEPA <- plot_volcano(
  df         = dea_Main_TEPA,
  title_text = "Differential Phosphorylation: TEPA vs Control in DIPG6&7",
  filename   = "Volcano_Main_TEPA_D6_D7_SITES.pdf"
)
# ==============================================================================


# ==============================================================================
#### STEP 7. KSEA Kinase activity GSEA - TEPA vs CONTROL####

library(KSEAapp)
library(tidyverse)

# 1. Load intergrated database
data(KSData)

# 2. Wrapper function to calculate KSEA, saving the table and generate barplot 
run_ksea_pipeline <- function(dea_df, comparison_name, p_cutoff = 0.05, top_n = 25) {
  
  # A. Clean and format input for KSApp
  ksea_input <- dea_df %>%
    filter(!is.na(gene_symbol) & gene_symbol != "") %>%
    filter(!is.na(P.Value) & !is.na(logFC) & is.finite(logFC)) %>%
    mutate(
      clean_residue = toupper(str_extract(site_id, "[STYsty][0-9]+"))
    ) %>%
    filter(!is.na(clean_residue)) %>%
    transmute(
      Protein      = as.character(sub("_.*", "", site_id)),
      Gene         = as.character(gene_symbol),
      Peptide      = "UNKNOWN",
      Residue.Both = as.character(clean_residue),
      pval         = as.numeric(P.Value),
      FC           = as.numeric(2^(logFC))
    ) %>%
    distinct(Gene, Residue.Both, .keep_all = TRUE) %>%
    as.data.frame()
  
  # B. Calculate KSEA score
  ksea_res <- KSEA.Scores(
    KSData = KSData,
    PX = ksea_input,
    NetworKIN = FALSE
  )
  
  # C. Clean and save the table
  ksea_clean <- as_tibble(ksea_res) %>%
    mutate(
      p.value     = as.numeric(as.vector(p.value)),
      z.score     = as.numeric(as.vector(z.score)),
      m           = as.numeric(as.vector(m)),
      FDR         = as.numeric(as.vector(FDR)),
      Kinase.Gene = as.character(Kinase.Gene),
      direction   = if_else(z.score > 0, "Up", "Down")
    ) %>%
    arrange(p.value)
  
  write_tsv(ksea_clean, paste0("KSEA_results_", comparison_name, ".tsv"))
  
  # D. Filter for plotting
  plot_ksea <- ksea_clean %>%
    filter(p.value < p_cutoff) %>%
    arrange(desc(abs(z.score))) %>%
    slice_head(n = top_n)
  
  if (nrow(plot_ksea) == 0) {
    message(paste0("No kinases with p < ", p_cutoff, " for: ", comparison_name))
    return(ksea_clean)
  }
  
  # E. Barplot
  p_bar <- ggplot(plot_ksea, aes(x = reorder(Kinase.Gene, z.score), y = z.score, fill = direction, alpha = m)) +
    geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
    scale_fill_manual(
      values = c("Down" = "#2B5C8F", "Up" = "#D73027"),
      labels = c("Inhibited (Down)", "Activated (Up)"),
      name = "Regulation"
    ) +
    scale_alpha_continuous(range = c(0.4, 1.0), name = "Substrates (m)") +
    coord_flip() +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", color = "grey40", linewidth = 0.5) +
    labs(
      title = paste0("Inferred Kinase Activity: ", gsub("_", " ", comparison_name)),
      subtitle = paste0("KSEA based on PhosphoSitePlus substrates (p < ", p_cutoff, ")"),
      x = "Kinase",
      y = "Kinase Activity Score (Z-score)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  ggsave(paste0("KSEA_Barplot_", comparison_name, ".pdf"), plot = p_bar, width = 8, height = 6.5)
  return(ksea_clean)
}

# 3. Execute the function for all the contrasts
ksea_D6_vs_D7   <- run_ksea_pipeline(dea_D6_vs_D7_Ctrl, "DIPG6 baseline (Up) vs DIPG7 baseline (Down)")
ksea_TEPA_D6    <- run_ksea_pipeline(dea_TEPA_in_D6,    "TEPA_in_DIPG6")
ksea_TEPA_D7    <- run_ksea_pipeline(dea_TEPA_in_D7,    "TEPA_in_DIPG7")
ksea_Interact   <- run_ksea_pipeline(dea_Interaction,   "Interaction")
ksea_Main_TEPA   <- run_ksea_pipeline(dea_Main_TEPA,   "Main_TEPA_D6_D7")

# ==============================================================================

# ==============================================================================
#### STEP 8. GSEA - TEPA vs CONTROL####

library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(tidyverse)


# 1. Retrieve and merge Hallmark + Reactome gene sets
h_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::mutate(source = "Hallmark") %>%
  dplyr::select(gs_name, gene_symbol, source)

reactome_sets <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") %>%
  dplyr::mutate(source = "Reactome") %>%
  dplyr::select(gs_name, gene_symbol, source)

combined_gene_sets <- bind_rows(h_sets, reactome_sets)

# TERM2GENE data frame for GSEA (column 1: pathway ID, column 2: gene)
combined_term2gene <- combined_gene_sets %>%
  dplyr::select(gs_name, gene_symbol)

# Lookup table to map pathway ID back to source (Hallmark vs Reactome)
pathway_source_map <- combined_gene_sets %>%
  dplyr::distinct(gs_name, source) %>%
  tibble::deframe()


# 2.  GSEA function for Hallmark + Reactome
run_gsea_hallmark_reactome <- function(dea_df, comparison_name, fdr_cutoff = 0.1, top_n_per_source = 15) {
  
  # A. Clean ranked vector using limma's t-statistic
  ranked_df <- dea_df %>%
    dplyr::filter(!is.na(gene_symbol) & gene_symbol != "" & !is.na(t)) %>%
    dplyr::arrange(desc(abs(t))) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
    dplyr::arrange(desc(t))
  
  gene_list <- setNames(ranked_df$t, ranked_df$gene_symbol)
  
  # B. Execute GSEA
  gsea_res <- GSEA(
    geneList      = gene_list,
    TERM2GENE     = combined_term2gene,
    minGSSize     = 10,
    maxGSSize     = 500,
    pvalueCutoff  = 0.2, # Permissive threshold to keep enough pathways for inspection
    pAdjustMethod = "BH",
    eps           = 1e-10,
    verbose       = FALSE
  )
  
  # C. Clean, annotate database source, and save the table
  res_clean <- as_tibble(gsea_res) %>%
    mutate(
      source = pathway_source_map[ID],
      # Clean prefixes and format underscores
      pathway_name = sub("^(HALLMARK_|REACTOME_)", "", ID),
      pathway_name = gsub("_", " ", pathway_name),
      # Truncate overly long pathway names for plotting readability
      pathway_name = stringr::str_trunc(pathway_name, width = 50),
      direction    = if_else(NES > 0, "Up", "Down")
    )
  
  write_tsv(res_clean, paste0("GSEA_Hallmark_Reactome_results_", comparison_name, ".tsv"))
  
  # D. Select top pathways per source (Hallmark vs Reactome)
  plot_data <- res_clean %>%
    filter(p.adjust < fdr_cutoff) %>%
    group_by(source) %>%
    slice_max(order_by = abs(NES), n = top_n_per_source, with_ties = FALSE) %>%
    ungroup()
  
  # Fallback to nominal p-value if FDR is too stringent
  if (nrow(plot_data) == 0) {
    plot_data <- res_clean %>%
      filter(pvalue < 0.05) %>%
      group_by(source) %>%
      slice_max(order_by = abs(NES), n = top_n_per_source, with_ties = FALSE) %>%
      ungroup()
    subtitle_text <- "Top NES per database (p < 0.05)"
  } else {
    subtitle_text <- paste0("Top NES per database (FDR < ", fdr_cutoff, ")")
  }
  
  if (nrow(plot_data) == 0) {
    message(paste0("No significant pathways for: ", comparison_name))
    return(gsea_res)
  }
  
  # E. Faceted Barplot across Hallmark and Reactome
  p_gsea <- ggplot(plot_data, aes(x = tidytext::reorder_within(pathway_name, NES, source), y = NES, fill = direction, alpha = setSize)) +
    geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
    tidytext::scale_x_reordered() +
    facet_wrap(~ source, scales = "free_y", ncol = 1) +
    scale_fill_manual(
      values = c("Down" = "#2B5C8F", "Up" = "#D73027"),
      name   = "Regulation"
    ) +
    scale_alpha_continuous(range = c(0.45, 1.0), name = "Gene Set Size") +
    coord_flip() +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(
      title    = paste0("GSEA Hallmark & Reactome: ", gsub("_", " ", comparison_name)),
      subtitle = subtitle_text,
      x        = NULL,
      y        = "Normalized Enrichment Score (NES)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle    = element_text(hjust = 0.5, size = 10, color = "grey30"),
      strip.background = element_rect(fill = "grey92", color = "grey40"),
      strip.text       = element_text(face = "bold", size = 11),
      legend.position  = "right",
      panel.grid.minor = element_blank()
    )
  
  ggsave(paste0("GSEA_Hallmark_Reactome_Barplot_", comparison_name, ".pdf"), plot = p_gsea, width = 9.5, height = 9)
  return(gsea_res)
}

# 3. Execute for all contrasts
gsea_D6_vs_D7   <- run_gsea_hallmark(dea_D6_vs_D7_Ctrl, "DIPG6 baseline (Up) vs DIPG7 baseline (Down)")
gsea_TEPA_D6    <- run_gsea_hallmark(dea_TEPA_in_D6,    "TEPA_in_DIPG6")
gsea_TEPA_D7    <- run_gsea_hallmark(dea_TEPA_in_D7,    "TEPA_in_DIPG7")
gsea_Interact   <- run_gsea_hallmark(dea_Interaction,   "Interaction_D6_vs_D7")
gsea_Main_TEPA  <- run_gsea_hallmark(dea_Main_TEPA,   "Main_TEPA_D6_D7")

# ==============================================================================









