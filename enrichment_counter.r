# Load several packages from CRAN and Bioconductor
require('pacman')
p_load(
data.table,
scales,
edgeR,
statmod,
poolr,
pheatmap,
svglite,
ggplot2,
ggrepel,
Rtsne,
pracma,
colourpicker,
RColorBrewer,
vegan,
tidyverse,
magrittr,
ggtext,
ggforce
)

# Read in a file containing curated names
curated_names <- fread("curated_names.tsv")

# Read in a BED file
aba_bed <- fread(
  "CP046654.1.bed",
  col.names = c(
    "chromosome",
    "left",
    "right",
    "locus_tag",
    "gene_name",
    "strand",
    "coding",
    "completeness"
  )
)

# Read in a file with counts and conditions
aba <- fread(
  "all_counts_seal.tsv.gz",
  header = FALSE,
  col.names = c("spacer", "count", "condition")
)

# Read in a file with spacer-related information
aba_key <- fread("aba_key.tsv")

# Read in a file with experimental design information
aba_design <- fread("ABA1_experimental_design.tsv", na.strings = "NA")

# Merge aba_bed and aba_key data frames
aba_genome <- aba_bed[aba_key[, .(spacer, type, locus_tag, y_pred, target, offset)], on = .(locus_tag)]

# define the experimental design space to only take into consideration "tubes"
aba_design <- aba_design[experiment == "tube"]

# Replace T with t and put parentheses around reps
# publication_design <- copy(aba_design)[, c("timing", "rep") := .(gsub("T", "t", timing), paste0("(", rep, ")"))]


# keep only the counts that are in the experimental design space
aba <- aba[condition %in% aba_design$condition]

# all conditions are induced except ABA1-1
aba_design[, induced := ifelse(condition == "ABA1-1", FALSE, TRUE)]


# combine drug and dose into one column
aba_design[, drug_dose := paste(drug, dose, sep = "_")]

# convert to factor
aba_design[, drug_dose := factor(drug_dose, levels = unique(drug_dose))]
aba_design[, condition := factor(condition, levels = unique(condition))]
aba_design[, timing := factor(timing, levels = unique(timing))]
aba_design[, induced := factor(induced, levels = unique(induced))]
aba_design[, rep := factor(rep, levels = unique(rep))]
aba_design[, drug := factor(drug, levels = unique(drug))]
aba_design[, dose := factor(dose, levels = unique(dose))]


# Create a design matrix for the groups
aba_permut <- model.matrix(~ 0 + timing : drug, data = aba_design) %>%
  set_rownames(aba_design$condition)

aba_permut <- aba_permut[, colSums(aba_permut != 0) > 0]

# get rid of "timing" and "drug_dose" in names of columns of aba_permut
colnames(aba_permut) <- gsub("timing|drug_dose|drug|dose", "", colnames(aba_permut))
colnames(aba_permut) <- gsub("\\(|\\)", "", colnames(aba_permut))
colnames(aba_permut) <- gsub(":", "_", colnames(aba_permut))


# convert single column into a table
# Convert data to a "wide" format with one column per condition
aba_grid <-
  data.table::dcast(
    aba,
    spacer ~ factor(condition, levels = unique(condition)),
    value.var = "count",
    fill = 0
  )

# Convert the data to a matrix, with spacer names as row names
aba_grid_matrix <- data.matrix(aba_grid[,-c("spacer")]) %>%
  set_rownames(aba_grid$spacer)

# Create a DGEList object
dge <- DGEList(counts = aba_grid_matrix, group = paste(aba_design$induced, aba_design$drug_dose), genes = row.names(aba_grid_matrix))

# Normalize
dge <- calcNormFactors(dge)

# Estimate dispersion
dge <- estimateGLMRobustDisp(dge, aba_permut)

# Fit the glmQLFTest
voom_fit <- voomLmFit(dge, aba_permut)


contrasts <- makeContrasts(
  T1 = T1_None - T0_None,
  T2 = T2_None - T0_None,
  T1_Rifampicin = T1_Rifampicin - T1_None,
  T1_Colistin = T1_Colistin - T1_None,
  T1_Col_Rif = T1_Colistin - T1_Rifampicin,
  T1_Imipenem = T1_Imipenem - T1_None,
  T1_Rif_Imi = T1_Rifampicin - T1_Imipenem,
  levels = aba_permut
)
##########################################################################################

# Create a single data.table with all the results, going through each contrast, one at a time
results <- lapply(colnames(contrasts), function(contrast) {
  fit <- contrasts.fit(voom_fit, contrasts = contrasts[, contrast])

  topTable(eBayes(fit), n = Inf) %>%
    # use_series(table) %>%
    data.table(keep.rownames = "spacer") %>%
    mutate(contrast = contrast)
}) %>%
  rbindlist()

targets <- fread("Ab_library.tsv", na.strings = "None")

results <- results %>% left_join(targets)


#############

all_string <- fread("STRG0060QIE.protein.enrichment.terms.v11.5.txt.gz") %>%
  mutate(locus_tag = str_replace(`#string_protein_id`, ".*\\.", "")) %>%
  unique()


gene_groups <- all_string %>%
      # filter(term %in% (all_string %>% group_by(term) %>% tally() %>% pull(unique(term)))) %>%
      group_by(category, term, description) %>%
      summarise(gene_count = n(), locus_tag = list(sort(unique(locus_tag)))) %>%
      mutate(locus_tag_group = vapply(locus_tag, paste, collapse = ",", FUN.VALUE = character(1)))


term_stats <- gene_groups %>% 
  unnest(locus_tag) %>% 
  inner_join(targets %>% 
  select(locus_tag) %>% unique()) %>% 
  group_by(term, gene_count, description) %>% 
  summarize(genes_targeted = n()) 

complete_terms <- term_stats %>% 
  filter(gene_count == genes_targeted)

# only perform enrichments where all genes are available
# gene_groups <- complete_terms %>% inner_join(gene_groups)

repeated_gene_groups <- gene_groups %>% 
  group_by(locus_tag) %>% 
  mutate(times_listed = n()) %>% 
  arrange(locus_tag) %>% 
  ungroup ()


# pick the best annotation for each locus_tag_group, i.e., highest in term, and the lowest in the category_rank
ranked_annotations <- repeated_gene_groups %>%
      group_by(locus_tag_group, category) %>%
      arrange(versionsort::ver_sort(term)) %>%
      slice(n()) %>%
      ungroup() %>%
      mutate(category_rank = case_when(
        category == "Biological Process (Gene Ontology)" ~ 1,
        category == "Molecular Function (Gene Ontology)" ~ 2,
        category == "Cellular Component (Gene Ontology)" ~ 3,
        category == "Protein Domains and Features (InterPro)" ~ 4,
        category == "Protein Domains (SMART)" ~ 5,
        category == "Protein Domains (Pfam)" ~ 6,
        category == "Annotated Keywords (UniProt)" ~ 7,
        category == "Reactome Pathways" ~ 8,
        category == "Subcellular localization (COMPARTMENTS)" ~ 9,
        category == "Local Network Cluster (STRING)" ~ 10,
        TRUE ~ NA_integer_
      )) %>%
      group_by(locus_tag_group) %>%
      filter(category_rank == min(category_rank))

enrichments <- ranked_annotations %>%
  ungroup() %>%
  distinct(locus_tag_group, .keep_all = TRUE) %>%
  select(-locus_tag_group) %>%
  unnest(locus_tag) %>%
  inner_join(term_stats) 


# Get the unique terms
unique_terms <- unique(enrichments$term)

target_spacers_for_terms <- term_stats %>%
  inner_join(enrichments, relationship = "many-to-many") %>%
  inner_join(targets, relationship = "many-to-many")

####################


#########################################################################################

# Split the spacer column by term
locus_tags_list <- split(target_spacers_for_terms$spacer, target_spacers_for_terms$term)

# Find the indices of each set of locus tags in rownames(dge)
gene_indices <- lapply(locus_tags_list, function(locus_tags) which(rownames(dge) %in% locus_tags))


v <- voomWithQualityWeights(dge, aba_permut, plot = TRUE)


v_targets <- v$E %>%
  data.table(keep.rownames = "spacer") %>%
  select(spacer) %>%
  left_join(
    targets %>%
      filter(locus_tag %in% all_string$locus_tag) %>%
      group_by(spacer) %>%
      filter(is.na(target) | target == "None" | (sp_dir != tar_dir & abs(as.numeric(offset)) == min(abs(as.numeric(offset))) & overlap == max(overlap)))
      %>%
      group_by(target)
  )

v_targets[y_pred == "None", y_pred := NA_integer_]

v_targets$y_pred <- as.numeric(v_targets$y_pred)

v_targets[is.na(target) | target == "None", weight := min(v_targets$y_pred, na.rm = TRUE)]

v_targets[spacer == target, weight := max(v_targets$y_pred, na.rm = TRUE)]

v_targets[mismatches >=1, weight := y_pred]

v_targets$weight <- rescale(as.numeric(v_targets$weight), to = c(1, 1000))

# Perform the competitive gene set test for all gene sets
all_sets <- lapply(colnames(contrasts), function(contrast_name) {
  contrast_column <- contrasts[, contrast_name]
  result <- camera(
    v,
    index = gene_indices, design = aba_permut,
    weights = v_targets$weight,
    inter.gene.cor = 0.05,
    contrast = contrast_column
  ) %>%
    data.table(keep.rownames = "term") %>%
    mutate(term = factor(term, levels = unique_terms), contrast = contrast_name)
  result
}) %>%
  do.call(rbind, .)

  #######################


all_sets <- all_sets %>% 
#   inner_join(enrichments) %>% 
#   inner_join(v_targets) %>%
#   inner_join(term_stats) %>% 
#   group_by(contrast, term, description) %>% 
#   nest(locus_tags = locus_tag) %>% 
#   group_by(locus_tags, contrast) %>% 
#   mutate(missing_genes = gene_count - genes_targeted) %>%
#   arrange(FDR, missing_genes) %>% 
#   slice(1) %>% 
#   ungroup %>% 
  rename(guide_count = NGenes) %>%
#   select(term, guide_count, Direction, PValue, FDR, contrast, description, genes_targeted, gene_count)  %>% unique()  %>% 
  data.table()

create_plot <- function(full_data, targets, enrichments, sets) {
  set_name <- deparse(substitute(sets))

  plot_data <- full_data %>%
      inner_join(aba_design) %>%
      filter(induced == TRUE) %>%
      # filter(
      #   drug_dose %like% "None"| 
      #   drug_dose %like% "Imipenem" |
      #   drug_dose %like% "Colistin_0.44" | 
      #   drug_dose %like% "Rifampicin_0.34") %>%

    rename(sample = condition) %>%
    # select(-type) %>%
    # mutate(imipenem = ifelse(is.na(imipenem), "Stock", imipenem), induced = ifelse(is.na(induced), "Stock", induced)) %>%
    group_by(sample) %>%
    # mutate(imipenem = factor(imipenem, levels = c("Stock", "0", "0.125", "0.25"))) %>%
    # mutate(induced = factor(induced, levels = c("Stock", "FALSE", "TRUE"))) %>%
    mutate(cpm = cpm(count)) %>%
    # filter(type %like% "perfect") %>%
    ungroup() %>%
    inner_join(targets, by = "spacer", relationship = "many-to-many") %>%
    inner_join(enrichments, relationship = "many-to-many") %>%
    inner_join(
      rbind(
        # sets %>% filter(Direction == "Down") %>% arrange(FDR) %>% head(6),
        # sets %>% filter(Direction == "Up") %>% arrange(FDR) %>% head(6)
        sets %>%
          arrange(FDR) %>%
          head(12)
      )
    ) %>%
    filter(FDR <= 0.05) %>%
    group_by(factor(drug_dose), factor(induced), term) %>%
    ungroup() %>%
    arrange(FDR) %>%
    # mutate(facet_title = paste(paste0("[", NGenes, " guides ", toupper(Direction), ": ", signif(FDR, 3), "]"))) %>%
    mutate(description = stringr::str_wrap(description, width = 30)) %>%
    mutate(facet_title = paste0("**", term, "**", " — ", Direction, "<br>", description)) %>%
    # mutate(facet_title = sub("([^ \n]+)", "**\\1**", facet_title)) %>%
    mutate(facet_title = gsub("\n", "<br>", facet_title)) %>%
    mutate(facet_title = paste(facet_title, paste0("**FDR** = ", signif(FDR, 2), ", *n* = ", guide_count), sep = "<br>")) %>%
    mutate(facet_title = factor(facet_title, levels = facet_title %>% unique()))

  ggplot(plot_data, aes(y = cpm, x = factor(timing), group = interaction(factor(drug_dose), factor(timing)))) +
    geom_tile(data = data.frame(timing= "T2"), aes(x = timing, y = 0), width = 1, height = Inf, fill = "grey50", alpha = 0.2, inherit.aes = FALSE) +
    geom_sina(aes(weight = as.numeric(weight), color = factor(drug_dose), size = weight), alpha = 0.5, shape = 16) +
    geom_violin(aes(weight = as.numeric(weight)), alpha = 0.25, draw_quantiles = c(0.25, 0.5, 0.75), position = "dodge") +
    facet_wrap(~facet_title, nrow = 3, scales = "free_y") +
    scale_size(range = c(0.25, 2.5)) +
    # use some nice colors for fill and color
    # scale_fill_manual(values = c("0" = "#1F78B4", "0.125" = "#FF7F00", "0.25" = "#E31A1C")) +
    # scale_color_manual(values = c("0" = "#1F78B4", "0.125" = "#FF7F00", "0.25" = "#E31A1C")) +
    ggtitle(set_name) +
    # bottom label should say induced and uninduced instead of true/false
    labs(x = NULL, y = "Counts per Million", color = "drug_dose (µg/mL)", fill = "drug_dose (µg/mL)", size = "Predicted Weight") +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      breaks = c(10^(0:5)),
      labels = scales::label_number(scale_cut = scales::cut_short_scale())
    ) +
    scale_x_discrete(labels = c("TRUE" = "Induced", "FALSE" = "Uninduced")) +
    theme_minimal() +
    theme(
      strip.text = element_markdown(),
      axis.text.x = element_text(size = rel(1.3), color = "black"),
      axis.title.y = element_text(size = rel(1.3), color = "black")
    )
}

create_plot(aba, v_targets, enrichments, all_sets[contrast == "T1_Col_Rif", ])





plot_results <- function(results, targets, enrichments, sets) {
  set_name <- deparse(substitute(sets))

  # Find the terms with the lowest FDR
  lowest_FDR_terms <- sets %>%
  filter(!contrast %in% c("T1", "T2")) %>%
    arrange(FDR) %>%
    pull(term) %>%
    unique() %>%
    head(5)

  # Filter the rows that have those terms
  filtered_sets <- sets %>%
    filter(term %in% lowest_FDR_terms) %>%
    arrange(FDR) %>% rename(winning_contrast = contrast) %>% unique()

  plot_data <- results %>%
    filter(!contrast %in% c("T1", "T2")) %>%
      inner_join(targets, relationship = "many-to-many") %>%
      inner_join(enrichments, relationship = "many-to-many") %>%
      inner_join(
        rbind(filtered_sets)  
      ) %>%
      # filter(FDR <= 0.05) %>%
      # group_by(factor(drug_dose), factor(induced), term) %>%
      ungroup() %>%
      arrange(FDR) %>%
      mutate(description = stringr::str_wrap(description, width = 30)) %>%
      mutate(facet_title = paste0("**", term, "**", " — ", Direction, "<br>", description)) %>%
      mutate(facet_title = gsub("\n", "<br>", facet_title)) %>%
      mutate(facet_title = paste(facet_title, paste0("**FDR** = ", signif(FDR, 2), ", *n* = ", guide_count), winning_contrast, sep = "<br>")) %>%
      mutate(facet_title = factor(facet_title, levels = facet_title %>% unique()))

  ggplot(plot_data, aes(y = logFC, x = factor(contrast))) +
    # geom_tile(data = data.frame(timing= "T2"), aes(x = timing, y = 0), width = 1, height = Inf, fill = "grey50", alpha = 0.2, inherit.aes = FALSE) +
    geom_sina(aes(weight = as.numeric(-log10(adj.P.Val)), color = factor(contrast), size = -log10(adj.P.Val)), alpha = 0.5, shape = 16) +
    geom_violin(aes(weight = as.numeric(-log10(adj.P.Val))), alpha = 0.35, draw_quantiles = c(0.25, 0.5, 0.75), position = "dodge") +
    facet_wrap(~facet_title, nrow = 3, scales = "free_y") +
    scale_size(range = c(0.25, 2.5)) +
    ggtitle(set_name) +
    labs(x = NULL, y = "Log Fold Change", color = "drug_dose (µg/mL)", fill = "drug_dose (µg/mL)", size = "Predicted Weight") +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    # scale_y_continuous(
    #   trans = scales::pseudo_log_trans(base = 10),
    #   breaks = c(10^(0:5)),
    #   labels = scales::label_number(scale_cut = scales::cut_short_scale())
    # ) +
    scale_x_discrete(labels = c("TRUE" = "Induced", "FALSE" = "Uninduced")) +
    theme_minimal() +
    theme(
      strip.text = element_markdown(),
      axis.text.x = element_text(size = rel(1.3), color = "black"),
      axis.title.y = element_text(size = rel(1.3), color = "black")
    )
}

plot_results(results, v_targets, enrichments, all_sets)






