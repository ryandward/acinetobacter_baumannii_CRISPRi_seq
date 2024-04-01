# Load several packages from CRAN and Bioconductor
require("pacman")
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

dge <- readRDS(file = "Chemical-Genomics/20240207_data_y_DGElist.rds")

colnames(dge$design) <- gsub("[^[:alnum:]]", "_", colnames(dge$design))
colnames(dge$design) <- gsub("___", " - ", colnames(dge$design))

dge$samples$group <- gsub("[^[:alnum:]]", "_", dge$samples$group)
dge$samples$group <- gsub("___", " - ", dge$samples$group)

contrast_levels <- readRDS(file = "Chemical-Genomics/contrast_levels.rds")
contrast_levels <- gsub("[^[:alnum:]]", "_", contrast_levels)
contrast_levels <- gsub("___", " - ", contrast_levels)

contrast_levels <- c(
  contrast_levels,
  "doripenem_T2_9 - meropenem_T2_13",
  "minocycline_T2_15 - tetracycline_T2_1",
  "minocycline_T2_15 - tigecycline_T2_2",
  "tetracycline_T2_1 - tigecycline_T2_2"
)

contrast_list <- setNames(as.list(contrast_levels), contrast_levels)

# Add the 'levels' argument to the list
contrast_list$levels <- dge$design

# Then use these contrast levels in the makeContrasts function
contrasts <- do.call(makeContrasts, contrast_list)


all_string <- fread("STRG0060QIE.protein.enrichment.terms.v11.5.txt.gz") %>%
  mutate(locus_tag = str_replace(`#string_protein_id`, ".*\\.", "")) %>%
  unique()

targets <- fread("Ab_library.tsv", na.strings = "None")


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
  ungroup()


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


#########################################################################################

most_representative_sets <- enrichments %>%
  inner_join(targets %>% select(spacer, locus_tag)) %>%
  group_by(term) %>%
  arrange(spacer) %>%
  mutate(spacers = paste(spacer, collapse = ",")) %>%
  select(spacers, term, genes_targeted) %>%
  unique() %>%
  group_by(spacers) %>%
  inner_join(enrichments %>% select(term, gene_count)) %>%
  unique() %>%
  mutate(pct_targ = genes_targeted / gene_count) %>%
  ungroup() %>%
  group_by(spacers) %>%
  filter(pct_targ == max(pct_targ)) %>%
  ungroup() %>%
  pull(term)

#########################################################################################

# Split the spacer column by term
sets_to_locus_tags <- split(target_spacers_for_terms$spacer, target_spacers_for_terms$term)

# Find the indices of each set of locus tags in rownames(dge)
sets_to_locus_tags_indices <- lapply(sets_to_locus_tags, function(locus_tags) which(rownames(dge) %in% locus_tags))

v <- voomWithQualityWeights(dge, dge$design, plot = TRUE)

# filter out guides that do not have bona-fide targets, but keep non-targeting guides
# i.e. guides that have no target, or guides that have a target in the same (wrong) direction
# overlap of spacer to gene should be maximal, i.e. 20

v_targets <- v$E %>%
  data.table(keep.rownames = "spacer") %>%
  select(spacer) %>%
  left_join(
    targets %>%
      filter(locus_tag %in% all_string$locus_tag) %>%
      group_by(spacer) %>%
      filter(
        is.na(target) |
          target == "None" |
          (
            sp_dir != tar_dir &
              abs(as.numeric(offset)) == min(abs(as.numeric(offset))) &
              overlap == max(overlap)
          )
      ) %>%
      group_by(target)
  )

# Assign the weight to the guides based on y_pred to be between 1 and 100
v_targets[y_pred == "None", y_pred := NA_integer_]

v_targets$y_pred <- as.numeric(v_targets$y_pred)

v_targets[is.na(target) | target == "None", weight := min(v_targets$y_pred, na.rm = TRUE)]

v_targets[spacer == target, weight := max(v_targets$y_pred, na.rm = TRUE)]

v_targets[mismatches >= 1, weight := y_pred]

v_targets$weight <- rescale(as.numeric(v_targets$weight), to = c(1, 100))

# Perform the competitive gene set test for all gene sets
all_sets <- lapply(colnames(contrasts), function(contrast_name) {
  contrast_column <- contrasts[, contrast_name]
  result <- camera(
    v,
    index = sets_to_locus_tags_indices,
    design = dge$design,
    weights = v_targets$weight,
    inter.gene.cor = 0.05,
    contrast = contrast_column
  ) %>%
    data.table(keep.rownames = "term") %>%
    mutate(term = factor(term, levels = unique_terms), contrast = contrast_name)
  result
}) %>%
  do.call(rbind, .)


all_sets <- all_sets %>%
  inner_join(enrichments) %>%
  inner_join(v_targets) %>%
  inner_join(term_stats) %>%
  group_by(contrast, term, description) %>%
  nest(locus_tags = locus_tag) %>%
  group_by(locus_tags, contrast) %>%
  mutate(missing_genes = gene_count - genes_targeted) %>%
  arrange(FDR, missing_genes) %>%
  ungroup() %>%
  rename(guide_count = NGenes) %>%
  select(term, guide_count, Direction, PValue, FDR, contrast, description, genes_targeted, gene_count) %>%
  unique() %>%
  data.table()


contrast_assignments <- contrasts %>%
  data.table(keep.rownames = "group") %>%
  melt(
    id.vars = "group",
    variable.name = "contrast",
    value.name = "assignment"
  ) %>%
  filter(assignment != 0)

group_assignments <- dge$design %>%
  data.table(keep.rownames = "sample") %>%
  melt(id.vars = "sample", variable.name = "group") %>%
  filter(value != 0) %>%
  select(-value)

original_data <- dge$counts %>%
  data.table(keep.rownames = "spacer") %>%
  melt(
    value.name = "count",
    id.vars = "spacer",
    variable.name = "sample"
  )

annotated_data <- dge$samples %>%
  data.table(keep.rownames = "sample") %>%
  inner_join(original_data) %>%
  group_by(sample) %>%
  mutate(cpm = 1e6 * count / sum(count))