source("enrichment_application.r")


library(dplyr)
library(ggplot2)
library(scales)

score_set <- function(selected_contrast_pattern, pattern = TRUE) {

  if(pattern) {
    these_contrasts <- all_sets |>
      filter(contrast %like% selected_contrast_pattern) |>
      pull(contrast)
  }

  else(
    these_contrasts <- all_sets |>
      filter(contrast %in% selected_contrast_pattern) |>
      pull(contrast)
  )

  all_sets |>
    filter(contrast %like% selected_contrast_pattern) |> # could relace with c("contrast1", "contrast2")
    select(Direction, term, description, FDR, genes_targeted, gene_count, contrast) |>
    mutate(FDR = signif(FDR, 2)) |>
    mutate(prop_targ = signif(genes_targeted / gene_count, 2)) |>
    mutate(sign = ifelse(Direction == "Up", 1, -1)) |>
    mutate(score = sign * -log10(FDR * prop_targ) * sqrt(genes_targeted)) |>
    arrange(desc(abs(score)))
}

# to score everything you could do
score_set("*") |>
  mutate(score = signif(score, 2)) |>
  dcast(term + description + genes_targeted + gene_count ~ contrast, value.var = "score")
# then write to clipboard and paste into excel with clipr::write_clip()

create_enrichment_plot <- function(selected_term, selected_contrast_patterns) {
  # Filter 'term_stats' using the 'selected_term' argument
  title <- term_stats %>%
    filter(term == selected_term) %>%
    pull(description)
  title <- paste(title, " (", selected_term, ")", sep = "")
  
  # Helper function to check if any of the patterns match the contrast
  matches_any_pattern <- function(contrast, patterns) {
    any(sapply(patterns, function(pattern) grepl(pattern, contrast)))
  }
  
  # Construct the plot with dynamic filtering based on 'selected_term' and 'selected_contrast_patterns'
  enrichment_plot <- contrast_assignments %>%
    inner_join(group_assignments) %>%
    inner_join(
      all_sets %>%
        filter(term == selected_term) %>%
        filter(sapply(contrast, matches_any_pattern, selected_contrast_patterns))
    ) %>%
        inner_join(contrast_assignments) %>%
        mutate(contrast = factor(contrast, levels = unique(contrast))) %>%
    mutate(
      label = paste(
        contrast,
        paste("FDR:", signif(FDR, 2)),
        paste(Direction, paste(paste(genes_targeted, gene_count, sep = "/"), " genes present", sep = "")),
        sep = "\n"
      )
    ) %>%
    mutate(label = factor(label, levels = unique(label))) %>%
    inner_join(enrichments) %>%
    inner_join(targets) %>%
    inner_join(annotated_data) %>%
    inner_join(v_targets) %>%
    arrange(assignment) %>%
    mutate(assignment = case_when(
      assignment == 1 ~ "Treatment",
      assignment == -1 ~ "Control"
    )) %>%
    mutate(group = factor(group, levels = unique(group))) %>%
    ggplot(aes(x = as.character(assignment), y = cpm)) +
    geom_tile(aes(alpha = factor(ifelse(FDR <= 0.05, "highlight", "no_highlight"))), width = Inf, height = Inf, fill = "light grey") +
    geom_sina(aes(weight = weight, size = y_pred), shape = 21, color = "darkgrey", lwd = 0.1) +
    geom_violin(aes(weight = weight), alpha = 0.0, draw_quantiles = c(0.25, 0.5, 0.75), lwd = 0.75) +
    scale_alpha_manual(values = c("highlight" = 0.00, "no_highlight" = 0.025), guide = FALSE) +
    scale_size(range = c(0.1, 3)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = c(10^(0:5)), labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    facet_wrap(~label) +
    labs(x = NULL, y = "Counts per Million") +
    ggtitle(title) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
    scale_size_continuous(range = c(1, 5), limits = c(0, 1)) +
    theme_minimal()

  print(enrichment_plot)
}

# variety of things for the cyclines
create_enrichment_plot("CL:1517", "cycline") # heme copper terminal oxidase
create_enrichment_plot("CL:1490", "cycline")
create_enrichment_plot("GO:0015920", "cycline")
create_enrichment_plot("GO:0045263", "cycline")
create_enrichment_plot("GO:0008137", "cycline")
create_enrichment_plot("CL:1490", "cycline")

create_enrichment_plot("GO:0140101", c("colistin", "polymyxin"))


create_enrichment_plot("GO:0071554", "penem") #direct target
create_enrichment_plot("CL:2527", "chlorhexidine" ) #direct target. minor differences in efficacy
create_enrichment_plot("GO:0003887", "floxacin" ) #direct target, different mechanism of action, different class of drug

create_enrichment_plot("CL:332", c("carvacrol", "thymol")) #NADH
create_enrichment_plot("GO:0042886", c("carvacrol", "thymol")) #potential target

create_enrichment_plot("CL:1490", c("fosfo", "indole"))
create_enrichment_plot("GO:0005216", "trimethoprim")

create_enrichment_plot("CL:535", "DMSO")

create_enrichment_plot("CL:2556", c("EDTA", "TPEN"))
create_enrichment_plot("CL:1265", c("EDTA", "TPEN"))
