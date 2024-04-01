source("enrichment_application.r")

# do each term separately
# unfortunately, this is not yet a function

this_term <- "CL:1517" 
# this_term <- "CL:1490"
# this_term <- "GO:0015920"
# this_term <- "GO:0045263"


title <- term_stats %>%
  filter(term == this_term) %>%
  pull(description)

title <- paste(title, " (", this_term, ")", sep = "")

enrichment_plot <- contrast_assignments %>%
  inner_join(
    group_assignments,
    relationship = "many-to-many"
  ) %>%
  inner_join(
    all_sets %>%
      filter(term == this_term) %>%
      filter(contrast %like% "cycline") %>%
      inner_join(contrast_assignments) %>%
      mutate(contrast = factor(contrast, levels = unique(contrast)))
  ) %>%
  mutate(
    label = paste(
      contrast,
      paste(
        "FDR:",
        signif(FDR, 2)
      ),
      paste(
        Direction,
        paste(
          paste(
            genes_targeted, gene_count,
            sep = "/"
          ), " genes present",
          sep = ""
        )
      ),
      sep = "\n"
    )
  ) %>%
  mutate(
    label = factor(label, levels = unique(label))
  ) %>%
  inner_join(
    enrichments
  ) %>%
  inner_join(
    targets
  ) %>%
  inner_join(
    annotated_data
  ) %>%
  inner_join(
    v_targets
  ) %>%
  arrange(
    assignment
  ) %>%
  mutate(
    assignment = case_when(
      assignment == 1 ~ "Treatment",
      assignment == -1 ~ "Control"
    )
  ) %>%
  mutate(
    group = factor(group, levels = unique(group))
  ) %>%
  ggplot(
    aes(x = as.character(assignment), y = cpm)
  ) +
  geom_tile(
    aes(alpha = factor(ifelse(FDR <= 0.05, "highlight", "no_highlight"))),
    width = Inf, height = Inf, fill = "light grey"
  ) +
  geom_sina(
    aes(
      weight = weight
    ),
    shape = 21,
    color = "darkgrey",
    lwd = 0.1
    # scale = "width"
  ) +
  geom_violin(
    aes(
      weight = weight
    ),
    alpha = 0.0,
    draw_quantiles = c(0.25, 0.5, 0.75),
    lwd = 0.75
  ) +
  scale_alpha_manual(
    values = c("highlight" = 0.00, "no_highlight" = 0.025), guide = FALSE
  ) +
  scale_size(
    range = c(0.1, 3)
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(10^(0:5)),
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  facet_wrap(~label) +
  labs(
    x = NULL,
    y = "Counts per Million"
  ) +
  ggtitle(title) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  scale_size_continuous(range = c(5, 1), limits = c(0, 1)) +
  theme_minimal()

print(enrichment_plot)
