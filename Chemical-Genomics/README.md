Perhaps the most important piece of this code is the visualization at the end, which is currently not a function. It is a series of data parsers.

```R


this_term <- "GO:0016676"

title <- term_stats %>%
  filter(term == this_term) %>%
  pull(description)

title = paste(title, " (", this_term, ")", sep = "")

enrichment_plot <- contrast_assignments %>%
  inner_join(
    group_assignments,
    relationship = "many-to-many"
  ) %>%
  inner_join(
    all_sets %>%
      filter(term == this_term) %>%
      arrange(FDR) %>%
      head(24)
  ) %>%
  arrange(FDR) %>%
  filter(FDR <= 0.0005) %>%
  mutate(label = paste(contrast, signif(FDR, 3), paste(Direction, paste("(", paste(genes_targeted,gene_count, sep = "/"), "genes in screen )")), sep = "\n")) %>%
  mutate(label = factor(label, levels = unique(label))) %>%
  inner_join(enrichments) %>%
  inner_join(annotated_data %>% inner_join(v_targets)) %>%
  arrange(assignment) %>%
  mutate(group = factor(group, levels = unique(group))) %>%
  ggplot(aes(x = as.character(assignment), y = cpm)) +
  geom_sina(aes(weight = as.numeric(weight), size = weight, color = group)) +
  geom_violin(aes(weight = as.numeric(weight)), alpha = 0.25, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_size(range = c(0.1, 3)) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(10^(0:5)),
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  facet_wrap(~label) +
  ggtitle(title)

plot(enrichment_plot)



```
