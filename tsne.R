stats <- fread("/home/ryandward/Downloads/Tu Anh/overall_stats.tsv")


library(Rtsne)


tsne_out <- melted_results %>% 
	filter(base %like% "None" & base %like% "T0" & shift %like% "T1") %>%
	filter(type == "perfect") %>%
	pivot_wider(
		id_cols = "spacer", 
		names_from = "condition", 
		values_from = "LFC.adj") %>% 
	column_to_rownames("spacer") %>%
	data.matrix %>%
	Rtsne(perplexity = 100)

tsne_dt <- melted_results %>% 
	filter(base %like% "None" & base %like% "T0" & shift %like% "T1") %>%
	filter(type == "perfect") %>%
	pivot_wider(
		id_cols = "spacer", 
		names_from = "condition", 
		values_from = "LFC.adj") %>%
	select(spacer) %>% 
	inner_join(
		melted_results %>% 
			select(Pathway, unique_name, spacer) %>% unique) %>% 
	cbind(
		tsne_out$Y %>% 
			as_tibble)

scm = function(palette = cols) {
	scale_color_manual(values = palette, na.value = alpha("light gray", 0.25))
}

tsne_dt %>% inner_join(aba_key) %>%
	mutate(
		Pathway = case_when(
			unique_name %like% "atp" ~ "ATP",
			Pathway == "Other" ~ NA_character_, 
			TRUE ~ Pathway)) %>% 
	ggplot(aes(x = V1, y = V2)) + 
	geom_point(aes(colour = Pathway), alpha = 0.5) +
	doc_theme +
	scm(rainbow(6))
