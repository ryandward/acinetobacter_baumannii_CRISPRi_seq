# Load data into melted_results
aba_all <- fread(
	"all_counts_seal.tsv.gz",
	header = FALSE,
	col.names = c("spacer", "count", "condition")
)

curated_names <- fread("curated_names.tsv")

aba_all <- aba_all %>% 
	group_by(condition) %>% 
	mutate(cpm = 10^6 * count / sum(count)) %>% 
	left_join(aba_key) %>% 
	left_join(
		curated_names %>% 
			rename("locus_tag" = "AB19606")
	) %>% 
	left_join(melted_results %>% select(unique_name, Pathway) %>% unique)

aba_multi_normalized <- aba_all %>% 
	select(spacer, count, locus_tag, type) %>% 
	inner_join(aba_design %>% filter(experiment == "tube")) %>% 
	left_join(curated_names %>% rename("locus_tag" = "AB19606")) %>% 
	select(spacer, count, locus_tag, timing, drug, dose, rep, type, unique_name) %>% 
	filter(drug %in% c("None", "Colistin", "Imipenem", "Rifampicin", "Meropenem")) %>%
	group_by(condition) %>% 
	mutate(cpm = 10^6*count/sum(count)) %>% 
	mutate(cpm.rel = cpm/mean(cpm[type == "control"])) %>% 
	ungroup %>% group_by(spacer) %>% 
	mutate(logfc = log2(cpm.rel/mean(cpm.rel[timing == "T0"]))) %>% 
	ungroup %>% group_by(condition) %>% 
	mutate(
		logfc.adj = logfc - mean(logfc[type == "control"]),
		fc.adj = 2^logfc.adj,
		treatment = case_when(drug == "None" ~ "none", TRUE ~ "drug"),
		treatment = factor(treatment, levels = c("none", "drug"))
	)

crossing(
	aba_multi_normalized %>% ungroup %>% select(drug, dose) %>% unique, 
	aba_multi_normalized %>% ungroup %>% filter(timing == "T0") %>% select(-drug, -dose) %>% unique
) %>%
	rbind(aba_multi_normalized) %>%
	filter(unique_name %like% "nuoB" & spacer == "ACAATTTGACGATCTTGCAA") %>%
	ggplot(aes(x = timing, y = fc.adj, group = interaction(unique_name, drug))) +
	geom_point(aes(colour = unique_name)) +
	stat_summary(
		fun.data = "mean_sdl",
		geom = "line",
		lwd = 0.75
	) +
	facet_grid(unique_name ~ drug + dose) +
	scale_colour_manual(
		values = c(
			"control" = "dark grey",
			"nuoB" = "#6A3D9A",
			"lpxC" = "#33A02C",
			"murA" = "#FF7F00",
			"glnS" = "#1F78B4"
		)
	)

