source("meta_counter.R")

p_load_current_gh("hrbrmstr/hrbrthemes")

operon_conversion <- fread("Operons/AB19606_operon_conversion.tsv")

aba_genome_operons <-
	aba_genome %>% 
	full_join(operon_conversion) 

aba_genome_operons_summary <- aba_genome_operons %>% 
	select(locus_tag, operon) %>% 
	unique %>% 
	group_by(operon) %>% 
	tally %>% 
	inner_join(aba_genome_operons) %>% 
	select(locus_tag, operon, n) %>% 
	inner_join(
		curated_names, 
		by = c("locus_tag" = "AB19606")) %>% 
	select(locus_tag, n, operon, unique_name) %>% 
	unique %>%
	arrange(desc(n), locus_tag)

highly_vulnerable <- -7.851051

old_vulnerability_classifications <-
	median_melted_results %>% 
	filter(base == "None_0_T0") %>%
	filter(type == "perfect") %>%
	mutate(
		Response = case_when(
			medLFC < highly_vulnerable  & FDR < 0.05 ~ "Highly Vulnerable",
			medLFC < -1 & FDR < 0.05 ~ "Vulnerable",
			medLFC > 1 & FDR < 0.05 ~ "Resistant",
			TRUE ~ "No Response")) %>%
	select(AB19606, condition, Response) %>%
	mutate(condition = case_when(
		condition == "None_0_T1 - None_0_T0" ~ "T1",
		condition == "None_0_T2 - None_0_T0" ~ "T2"))	%>%
	filter(!is.na(condition)) 

operon_analysis_vulnerability <- 
	old_vulnerability_classifications %>% 
	inner_join(
		aba_genome_operons_summary, 
		by = c("AB19606" = "locus_tag"))

doc_theme <- theme_ipsum(
	base_family = "Arial", 
	caption_margin = 12,
	axis_title_size = 12,
	axis_col = "black")

#######################
# all operon

melted_results %>% 
	filter(base == "None_0_T0") %>%
	filter(shift %in% c("None_0_T1", "None_0_T2")) %>%
	filter(type == "perfect") %>%
	inner_join(
		operon_conversion %>% 
			group_by(operon) %>% 
			tally %>% 
			arrange(desc(n)) %>% 
			inner_join(operon_conversion), 
				by = c("AB19606" = "locus_tag")) %>% 
	{
		ggplot(data = ., aes(n, LFC, group = n)) +
			scale_x_continuous(
		breaks = seq(
			min(.$n), 
			max(.$n), 
			by = 1)) 
		} +
	doc_theme +
	theme(panel.grid.minor = element_blank()) +
	geom_boxplot() +
	facet_grid(~shift)

##############
# nuo operon

melted_results %>% 
	filter(base == "None_0_T0") %>%
	filter(shift %in% c("None_0_T1", "None_0_T2")) %>%
	# filter(type == "perfect") %>%
	inner_join(
		operon_conversion %>% 
			group_by(operon) %>% 
			tally %>% 
			arrange(desc(n)) %>% 
			inner_join(operon_conversion), 
		by = c("AB19606" = "locus_tag")) %>% 
	filter(operon == "TU286Y-1151") %>% 
	arrange(operon, AB19606, offset) %>% 
	inner_join(
		aba_genome %>% 
			select(
				locus_tag, left, right) %>%
			rename("AB19606" = locus_tag) %>% unique) %>%
	mutate(macro_offset = left + offset) %>%
	ggplot(aes(x = macro_offset, y = LFC.adj, colour = unique_name, size = y_pred, alpha = y_pred)) + 
	geom_point() +
	facet_grid(~shift) +
	doc_theme

########################

melted_results %>% 
	filter(base == "None_0_T0") %>%
	filter(shift %in% c("None_0_T1", "None_0_T2")) %>%
	filter(type == "perfect") %>%
	inner_join(
		operon_conversion %>% 
			group_by(operon) %>% 
			tally %>% 
			arrange(desc(n)) %>% 
			inner_join(operon_conversion), 
		by = c("AB19606" = "locus_tag")) %>% 
	filter(operon == "TU286Y-1151") %>% 
	arrange(operon, AB19606, offset) %>% 
	inner_join(
		aba_genome %>% 
			select(
				locus_tag, left, right) %>%
			rename(
				"AB19606" = locus_tag,
				"left gene coordinate" = left) %>% unique) %>%
	mutate(macro_offset = `left gene coordinate` + offset) %>%
	ggplot(aes(x = `left gene coordinate`, y = LFC.adj, fill = unique_name)) + 
	geom_boxplot() +
	facet_grid(facets = c("type", "shift")) +
	doc_theme

########################

melted_results %>% 
	filter(
		condition == "Colistin_0.44_T1 - None_0_T1" |
		condition == "Colistin_0.44_T2 - None_0_T2") %>%
	filter(type == "perfect") %>%
	inner_join(
		operon_conversion %>% 
			group_by(operon) %>% 
			tally %>% 
			arrange(desc(n)) %>% 
			inner_join(operon_conversion), 
		by = c("AB19606" = "locus_tag")) %>% 
	filter(operon == "TU286Y-1151") %>% 
	arrange(operon, AB19606, offset) %>% 
	inner_join(
		aba_genome %>% 
			select(
				locus_tag, left, right) %>%
			rename(
				"AB19606" = locus_tag,
				"left gene coordinate" = left) %>% unique) %>%
	mutate(macro_offset = `left gene coordinate` + offset) %>%
	ggplot(aes(x = `left gene coordinate`, y = LFC.adj, fill = unique_name)) + 
	geom_boxplot() +
	facet_grid(facets = c("type", "condition")) +
	doc_theme
