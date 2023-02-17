source("meta_counter.R")

p_load_current_gh("hrbrmstr/hrbrthemes")
p_load(ggallin, gtools)

operon_conversion <- fread("Operons/AB19606_operon_conversion.tsv")

conflict_prefer("extract", "magrittr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

operon_details <- 
	curated_names %>% 
	slice(mixedorder(AB19606)) %>%
	inner_join(
		operon_conversion, by = c("AB19606" = "locus_tag")) %>% 
	group_by(operon) %>% 
	arrange(AB19606) %>%
	summarise(
		operon_genes = paste(unique_name, collapse = ", "),
		operon_AB19606 = paste(AB19606, collapse = ", "),
		operon_AB030 = paste(AB030, collapse = ", ")) 

aba_genome_operons <-
	aba_genome %>% 
	full_join(operon_conversion) 

aba_genome_operons_summary <- aba_genome_operons %>% 
	select(locus_tag, operon) %>% 
	unique %>% 
	group_by(operon) %>% 
	tally(name = "total_size") %>% 
	inner_join(aba_genome_operons) %>% 
	select(locus_tag, operon, total_size) %>% 
	inner_join(
		curated_names, 
		by = c("locus_tag" = "AB19606")) %>% 
	select(locus_tag, total_size, operon, unique_name) %>% 
	unique %>%
	slice(mixedorder(locus_tag))  

operon_pathways <- melted_results %>% 
	filter(type == "perfect") %>%
	select(AB19606, Pathway) %>% 
	rename("locus_tag" = AB19606) %>% 
	unique %>% 
	inner_join(aba_genome_operons) %>% 
	unique %>% 
	select(operon, Pathway, locus_tag) %>% 
	unique %>%
	slice(mixedorder(locus_tag)) %>% 
	group_by(operon) %>% 
	summarise(essential_size = n(), Pathways = paste(unique(Pathway), collapse = ", "))

operon_details <- operon_details %>% inner_join(operon_pathways) %>% 
	inner_join(aba_genome_operons_summary %>% select(operon, total_size) %>% unique) %>% 
	select(operon, operon_genes, operon_AB19606, operon_AB030, total_size, essential_size, Pathways) %>% 
	slice(mixedorder(operon)) %>%
	mutate(operon = gsub("-", "_", operon)) %>% 
	slice(mixedorder(operon)) %>% 
	mutate(operon = gsub("_", "-", operon))

operon_details %>% fwrite("operon_details.tsv", sep = "\t")

##########################################################################################

operon_median_results <- melted_results %>% 
	filter(type == "perfect") %>% 
	filter(condition %in% interest$condition) %>%
	left_join(aba_genome_operons_summary) %>% 
	group_by(condition, operon) %>% 
	summarise(operon_mLFC = median(LFC.adj), FDR = stouffer(FDR)$p) 

gene_median_results <- melted_results %>%
	filter(type == "perfect") %>%
	filter(condition %in% interest$condition) %>%
	group_by(condition, unique_name, AB19606, AB030) %>% 
	summarise(gene_mLFC = median(LFC.adj), FDR = stouffer(FDR)$p)

melted_results %>% 
	filter(type == "perfect") %>% 
	filter(condition %in% interest$condition) %>%
	left_join(aba_genome_operons_summary) %>% 
	inner_join(aba_genome %>% select(left, right, spacer), by = "spacer") 

aba_genome_operons_summary <- melted_results %>% 
	filter(type == "perfect") %>% 
	filter(condition %in% interest$condition) %>%
	left_join(aba_genome_operons_summary) %>% 
	inner_join(aba_genome %>% select(left, right, spacer, strand), by = "spacer") %>% 
	group_by(operon, strand) %>% 
	summarise(tss = case_when(
		strand == "+" ~ min(left), 
		strand == "-" ~ max(right))) %>% unique %>% 
	inner_join(aba_genome_operons_summary)

imipenem_operons <- melted_results %>% 
	filter(type == "perfect") %>% 
	filter(condition %in% interest$condition) %>%
	left_join(aba_genome_operons_summary) %>% 
	filter(condition %like% "Imipenem" & condition %like% "T1" & condition %like% "0.09" & (Pathway %like% "tRNA" | Pathway %like% "PG") & unique_name != "ftsN") %>% 
	select(operon) %>% 
	unique 

melted_results %>% 
	filter(type == "perfect") %>% 
	filter(condition %in% interest$condition) %>%
	left_join(aba_genome_operons_summary) %>% 
	filter(condition %like% "Imipenem" & condition %like% "T1" & condition %like% "0.09" & operon %in% imipenem_operons$operon & unique_name != "ftsN") %>%
	inner_join(operon_details %>% mutate(operon_genes = gsub(", ", "\n", operon_genes))) %>%
	ggplot(aes(x = tss, y = LFC, group = tss)) + 
	geom_boxplot(width = 50000, aes(fill = Pathway, colour = Pathway), outlier.shape = NA) + 
	geom_abline(slope = 0) + 
	geom_jitter(aes(colour = Pathway), size = 2, alpha = 0.5, height = 0, width = 25000) +
	geom_label_repel(
		force = 5,
		aes(label = operon_genes),
		box.padding = 1.0,
		point.padding = 1.5,
		stat = "summary",
		max.iter = 1000000000,
		min.segment.length = 0,
		segment.curvature = -0.25,
		segment.angle = 90,
		fun = median)  + 
	doc_theme +
	scale_colour_manual(
		values = c(
			"Other" = "black",
			"NADH" = "#6A3D9A",
			"LOS" = "#33A02C",
			"Cell Wall/PG" = "#FF7F00",
			"tRNA Ligase" = "#1F78B4")) +
	scale_fill_manual(
		values = c(
			"Other" = "gray",
			"NADH" = "#CAB2D6",
			"LOS" = "#B2DF8A",
			"Cell Wall/PG" = "#FDBF6F",
			"tRNA Ligase" = "#A6CEE3")) 

##########################################################################################



operon_median_results %>%
	filter(condition %in%  c("None_0_T1 - None_0_T0", "None_0_T2 - None_0_T0")) %>%
	inner_join(operon_pathways %>% filter(Pathway %like% "Ribosome")) %>%
	mutate(Pathway = case_when(
		Pathway == "Ribosome" ~ "Ribosome",
		TRUE ~ NA_character_)) %>%
	ggplot(
		aes(x = operon_mLFC,
				y = FDR,
				colour = Pathway)) +
	geom_point(aes(size = operon_pathway_size)) +
	# geom_point(data = . %>% filter(Pathway == "Ribosome"), size = 3) +
	geom_hline(yintercept = 0.05,
						 linetype = "dashed",
						 color = "red",
						 lwd = 1) +
	geom_vline(xintercept = -1,
						 linetype = "dashed",
						 color = "red", 
						 lwd = 1) +
	geom_vline(xintercept = 0,
						 linetype = "solid",
						 color = "black",
						 lwd = 1) +
	doc_theme +
	scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
	scale_colour_manual(values = c("dark red"), na.value = "grey") +
	# theme(legend.position = "none") +
	geom_text_repel(
		data = . %>% filter(Pathway == "Ribosome" & FDR < 0.05 & operon_mLFC < -1),
		aes(label = operon),
		min.segment.length = 0,
		box.padding = 1,
		point.padding = .25,
		# parse = TRUE,
		max.overlaps = Inf,
		colour = "black") +
	facet_grid(~condition) -> to_plot

print(to_plot)

