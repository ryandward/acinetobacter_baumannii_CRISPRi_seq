source("operon_analytics.R")

library(ggplot2)
library(dplyr)
library(tidyr)

gene_pathways <- melted_results %>% 
	filter(type == "perfect") %>%
	select(AB19606, Pathway) %>% 
	rename("locus_tag" = "AB19606") %>% 
	unique %>% 
	rename(`Gene Pathway` = Pathway) %>% 
	inner_join(curated_names %>% rename(locus_tag = AB19606)) %>% select(locus_tag, `Gene Pathway`, unique_name)

gene_pathways[unique_name %like% "nuo" | unique_name %like% "cyo" | unique_name %like% "atp" | unique_name %like% "sdh", `Gene Pathway` := "Ox Phos"]

operon_details %>% 
	separate_rows(operon_AB19606, sep = ", ") %>%
	rename(AB19606 = operon_AB19606) %>%
	select(operon, strand, tss, condensed_name, long_name, Pathways, AB19606) %>%
	inner_join(melted_results %>% filter(type == "perfect") %>% filter(condition == "Imipenem_0.09_T1 - None_0_T1")) %>% 
	inner_join(gene_pathways %>% rename(AB19606 = locus_tag))  %>%
	mutate(direct_target = case_when(unique_name == 'ftsI' ~ TRUE, TRUE ~ FALSE)) %>%
	filter(Pathways %like% "tRNA" | Pathways %like% "PG" | long_name == "dapA") %>%
	mutate(Pathways = gsub("Cell Wall/PG", "Cell Wall/Division", Pathways)) %>%
	mutate(`Gene Pathway` = gsub("Cell Wall/PG", "Cell Wall/Division", `Gene Pathway`)) %>%
	select(long_name, Pathways, LFC.adj, tss, `Gene Pathway`, direct_target, FDR) %>%
	mutate(tss = factor(tss)) %>%
	ggplot(aes(x = tss, y = LFC.adj, fill = Pathways)) +
	geom_boxplot(alpha = 0.35, outlier.shape = NA) +
	doc_theme +
	theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
	labs(x = "TSS", y = "LFC.adj", fill = "Pathways") +
	geom_jitter(alpha = 0.5, aes(colour = `Gene Pathway`, size = FDR)) +
	scale_size_continuous(trans = reverse_trans() %of% log10_trans(),
												breaks = trans_breaks("log10", function(x) 10^-x),
												labels = trans_format("log10", math_format(10^.x))) +
	
	geom_abline(slope = 0) + 
	scale_fill_manual(
		values = c(
			"Other" = "dark grey",
			"Ribosome" = "#E31A1C",
			"Ribosome+Other" = "#FB9A99",
			"Ox Phos" = "#6A3D9A",
			"Ox Phos+Other" = "#CAB2D6",
			"LOS" = "#33A02C",
			"LOS+Other" = "#B2DF8A",
			"Cell Wall/Division" = "#FF7F00",
			"Cell Wall/Division+Other" = "#FDBF6F",
			"tRNA Ligase" = "#1F78B4",
			"tRNA Ligase+Other" = "#A6CEE3")) +
	scale_color_manual(
		values = c(
			"Other" = "dark grey",
			"Ribosome" = "#E31A1C",
			"Ribosome+Other" = "#FB9A99",
			"Ox Phos" = "#6A3D9A",
			"Ox Phos+Other" = "#CAB2D6",
			"LOS" = "#33A02C",
			"LOS+Other" = "#B2DF8A",
			"Cell Wall/Division" = "#FF7F00",
			"Cell Wall/Division+Other" = "#FDBF6F",
			"tRNA Ligase" = "#1F78B4",
			"tRNA Ligase+Other" = "#A6CEE3")) +
	geom_label_repel(
		force = 5,
		fill = alpha("white", 0.5),
		color = "black",
		aes(label = long_name),
		box.padding = 1.0,
		point.padding = 1.5,
		stat = "summary",
		max.iter = 1000000000,
		min.segment.length = 0,
		segment.curvature = -0.25,
		# segment.angle = 90,
		fun = median) 
	