median_melted_results %>% 
	inner_join(gene_pathways %>% rename(AB19606 = locus_tag)) %>%
	filter(condition %in% c("Imipenem_0.09_T1 - None_0_T1") & type == "perfect") %>%
	mutate(Pathway = case_when(
		`Gene Pathway` %like% "tRNA" ~ "tRNA Ligase",
		`Gene Pathway` %like% "PG" ~ "Cell Wall/Division",
		# `Gene Pathway` %like% "Ox Phos" ~ "Ox Phos",
		TRUE ~ NA_character_),
		Significance = ifelse((medLFC < -1 & FDR < 0.05) | (medLFC > 1 & FDR < 0.05), "Significant", "Not Significant")) %>%
	ggplot(aes(x = medLFC, y = FDR, colour = `Pathway`, fill = `Pathway`, alpha = Significance)) +
	geom_point(shape = 21, size = 5, 
						 data = . %>% filter(!Pathway %in% c("tRNA Ligase", "Cell Wall/Division", "Ox Phos"))) +
	geom_point(shape = 21, size = 5, 
						 data = . %>% filter(Pathway %in% c("tRNA Ligase", "Cell Wall/Division", "Ox Phos"))) +
	geom_hline(yintercept = 0.05,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = -1,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = 1,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = 0,
						 linetype = "solid",
						 color = "black",
						 lwd = 1) +
	doc_theme +
	scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
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
			"tRNA Ligase+Other" = "#A6CEE3"),
		na.value = "grey") +
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
			"tRNA Ligase+Other" = "#A6CEE3"),
		na.value = "grey") +
	scale_alpha_manual(
		values = c(
			"Significant" = 0.85,
			"Not Significant" = 0.15)) +
	geom_text_repel(data = . %>% filter(abs(medLFC) >= 1 & FDR <= 0.05 & !is.na(Pathway) ), 
									aes(label = unique_name),
									segment.color = 'black', 
									segment.size = 0.5, 
									min.segment.length = 0.25,
									force = 25,
									segment.linetype = 'solid',
									colour = "black") +
	
	
	facet_wrap(~condition) -> to_plot

print(to_plot)
