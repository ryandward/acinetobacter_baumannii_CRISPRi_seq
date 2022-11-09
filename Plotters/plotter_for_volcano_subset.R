library(pacman)

p_load(tidyverse, ggplot2, data.table, ggrepel, hrbrthemes, viridis, ggallin)

#https://cran.r-project.org/web/packages/ggallin/README.html

melted_results <- fread("Results/melted_results.tsv.gz", sep = "\t")
median_melted_results <- fread("Results/median_melted_results.tsv.gz", sep = "\t")

interest <- fread("interest.tsv", sep = "\t")

doc_theme <- theme_ipsum(
	base_family = "Arial", 
	caption_margin = 12,
	axis_title_size = 12,
	axis_col = "black")

median_melted_results %>%
	filter(condition %in% c("None_0_T1 - None_0_T0", "None_0_T2 - None_0_T0") & type == "perfect") %>%
	mutate(condition = case_when(
		condition %like% "T1" ~ "T1",
		condition %like% "T2" ~ "T2"
	)) %>%
	filter(condition %like% "T1") %>%
	mutate(Pathway = case_when(
		Pathway == "Ribosome" ~ "Ribosome",
		Pathway == "LOS" ~ "LOS",
		unique_name %like% "nuo" ~ "NDH-1",
		Pathway %like% "Cell Wall" ~ "PG/Division",
		TRUE ~ "Other")) %>%
	ggplot(
		aes(x = medLFC,
				y = FDR)) +
	geom_point(
		data = . %>% filter(FDR >= 0.05 | abs(medLFC) <= 1),
		colour = "light grey") +
	geom_point(
		data = . %>% filter(FDR < 0.05 & abs(medLFC) > 1),
		aes(colour = Pathway), size = 3, alpha = 0.85) +
	geom_hline(yintercept = 0.05,
						 linetype = "dashed",
						 color = "black",
						 lwd = 1) +
	geom_vline(xintercept = -1,
						 linetype = "dashed",
						 color = "black", 
						 lwd = 1) +
	geom_vline(xintercept = 0,
						 linetype = "solid",
						 color = "black",
						 lwd = 1) +
	# geom_text_repel(
	# 	data = . %>% filter(unique_name == "GO593_00515"),
	# 	aes(label = gene_name_stylized),
	# 	min.segment.length = 0,
	# 	parse = TRUE,
	# 	max.overlaps = Inf,
	# 	colour = "black") +
	doc_theme +
	scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
	scale_colour_manual(
		values = c(
			"Ribosome" = "#FB9A99",
			"NDH-1" = "#CAB2D6",
			"LOS" = "#B2DF8A",
			"PG/Division" = "#FDBF6F",
			"Other" = "#808080")) +
	# scale_alpha_manual(
	# 	values = c(
	# 		"T1.Ribosome" = 1, 
	# 		"T2.Ribosome" = 1,
	# 		"T1.NDH-1" = 1,
	# 		"T2.NDH-1" = 1,
	# 		"T1.LOS" = 1,
	# 		"T2.LOS" = 1,
	# 		"T1.PG/Division" = 1,
	# 		"T2.PG/Division" = 1,
	# 		na.value = 0.25)) +
	theme(legend.position = "bottom") +
	facet_wrap(~condition, scales = "free") -> to_plot

print(to_plot)


#T1 Ribosomes


median_melted_results %>%
	filter(condition ==  "None_0_T2 - None_0_T0" & type == "perfect") %>%
	mutate(Pathway = case_when(
		Pathway == "Ribosome" ~ "Ribosome",
		TRUE ~ NA_character_)) %>%
	ggplot(
		aes(x = medLFC,
				y = FDR,
				colour = Pathway)) +
	geom_point() +
	geom_point(data = . %>% filter(Pathway == "Ribosome"), size = 3) +
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
	theme(legend.position = "none") +
	geom_text_repel(
		data = . %>% filter(Pathway == "Ribosome" & FDR < 0.05 & medLFC < -1),
		aes(label = gene_name_stylized),
		min.segment.length = 0,
		parse = TRUE,
		max.overlaps = Inf,
		colour = "black") +
	theme(legend.position="none")-> to_plot

print(to_plot)

	
	#T1 nuo
	
median_melted_results %>%
	filter(condition ==  "None_0_T2 - None_0_T0" & type == "perfect") %>%
	mutate(Pathway = case_when(
		Pathway == "NADH" ~ "NADH",
		TRUE ~ NA_character_)) %>%
	ggplot(
		aes(x = medLFC,
				y = FDR,
				colour = Pathway)) +
	geom_point() +
	geom_point(data = . %>% filter(Pathway == "NADH"), size = 3) +
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
	scale_colour_manual(values = c("dark cyan"), na.value = "grey") +
	theme(legend.position = "none") +
	geom_text_repel(
		data = . %>% filter(Pathway == "NADH" & FDR < 0.05),
		aes(label = gene_name_stylized),
		min.segment.length = 0,
		parse = TRUE,
		max.overlaps = Inf,
		colour = "black")+
	theme(legend.position="none") -> to_plot

print(to_plot)

#T2 colistin NADH

median_melted_results %>%
	filter(condition ==  "Colistin_0.44_T2 - None_0_T2" & type == "perfect") %>%
	mutate(Pathway = case_when(
		# Pathway == "NADH" ~ "NADH",
		TRUE ~ NA_character_)) %>%
	ggplot(
		aes(x = medLFC,
				y = FDR,
				colour = Pathway)) +
	geom_point() +
	geom_point(data = . %>% filter(Pathway == "NADH"), size = 3) +
	geom_hline(yintercept = 0.05,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = 0,
						 linetype = "dotted",
						 color = "black") +
	doc_theme +
	scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
	scale_colour_manual(values = c("dark cyan"), na.value = "grey") +
	theme(legend.position = "none") +
	geom_text_repel(
		data = . %>% filter(Pathway == "NADH"),
		aes(label = gene_name_stylized),
		parse = TRUE,
		colour = "black") -> to_plot

print(to_plot)

#T2 colistin LOS

median_melted_results %>%
	filter(condition ==  "Colistin_0.44_T2 - None_0_T2" & type == "perfect") %>%
	mutate(Pathway = case_when(
		Pathway == "LOS" & !(unique_name %like% "lpt")  ~ "LOS",
		TRUE ~ NA_character_)) %>%
	ggplot(
		aes(x = medLFC,
				y = FDR,
				colour = Pathway)) +
	geom_point() +
	geom_point(data = . %>% filter(Pathway == "LOS"), size = 3) +
	geom_hline(yintercept = 0.05,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = 0,
						 linetype = "dotted",
						 color = "black") +
	doc_theme +
	scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
	scale_colour_manual(values = c("purple"), na.value = "grey") +
	theme(legend.position = "none") +
	geom_text_repel(
		data = . %>% filter(Pathway == "LOS"),
		aes(label = gene_name_stylized),
		min.segment.length = 0,
		parse = TRUE,
		max.overlaps = Inf,
		colour = "black") -> to_plot

print(to_plot)


#T2 rifampicin NADH

median_melted_results %>%
	filter(condition ==  "Rifampicin_0.34_T2 - None_0_T2" & type == "perfect") %>%
	mutate(Pathway = case_when(
		Pathway == "NADH" ~ "NADH",
		TRUE ~ NA_character_)) %>%
	ggplot(
		aes(x = medLFC,
				y = FDR,
				colour = Pathway)) +
	geom_point() +
	geom_point(data = . %>% filter(Pathway == "NADH"), size = 3) +
	geom_hline(yintercept = 0.05,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = 0,
						 linetype = "dotted",
						 color = "black") +
	doc_theme +
	scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
	scale_colour_manual(values = c("dark cyan"), na.value = "grey") +
	theme(legend.position = "none") +
	geom_text_repel(
		data = . %>% filter(Pathway == "NADH"),
		aes(label = gene_name_stylized),
		parse = TRUE,
		colour = "black") -> to_plot

print(to_plot)

#T2 rifampicin LOS

median_melted_results %>%
	filter(condition ==  "Rifampicin_0.34_T2 - None_0_T2" & type == "perfect") %>%
	mutate(Pathway = case_when(
		Pathway == "LOS" & !(unique_name %like% "lpt")  ~ "LOS",
		TRUE ~ NA_character_)) %>%
	ggplot(
		aes(x = medLFC,
				y = FDR,
				colour = Pathway)) +
	geom_point() +
	geom_point(data = . %>% filter(Pathway == "LOS"), size = 3) +
	geom_hline(yintercept = 0.05,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = 0,
						 linetype = "dotted",
						 color = "black") +
	doc_theme +
	scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
	scale_colour_manual(values = c("purple"), na.value = "grey") +
	theme(legend.position = "none") +
	geom_text_repel(
		data = . %>% filter(Pathway == "LOS"),
		aes(label = gene_name_stylized),
		min.segment.length = 0,
		parse = TRUE,
		max.overlaps = Inf,
		colour = "black") -> to_plot

print(to_plot)




################################################################################
# cell shape, high dose, carbapenem drugs 
CL707 <- fread("CL707_by_name.tsv")

median_melted_results %>%
	filter(condition %in% interest$condition) %>%
	filter(condition %like% "Meropenem" | condition %like% "Imipenem") %>%
	filter(condition %like% "0.09" | condition %like% "0.17") %>%
	filter(type == "perfect") %>%
	mutate(Pathway = case_when(
		unique_name %in% CL707$unique_name & FDR < 0.05 & abs(medLFC) > 1 ~ "Cell shape",
		TRUE ~ NA_character_)) %>%
	ggplot(
		aes(x = medLFC,
				y = FDR,
				colour = Pathway)) +
	geom_point(data = . %>% filter(is.na(Pathway))) +
	geom_point(data = . %>% filter(!is.na(Pathway)), size = 2) +
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
	geom_label_repel(
		data = . %>% filter(Pathway == "Cell shape"),
		aes(label = gene_name_stylized),
		min.segment.length = 0,
		label.size = 0.15,
		parse = TRUE,
		max.overlaps = Inf,
		colour = "black") +
	doc_theme +
	scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
	scale_colour_manual(values = c("#E31A1C", "#FB9A99"), na.value = "grey") +
	theme(legend.position = "none") +
	theme(legend.position = "none") +
	facet_wrap(~condition) -> to_plot

print(to_plot)


################################################################################
# Colistin/Rifampcin

median_melted_results %>%
	filter(condition %in% interest$condition) %>%
	filter(condition %like% "Colistin" | condition %like% "Rifampicin") %>%
	filter(type == "perfect") %>%
	filter(condition %like% "T2") %>%
	mutate(Pathway = case_when(
		Pathway %like% "LOS" ~ "LOS",
		unique_name %like% "nuo" ~ "NDH-1",
		TRUE ~ NA_character_)) %>%
	ggplot(
		aes(x = medLFC,
				y = FDR,
				colour = Pathway)) +
	geom_point(data = . %>% filter(is.na(Pathway))) +
	geom_point(data = . %>% filter(!is.na(Pathway)), size = 2) +
	geom_hline(yintercept = 0.05,
						 linetype = "dashed",
						 color = "black",
						 lwd = 1) +
	geom_vline(xintercept = -1,
						 linetype = "dashed",
						 color = "black", 
						 lwd = 1) +
	geom_vline(xintercept = 1,
						 linetype = "dashed",
						 color = "black", 
						 lwd = 1) +
	geom_vline(xintercept = 0,
						 linetype = "solid",
						 color = "black",
						 lwd = 1) +
	geom_label_repel(
		data = . %>% filter(Pathway == "LOS" | Pathway == "NDH-1"),
		aes(label = gene_name_stylized),
		min.segment.length = 0,
		label.size = 0.15,
		parse = TRUE,
		max.overlaps = Inf,
		colour = "black") +
	doc_theme +
	scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
	scale_colour_manual(values = c("#33A02C", "#6A3D9A"), na.value = "grey") +
	facet_wrap(~condition, scales = "free") -> to_plot

print(to_plot)


