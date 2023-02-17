source("meta_counter.R")

melted_results[
	unique_name %like% "^sdh" | unique_name %like% "^atp" | unique_name %like% "^cyo" | unique_name %like% "^nuo",
	Pathway := "Ox Phos"]

# Load several packages from CRAN and Bioconductor
require('pacman')
p_load(
	conflicted,
	data.table,
	scales,
	edgeR,
	statmod,
	poolr,
	pheatmap,
	svglite,
	ggplot2,
	ggrepel,
	Rtsne,
	pracma,
	colourpicker,
	RColorBrewer,
	vegan,
	tidyverse,
	magrittr,
	ggallin,
	gtools)

doc_theme <- theme_ipsum(base_family = "Arial", caption_margin = 12, axis_title_size = 12, axis_col = "black")

p_load_current_gh("hrbrmstr/hrbrthemes")

conflict_prefer("extract", "magrittr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Unfortunately, it's hard to anticipate lots of kinds of missingness.

operon_conversion <- fread("Operons/AB19606_operon_conversion.tsv", na.strings = "")
operon_conversion <- operon_conversion %>% 
	rbind(curated_names %>% filter(!AB19606 %in% operon_conversion$locus_tag) %>% select(AB19606) %>% rename(locus_tag = AB19606), fill = TRUE)
operon_conversion <- operon_conversion %>% slice(mixedorder(locus_tag)) %>% filter(!is.na(locus_tag))
operon_conversion[is.na(operon), operon := locus_tag]

operon_details <- curated_names %>% 
	inner_join(
		aba_genome %>% rename("AB19606" = "locus_tag") %>%
			select(left, right, strand, AB19606) %>% unique) %>%
	inner_join(operon_conversion %>% rename("AB19606" = "locus_tag")) %>%
	arrange(left) %>%
	group_by(operon) %>%
	summarise(
		strand = unique(strand),
		tss = case_when(
			strand == "+" ~ min(left), 
			strand == "-" ~ max(right)),
		operon_genes = case_when(
			strand == "+" ~ paste(unique_name, collapse = ", "), 
			strand == "-" ~ paste(rev(unique_name), collapse = ", ")),
		operon_AB19606 = case_when(
			strand == "+" ~ paste(AB19606, collapse = ", "), 
			strand == "-" ~ paste(rev(AB19606), collapse = ", ")),			
		operon_AB030 = case_when(
			strand == "+" ~ paste(AB030, collapse = ", "), 
			strand == "-" ~ paste(rev(AB030), collapse = ", ")),
		essential_size = n())	%>% 
	inner_join(operon_conversion %>% group_by(operon) %>% tally(name = "total_size")) %>%
	arrange(tss)

aba_genome_operons <-
	aba_genome %>% 
	full_join(operon_conversion) 

operon_pathways <- melted_results %>% 
	filter(type == "perfect") %>%
	select(AB19606, Pathway) %>% 
	rename("locus_tag" = "AB19606") %>% 
	unique %>% 
	inner_join(aba_genome_operons) %>% 
	unique %>% 
	select(operon, Pathway, locus_tag) %>% 
	unique %>% 
	mutate(
		Pathway = factor(Pathway), 
		Pathway = factor(Pathway, levels = c(levels(Pathway)[levels(Pathway) != "Other"], "Other"))) %>%
	slice(mixedorder(Pathway)) %>%
	group_by(operon) %>% 
	summarise(Pathways = paste(unique(Pathway), collapse = "+"))

operon_details <- operon_details %>% inner_join(operon_pathways) %>% arrange(tss)

operon_details %>% fwrite("operon_details.tsv", sep = "\t")

operon_details <- fread("operon_details.tsv")

operon_details <- operon_details %>%
	mutate(
		Pathways = factor(Pathways),
		Pathways = factor(Pathways, levels = c(levels(Pathways)[levels(Pathways) != "Other"], "Other")))

##########################################################################################

operon_median_results <- melted_results %>% 
	filter(type == "perfect") %>% 
	left_join(operon_conversion, by = c("AB19606" = "locus_tag")) %>% 
	filter(condition %in% interest$condition) %>%
	group_by(condition, operon) %>% 
	summarise(operon_mLFC = median(LFC.adj), FDR = stouffer(FDR)$p) 

gene_median_results <- melted_results %>%
	filter(type == "perfect") %>%
	filter(condition %in% interest$condition) %>%
	group_by(condition, unique_name, AB19606, AB030) %>% 
	summarise(gene_mLFC = median(LFC.adj), FDR = stouffer(FDR)$p)

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
			"Ox Phos" = "#6A3D9A",
			"LOS" = "#33A02C",
			"Cell Wall/PG" = "#FF7F00",
			"tRNA Ligase" = "#1F78B4")) +
	scale_fill_manual(
		values = c(
			"Other" = "gray",
			"Ox Phos" = "#CAB2D6",
			"LOS" = "#B2DF8A",
			"Cell Wall/PG" = "#FDBF6F",
			"tRNA Ligase" = "#A6CEE3")) 

##########################################################################################



operon_median_results %>%
	filter(condition %in% c("None_0_T1 - None_0_T0", "None_0_T2 - None_0_T0")) %>%
	group_by(condition) %>% 
	arrange(operon_mLFC) %>% mutate(operon_mLFC_ix = row_number())%>%
	arrange(FDR) %>% mutate(FDR_ix = row_number()) %>%
	inner_join(operon_details) %>%
	mutate(FDR = case_when(FDR != 0 ~ FDR, FDR == 0 ~ min(FDR[FDR != 0]))) %>%
	mutate(`Transcription Unit` = case_when(essential_size == 1 ~ operon_genes, TRUE ~ operon_genes)) %>%
	ggplot(
		aes(x = operon_mLFC,
				y = FDR,
				colour = Pathways)) +
	geom_point(aes(size = essential_size)) +
	geom_hline(yintercept = 0.05,
						 linetype = "dashed",
						 color = "dark grey",
						 lwd = 1) +
	geom_vline(xintercept = -1,
						 linetype = "dashed",
						 color = "dark grey", 
						 lwd = 1) +
	geom_vline(xintercept = 0,
						 linetype = "solid",
						 color = "dark grey",
						 lwd = 1) +
	doc_theme +
	scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
	geom_text_repel(
		max.iter = 1000000000,
		data = . %>% filter(operon_mLFC_ix <= 15 | FDR_ix <= 10),
		aes(label = `Transcription Unit`),
		force = 7.5,
		min.segment.length = 0,
		box.padding = 2,
		point.padding = .25,
		# parse = TRUE,
		size = 3,
		max.overlaps = Inf,
		colour = "black") +
	facet_wrap(~ condition, scales = "free") +
	scale_colour_manual(
		values = c(
			"Other" = "gray",
			"Ribosome" = "#E31A1C",
			"Ribosome+Other" = "#FB9A99",
			"Ox Phos" = "#6A3D9A",
			"Ox Phos+Other" = "#CAB2D6",
			"LOS" = "#33A02C",
			"LOS+Other" = "#B2DF8A",
			"Cell Wall/PG" = "#FF7F00",
			"Cell Wall/PG+Other" = "#FDBF6F",
			"tRNA Ligase" = "#1F78B4",
			"tRNA Ligase+Other" = "#A6CEE3")) 


operon_median_results %>%
	filter(condition %in% c("Colistin_0.44_T1 - None_0_T1", "Colistin_0.44_T2 - None_0_T2")) %>%
	group_by(condition) %>% 
	arrange(operon_mLFC) %>% mutate(operon_mLFC_ix = row_number()) %>%
	arrange(desc(operon_mLFC)) %>% mutate(operon_mLFC_desc_ix = row_number()) %>%
	arrange(FDR) %>% mutate(FDR_ix = row_number()) %>%
	inner_join(operon_details) %>%
	mutate(FDR = case_when(FDR != 0 ~ FDR, FDR == 0 ~ min(FDR[FDR != 0]))) %>%
	mutate(`Transcription Unit` = case_when(essential_size == 1 ~ operon_genes, TRUE ~ operon_genes)) %>%
	ggplot(
		aes(x = operon_mLFC,
				y = FDR,
				colour = Pathways)) +
	geom_point(aes(size = essential_size)) +
	geom_hline(yintercept = 0.05,
						 linetype = "dashed",
						 color = "dark grey",
						 lwd = 1) +
	# geom_vline(xintercept = -1,
	# 					 linetype = "dashed",
	# 					 color = "dark grey", 
	# 					 lwd = 1) +
	geom_vline(xintercept = 0,
						 linetype = "solid",
						 color = "dark grey",
						 lwd = 1) +
	doc_theme +
	scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
	geom_text_repel(
		max.iter = 1000000000,
		data = . %>% filter(FDR < 0.05 & (operon_mLFC_ix <= 10 & operon_mLFC < -0.5 | (operon_mLFC_desc_ix <= 10 & operon_mLFC > 0.5 ) | FDR_ix <= 10)),
		aes(label = `Transcription Unit`),
		force = 7.5,
		min.segment.length = 0,
		box.padding = 2,
		point.padding = .25,
		# parse = TRUE,
		size = 3,
		max.overlaps = Inf,
		colour = "black") +
	facet_wrap(~ condition, scales = "free") +
	scale_colour_manual(
		values = c(
			"Other" = "gray",
			"Ribosome" = "#E31A1C",
			"Ribosome+Other" = "#FB9A99",
			"Ox Phos" = "#6A3D9A",
			"Ox Phos+Other" = "#CAB2D6",
			"LOS" = "#33A02C",
			"LOS+Other" = "#B2DF8A",
			"Cell Wall/PG" = "#FF7F00",
			"Cell Wall/PG+Other" = "#FDBF6F",
			"tRNA Ligase" = "#1F78B4",
			"tRNA Ligase+Other" = "#A6CEE3")) 


operon_median_results %>%
	filter(condition %in% c("Rifampicin_0.34_T1 - None_0_T1", "Rifampicin_0.34_T2 - None_0_T2")) %>%
	group_by(condition) %>% 
	arrange(operon_mLFC) %>% mutate(operon_mLFC_ix = row_number()) %>%
	arrange(desc(operon_mLFC)) %>% mutate(operon_mLFC_desc_ix = row_number()) %>%
	arrange(FDR) %>% mutate(FDR_ix = row_number()) %>%
	inner_join(operon_details) %>%
	mutate(FDR = case_when(FDR != 0 ~ FDR, FDR == 0 ~ min(FDR[FDR != 0]))) %>%
	mutate(`Transcription Unit` = case_when(essential_size == 1 ~ operon_genes, TRUE ~ operon_genes)) %>%
	ggplot(
		aes(x = operon_mLFC,
				y = FDR,
				colour = Pathways)) +
	geom_point(aes(size = essential_size)) +
	geom_hline(yintercept = 0.05,
						 linetype = "dashed",
						 color = "dark grey",
						 lwd = 1) +
	# geom_vline(xintercept = -1,
	# 					 linetype = "dashed",
	# 					 color = "dark grey", 
	# 					 lwd = 1) +
	geom_vline(xintercept = 0,
						 linetype = "solid",
						 color = "dark grey",
						 lwd = 1) +
	doc_theme +
	scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
	geom_text_repel(
		max.iter = 1000000000,
		data = . %>% filter(FDR < 0.05 & (operon_mLFC_ix <= 10 & operon_mLFC < -0.5 | (operon_mLFC_desc_ix <= 10 & operon_mLFC > 0.5 ) | FDR_ix <= 10)),
		aes(label = `Transcription Unit`),
		force = 7.5,
		min.segment.length = 0,
		box.padding = 2,
		point.padding = .25,
		# parse = TRUE,
		size = 3,
		max.overlaps = Inf,
		colour = "black") +
	facet_wrap(~ condition, scales = "free") +
	scale_colour_manual(
		values = c(
			"Other" = "gray",
			"Ribosome" = "#E31A1C",
			"Ribosome+Other" = "#FB9A99",
			"Ox Phos" = "#6A3D9A",
			"Ox Phos+Other" = "#CAB2D6",
			"LOS" = "#33A02C",
			"LOS+Other" = "#B2DF8A",
			"Cell Wall/PG" = "#FF7F00",
			"Cell Wall/PG+Other" = "#FDBF6F",
			"tRNA Ligase" = "#1F78B4",
			"tRNA Ligase+Other" = "#A6CEE3")) 



operon_median_results %>%
	filter(condition %in% c("Imipenem_0.09_T1 - None_0_T1", "Imipenem_0.09_T2 - None_0_T2")) %>%
	group_by(condition) %>% 
	arrange(operon_mLFC) %>% mutate(operon_mLFC_ix = row_number()) %>%
	arrange(desc(operon_mLFC)) %>% mutate(operon_mLFC_desc_ix = row_number()) %>%
	arrange(FDR) %>% mutate(FDR_ix = row_number()) %>%
	inner_join(operon_details) %>%
	mutate(FDR = case_when(FDR != 0 ~ FDR, FDR == 0 ~ min(FDR[FDR != 0]))) %>%
	mutate(`Transcription Unit` = case_when(essential_size == 1 ~ operon_genes, TRUE ~ operon_genes)) %>%
	ggplot(
		aes(x = operon_mLFC,
				y = FDR,
				colour = Pathways)) +
	geom_point(aes(size = essential_size)) +
	geom_hline(yintercept = 0.05,
						 linetype = "dashed",
						 color = "dark grey",
						 lwd = 1) +
	# geom_vline(xintercept = -1,
	# 					 linetype = "dashed",
	# 					 color = "dark grey", 
	# 					 lwd = 1) +
	geom_vline(xintercept = 0,
						 linetype = "solid",
						 color = "dark grey",
						 lwd = 1) +
	doc_theme +
	scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
	geom_text_repel(
		max.iter = 1000000000,
		data = . %>% filter(FDR < 0.05 & (operon_mLFC_ix <= 10 & operon_mLFC < -0.5 | (operon_mLFC_desc_ix <= 10 & operon_mLFC > 0.5 ) | FDR_ix <= 10)),
		aes(label = `Transcription Unit`),
		force = 7.5,
		min.segment.length = 0,
		box.padding = 2,
		point.padding = .25,
		# parse = TRUE,
		size = 3,
		max.overlaps = Inf,
		colour = "black") +
	facet_wrap(~ condition, scales = "free") +
	scale_colour_manual(
		values = c(
			"Other" = "gray",
			"Ribosome" = "#E31A1C",
			"Ribosome+Other" = "#FB9A99",
			"Ox Phos" = "#6A3D9A",
			"Ox Phos+Other" = "#CAB2D6",
			"LOS" = "#33A02C",
			"LOS+Other" = "#B2DF8A",
			"Cell Wall/PG" = "#FF7F00",
			"Cell Wall/PG+Other" = "#FDBF6F",
			"tRNA Ligase" = "#1F78B4",
			"tRNA Ligase+Other" = "#A6CEE3")) 


operon_median_results %>%
	filter(condition %in% c("Meropenem_0.17_T1 - None_0_T1", "Meropenem_0.17_T2 - None_0_T2")) %>%
	group_by(condition) %>% 
	arrange(operon_mLFC) %>% mutate(operon_mLFC_ix = row_number()) %>%
	arrange(desc(operon_mLFC)) %>% mutate(operon_mLFC_desc_ix = row_number()) %>%
	arrange(FDR) %>% mutate(FDR_ix = row_number()) %>%
	inner_join(operon_details) %>%
	mutate(FDR = case_when(FDR != 0 ~ FDR, FDR == 0 ~ min(FDR[FDR != 0]))) %>%
	mutate(`Transcription Unit` = case_when(essential_size == 1 ~ operon_genes, TRUE ~ operon_genes)) %>%
	ggplot(
		aes(x = operon_mLFC,
				y = FDR,
				colour = Pathways)) +
	geom_point(aes(size = essential_size)) +
	geom_hline(yintercept = 0.05,
						 linetype = "dashed",
						 color = "dark grey",
						 lwd = 1) +
	# geom_vline(xintercept = -1,
	# 					 linetype = "dashed",
	# 					 color = "dark grey", 
	# 					 lwd = 1) +
	geom_vline(xintercept = 0,
						 linetype = "solid",
						 color = "dark grey",
						 lwd = 1) +
	doc_theme +
	scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
	geom_text_repel(
		max.iter = 1000000000,
		data = . %>% filter(FDR < 0.05 & (operon_mLFC_ix <= 10 & operon_mLFC < -0.5 | (operon_mLFC_desc_ix <= 10 & operon_mLFC > 0.5 ) | FDR_ix <= 10)),
		aes(label = `Transcription Unit`),
		force = 7.5,
		min.segment.length = 0,
		box.padding = 2,
		point.padding = .25,
		# parse = TRUE,
		size = 3,
		max.overlaps = Inf,
		colour = "black") +
	facet_wrap(~ condition, scales = "free") +
	scale_colour_manual(
		values = c(
			"Other" = "gray",
			"Ribosome" = "#E31A1C",
			"Ribosome+Other" = "#FB9A99",
			"Ox Phos" = "#6A3D9A",
			"Ox Phos+Other" = "#CAB2D6",
			"LOS" = "#33A02C",
			"LOS+Other" = "#B2DF8A",
			"Cell Wall/PG" = "#FF7F00",
			"Cell Wall/PG+Other" = "#FDBF6F",
			"tRNA Ligase" = "#1F78B4",
			"tRNA Ligase+Other" = "#A6CEE3")) 


