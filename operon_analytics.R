# source("meta_counter.R")

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
	pracma,
	RColorBrewer,
	tidyverse,
	magrittr,
	ggallin,
	gtools)

p_load_current_gh("hrbrmstr/hrbrthemes")

conflict_prefer("extract", "magrittr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")

doc_theme <- theme_ipsum(base_family = "Arial", caption_margin = 12, axis_title_size = 12, axis_col = "black")

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
		`Essentials in operon` = n())	%>% 
	inner_join(operon_conversion %>% group_by(operon) %>% tally(name = "total_size")) %>%
	arrange(tss)


##########################################################################################
# excessively hard way to find the condensed names of an operon

operonize = function(x, .sep = "")
{
	concat = paste(x, collapse = .sep)
	substring(concat, 1L, cumsum(c(nchar(x[[1L]]), nchar(x[-1L]) + nchar(.sep))))
}

short_names <- curated_names %>% 
	inner_join(
		aba_genome %>% rename("AB19606" = "locus_tag") %>%
			select(left, right, strand, AB19606) %>% unique) %>%
	separate(unique_name, c("prefix", "suffix"), sep = "(?<=[a-z]{3})") %>%	
	mutate(suffix = case_when(is.na(suffix) ~ "", TRUE ~ suffix)) %>% 
	inner_join(operon_conversion %>% rename("AB19606" = "locus_tag")) %>%
	group_by(operon) %>% 
	arrange(case_when(strand == "+" ~ left, strand == "-" ~ desc(right)))


# Convert data to data.table
short_names <- data.table(short_names)

# Calculate the lengths of runs with identical operon and prefix
groups <- rle(paste(short_names$operon, short_names$prefix, sep = "_"))$lengths

# Add the group numbers to the data table
short_names[, group := rep(seq_along(groups), groups)]

# View the result
short_names <- short_names %>%
	group_by(group) %>%
	mutate(
		suffices = case_when(
			prefix == lag(prefix) ~ operonize(suffix, ""),
			TRUE ~ suffix)) %>%
	group_by(group) %>% 
	filter(nchar(suffices) == max(nchar(suffices))) %>% 
	mutate(operon_substr = paste0(prefix, suffices)) %>% 
	ungroup %>% group_by(operon) %>% 
	mutate(long_name = paste(operon_substr, collapse = "-")) %>%
	select(operon, long_name) %>%
	unique %>% 
	mutate(condensed_name = str_replace(long_name, "^([^-]+-[^-]+)-.*$", "\\1…")) 


##########################################################################################

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

##########################################################################################

operon_details <- operon_details %>% inner_join(short_names) %>%
	select(operon, strand, tss, long_name, condensed_name, Pathways, operon_genes, operon_AB19606, operon_AB030, `Essentials in operon`, total_size)


##########################################################################################

operon_details %>% fwrite("operon_details.tsv", sep = "\t")

curated_names_operons_pathways <- operon_conversion %>% 
	inner_join(operon_pathways) %>% 
	rename(AB19606 = locus_tag) %>% 
	inner_join(curated_names) %>% 
	inner_join(operon_details)

curated_names_operons_pathways %>% fwrite("curated_names_operons_pathways.tsv", sep = "\t")

operon_details <- fread("operon_details.tsv")

operon_details <- operon_details %>%
	mutate(
		Pathways = factor(Pathways),
		Pathways = factor(Pathways, levels = c(levels(Pathways)[levels(Pathways) != "Other"], "Other")))

##########################################################################################

operon_median_results <- melted_results %>% 
	filter(type == "perfect") %>% 
	left_join(operon_conversion, by = c("AB19606" = "locus_tag")) %>% 
	# filter(condition %in% interest$condition) %>%
	group_by(condition, operon) %>% 
	summarise(operon_mLFC = median(LFC.adj), FDR = stouffer(FDR)$p) %>%
	inner_join(operon_details %>% select(operon, long_name, condensed_name)) %>%
	select(condition, long_name, condensed_name, operon_mLFC, FDR, operon)

operon_median_results %>% fwrite("Results/operon_median_results.tsv", sep = "\t")

gene_median_results <- melted_results %>%
	filter(type == "perfect") %>%
	# filter(condition %in% interest$condition) %>%
	group_by(condition, unique_name, AB19606, AB030) %>% 
	summarise(gene_mLFC = median(LFC.adj), FDR = stouffer(FDR)$p)

gene_median_results %>% fwrite("Results/gene_median_results.tsv", sep = "\t")

##########################################################################################
# function to plot operons

plot_operon <- function(
		operon_median_results, 
		conditions, 
		use_condensed_name = FALSE, 
		max_left = 10, 
		max_right = 10, 
		max_top = 10, 
		min_sig = 1e-100) {
	operon_median_results %>% 
		mutate(FDR = case_when(FDR < min_sig ~ min_sig, TRUE ~ FDR)) %>%
		mutate(Significance = case_when(
			FDR < 0.05 & abs(operon_mLFC) >= 1 ~ "Significant",
			TRUE ~ "Not Significant")) %>%
		filter(condition %in% conditions) %>%
		group_by(condition) %>% 
		arrange(operon_mLFC) %>% mutate(operon_mLFC_ix = row_number()) %>%
		arrange(desc(operon_mLFC)) %>% mutate(operon_mLFC_desc_ix = row_number()) %>%
		arrange(FDR) %>% mutate(FDR_ix = row_number()) %>%
		inner_join(operon_details) %>%
		mutate(FDR = case_when(FDR != 0 ~ FDR, FDR == 0 ~ min(FDR[FDR != 0]))) %>%
		{
			if (use_condensed_name) {
				mutate(., `Transcription Unit` = condensed_name)
			} else {
				mutate(., `Transcription Unit` = long_name)
			}
		} %>%
		ggplot(
			aes(x = operon_mLFC,
					y = FDR,
					fill = Pathways)) +
		geom_point(data = . %>% filter(Significance == "Significant"), aes(size = `Essentials in operon`, fill = Pathways, color = Significance), alpha = .65, stroke = 0.5, shape = 21) +
		geom_point(data = . %>% filter(Significance == "Not Significant"), aes(size = `Essentials in operon`, fill = Pathways, color = Significance), alpha = 0.35, stroke = 0, shape = 21) +
		
		geom_hline(yintercept = 0.05,
							 linetype = "dashed",
							 color = "#5A5A5A",
							 lwd = 0.5) +
		geom_vline(xintercept = -1,
							 linetype = "dashed",
							 color = "#5A5A5A", 
							 lwd = 0.5) +
		geom_vline(xintercept = 1,
							 linetype = "dashed",
							 color = "#5A5A5A", 
							 lwd = 0.5) +
		geom_vline(xintercept = 0,
							 linetype = "solid",
							 color = "#5A5A5A",
							 lwd = 0.5) +
		doc_theme +
		scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
		geom_label_repel(
			fill = alpha(c("white"), 0.65),
			max.iter = 1000000000,
			data = . %>% filter(
				FDR < 0.05 & 
					(operon_mLFC_ix <= max_left & operon_mLFC < -0.5 | 
					 	(operon_mLFC_desc_ix <= max_right & operon_mLFC > 0.5 ) | 
					 	FDR_ix <= max_top)),
			aes(label = `Transcription Unit`),
			force = 7.5,
			segment.size = 0.15,
			min.segment.length = 0,
			box.padding = 2,
			point.padding = .25,
			size = 3,
			max.overlaps = Inf,
			colour = "black") +
		facet_wrap(~ condition, scales = "free_x") +
		scale_fill_manual(
			values = c(
				"Other" = "grey",
				"Ribosome" = "#E31A1C",
				"Ribosome+Other" = "#E31A1C",
				"Ox Phos" = "#6A3D9A",
				"Ox Phos+Other" = "#6A3D9A",
				"LOS" = "#33A02C",
				"LOS+Other" = "#33A02C",
				"Cell Wall/PG" = "#FF7F00",
				"Cell Wall/PG+Other" = "#FF7F00",
				"tRNA Ligase" = "#1F78B4",
				"tRNA Ligase+Other" = "#1F78B4")) +
		# scale_fill_manual(
		# 	values = c(
		# 		"Other" = "grey",
		# 		"Ribosome" = "#E31A1C",
		# 		"Ribosome+Other" = "#FB9A99",
		# 		"Ox Phos" = "#6A3D9A",
		# 		"Ox Phos+Other" = "#CAB2D6",
		# 		"LOS" = "#33A02C",
		# 		"LOS+Other" = "#B2DF8A",
		# 		"Cell Wall/PG" = "#FF7F00",
		# 		"Cell Wall/PG+Other" = "#FDBF6F",
		# 		"tRNA Ligase" = "#1F78B4",
		# 		"tRNA Ligase+Other" = "#A6CEE3")) +
		guides(fill = guide_legend(
			override.aes = list(shape = 21, size = 5))) + 
		guides(alpha = guide_legend(
			override.aes = list(shape = 21, size = 5, fill = "black"))) +
		scale_alpha_manual(
			values = c(
				"Significant" = 0.65,
				"Not Significant" = 0.35)) +
		scale_color_manual(
			values = c(
				"Significant" = "black",
				"Not Significant" = alpha("white", 0))) +
		scale_linewidth_manual(
			values = c(
				"Significant" = 1,
				"Not Significant" = 0)) +
		scale_size_area(breaks = c(1, 5, 10, 13))
}

##########################################################################################
plot_operon(
	operon_median_results, 
	c("None_0_T1 - None_0_T0"), 
	use_condensed_name = TRUE, 
	max_left = 5, max_right = 5, max_top = 10, min_sig = 1e-150)

plot_operon(
	operon_median_results, 
	c("None_0_T2 - None_0_T0"), 
	use_condensed_name = TRUE, 
	max_left = 5, max_right = 5, max_top = 10, min_sig = 1e-150)

plot_operon(operon_median_results, c("None_0_T2 - None_0_T0"), 15, 5, 10)

plot_operon(operon_median_results, c("None_0_T2 - None_0_T1"), 15, 5, 10)


plot_operon(operon_median_results, c("Colistin_0.44_T1 - None_0_T1", "Colistin_0.44_T2 - None_0_T2"))

plot_operon(operon_median_results, c("Rifampicin_0.34_T2 - None_0_T2", "Colistin_0.44_T2 - None_0_T2"), TRUE, 7, 7, 5, 1e-100)

plot_operon(operon_median_results, c("Rifampicin_0.34_T2 - Colistin_0.44_T2"), FALSE, 10, 10, 10, 1e-150)

plot_operon(operon_median_results, c("Meropenem_0.17_T2 - None_0_T2", "Imipenem_0.09_T2 - None_0_T2"), FALSE, 10, 10, 10, 1e-150)

plot_operon(operon_median_results, c("Imipenem_0.09_T1 - None_0_T2"), FALSE, 10, 10, 5, 1e-150)

plot_operon(operon_median_results, c("Rifampicin_0.34_T1 - None_0_T1", "Rifampicin_0.34_T2 - None_0_T2"))

plot_operon(operon_median_results, c("Meropenem_0.17_T1 - None_0_T1", "Meropenem_0.17_T2 - None_0_T2"))

plot_operon(operon_median_results, c("Meropenem_0.17_T1 - None_0_T1", "Meropenem_0.17_T2 - None_0_T2"))

plot_operon(operon_median_results, c("Imipenem_0.09_T1 - None_0_T1", "Imipenem_0.09_T2 - None_0_T2"))

plot_operon(operon_median_results, c("Imipenem_0.09_T2 - Imipenem_0.06_T2"))

plot_operon(operon_median_results, c("Meropenem_0.17_T2 - Meropenem_0.11_T2"))

plot_operon(operon_median_results, c("Colistin_0.44_T2 - Colistin_0.44_T1"))

plot_operon(operon_median_results, c("Rifampicin_0.34_T2 - Rifampicin_0.34_T1"))

