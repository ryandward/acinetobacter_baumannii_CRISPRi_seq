library(pacman)

# source("meta_counter.R")

p_load(ggplot2, data.table, pheatmap, ggrepel, hrbrthemes, viridis)

doc_theme <- theme_ipsum(
	base_family = "Arial", 
	caption_margin = 12,
	axis_title_size = 12,
	axis_col = "black")

melted_results <- fread("Results/melted_results.tsv.gz", sep = "\t")

interest <- fread("interest.tsv", sep = "\t")

melted_results_of_interest <- copy(melted_results[condition %in% interest$condition])

melted_results_of_interest[, condition := gsub(" - ", "_vs_", condition)]

melted_results_of_interest[, FC.adj := 2^LFC.adj]

chemicals <-
	data.table(condition = tstrsplit(contrast_levels, " - ")[[1]])

chemicals <- unique(chemicals[!condition %like% "None"])

for (exp_shift in chemicals$condition) {
	
	
	exp_baseline <- melted_results_of_interest[shift %in% exp_shift, unique(base)]
	
	comp_shift  <- exp_baseline
	
	comp_baseline <- melted_results_of_interest[shift %in% exp_baseline, unique(base)]
	
	significant_spacers <- melted_results_of_interest[shift == exp_shift, unique(spacer)]
	
	pathway_analytics <-
		dcast(
			melted_results_of_interest[
				spacer %in% significant_spacers &
					((shift == exp_shift & base == exp_baseline) | 
					 	(shift == comp_shift & base == comp_baseline))],
			spacer + Pathway + AB19606 + unique_name + type + y_pred ~ shift + base,
			value.var = "LFC.adj")
	

	pathway_analytics_all <- copy(pathway_analytics)
	
	pathway_analytics_all[, Pathway := "All Genes"]
	
	pathway_analytics <- rbind(pathway_analytics, pathway_analytics_all)
	
	pathway_analytics <- pathway_analytics[type == "mismatch"]
	

	
	pathway_analytics[, Pathway := factor(Pathway, levels = unique(Pathway))]
	
	pathway_analytics_points <-
		pathway_analytics[
			Pathway != "All Genes" | 
				(Pathway == "All Genes" & 
				 	!spacer %in% pathway_analytics[
				 		Pathway!= "All Genes", spacer]) ]
	
	pathway_analytics <- pathway_analytics %>%
		 mutate(Pathway = case_when(
			Pathway == "Ribosome" ~ "Ribosome",
			Pathway == "LOS" ~ "LOS",
			unique_name %like% "nuo" ~ "NDH-1",
			Pathway %like% "Cell Wall" ~ "PG/Division",
			Pathway %like% "tRNA" ~ "tRNA Ligase",
			Pathway %like% "All Genes" ~ "All Genes"))
	
	pathway_analytics <- pathway_analytics[
		Pathway %in% c(
			"All Genes",
			"tRNA Ligase",
			"NDH-1",
			"PG/Division",
			"Ribosome"
			)]
	
	hw_sp <-
		ggplot(
			pathway_analytics,
			aes_string(
				x = paste(comp_shift, comp_baseline, sep = "_"),
				# x = "y_pred",
				y = paste(exp_shift, exp_baseline, sep = "_"),
				fill = "Pathway",
				color = "Pathway")) +
		# Put "All Genes" in the background
		geom_point(data = pathway_analytics[Pathway == "All Genes"]) +
		geom_point(data = pathway_analytics[Pathway != "All Genes"])
	
	
	hw_sp <-
		hw_sp +
		geom_smooth(
			# Put "All Genes" in the background
			data = pathway_analytics[Pathway == "All Genes"],
			method = glm,
			se = TRUE,
			fullrange = TRUE,
			formula = 'y~x') +
		geom_smooth(
			data = pathway_analytics[Pathway != "All Genes"],
			method = glm,
			se = TRUE,
			fullrange = TRUE,
			formula = 'y~x') +
		theme(legend.position = "bottom") +
		# scale_fill_viridis(discrete = TRUE, alpha = 0.5, direction = -1) +
		# scale_color_viridis(discrete = TRUE, alpha = 0.7, direction = -1) +
		ggtitle(paste(exp_shift, "vs", comp_shift)) +
		doc_theme + 
		scale_fill_manual(
			values = c(
				"All Genes" = alpha("grey", 0.25),
				"Ribosome" = "#E31A1C",
				"NDH-1" = "#CAB2D6",
				"LOS" = "#B2DF8A",
				"PG/Division" = "#FDBF6F",
				"tRNA Ligase" = "#A6CEE3"
				)) +
		scale_colour_manual(
			values = c(
				"All Genes" = alpha("grey", 0.25),
				"Ribosome" = "#E31A1C",
				"NDH-1" = "#CAB2D6",
				"LOS" = "#B2DF8A",
				"PG/Division" = "#FDBF6F",
				"tRNA Ligase" = "#A6CEE3"
				)) 
		theme(legend.position = "none")
	
	print(hw_sp)
	
	# save as 750x500

	}
