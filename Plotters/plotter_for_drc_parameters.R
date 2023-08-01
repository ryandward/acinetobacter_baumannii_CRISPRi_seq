source("drc_logistic_functions.R")
source("packages.R")

melted_results <- fread("Results/melted_results.tsv.gz", sep = "\t")
median_melted_results <- fread("Results/median_melted_results.tsv.gz", sep = "\t")

if (!exists("full_results")) {
	if (exists("file_names_full")) {
		full_results <- read_results(file_names_full)
	} else {
		message("file_names_full not found.")
	}
}

if (!exists("reduced_results")) {
	if (exists("file_names_reduced")) {
		reduced_results <- read_results(file_names_reduced)
	} else {
		message("file_names_reduced not found.")
	}
}

process_data <- function(data_list, data_name) {
	data <- data_list[[data_name]]
	data %>%
		mutate(
			Timing = case_when(
				Condition %like% "T1" ~ "T1",
				Condition %like% "T2" ~ "T2"),
			Drug = case_when(
				Condition %like% "Rifampicin" ~ "Rifampicin",
				Condition %like% "Colistin" ~ "Colistin",
				Condition %like% "Meropenem" ~ "Meropenem",
				Condition %like% "Imipenem" ~ "Imipenem",
				Condition %like% "^None_0" ~ "Inducer only"),
			Drug = factor(Drug, levels = c("Inducer only", "Colistin", "Rifampicin", "Meropenem", "Imipenem"))
		)
}

# plot.fit_predictions <- process_data(full_results, "fit_predictions")
# plot.fit_points <- process_data(full_results, "fit_points")
# plot.vulnerability <- full_results[['vuln.summary']]

plot_gene_dose_effect <- function(results_data, unique_names, conditions, colors, bands = FALSE, hormesis = NULL) {
	plot.fit_predictions <- process_data(results_data, "fit_predictions")
	plot.fit_points <- process_data(results_data, "fit_points")
	plot.vulnerability <- results_data[['vuln.summary']]
	
	plot.genes <- unique_names
	plot.conditions <- conditions
	
	plot.fit_predictions <- plot.fit_predictions %>%
		filter(Gene %in% unique_names) %>%
		mutate(Gene = factor(Gene, levels = unique(unique_names)))
	
	plot.fit_points <- plot.fit_points %>% 
		filter(Gene %in% unique_names) %>%
		mutate(Gene = factor(Gene, levels = unique(unique_names)))
	
	plot.parameters <- plot.vulnerability %>%
		filter(unique_name %in% plot.genes) %>% 
		filter(condition %in% plot.conditions) %>%
		filter(vuln.p <= 0.05) %>%
		rename(Condition = condition, Gene = unique_name)
	
	if (!is.null(hormesis) && nrow(hormesis) > 0) {
		hormesis_lines <- hormesis %>%
			filter(unique_name %in% unique_names & condition %in% conditions) %>%
			mutate(
				condition = factor(condition, levels = conditions),
				max_response_fitted = ifelse(highest_fitted > response.max, highest_fitted, lowest_fitted),
				max_dose_pred = ifelse(highest_fitted > response.max, highest_y_pred, lowest_y_pred)
			)
	} else {
		hormesis_lines <- NULL
	}
	
	plot.graphic <- plot.fit_predictions %>%
		filter(Gene %in% plot.genes & Condition %in% plot.conditions) %>%
		mutate(label = gsub("", "", Condition)) %>%
		ggplot() +
		geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.75) +
		geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.75) +
		geom_vline(data = plot.parameters, aes(xintercept = as.numeric(vuln.kd_50)), linetype = "dashed", colour = "black", linewidth = 0.25) +
		geom_hline(data = plot.parameters, aes(yintercept = as.numeric(vuln.est)), colour = "black", linetype = "dashed", linewidth = 0.25) +
		geom_point(data = plot.parameters, aes(x = as.numeric(vuln.kd_50), y = as.numeric(vuln.est), color = Gene), shape = 1, size = 6, stroke = 1.5) +
		geom_line(alpha = 0.50, size = 3, aes(x = y_pred, y = .fitted, color = Gene)) +
		geom_point(data = plot.fit_points %>% filter(Gene %in% plot.genes & Condition %in% plot.conditions), shape = 20, size = 3.5, aes(x = y_pred, y = LFC.adj, color = Gene))
	
	if (bands) {
		plot.graphic <- plot.graphic + geom_ribbon(data = plot.fit_predictions %>% filter(Gene %in% plot.genes & Condition %in% plot.conditions), alpha = 0.25, aes(x = y_pred, y = .fitted, ymin = .lower, ymax = .upper, fill = Gene))
	}
	
	if (!is.null(hormesis_lines)) {
		plot.graphic <- plot.graphic +
			geom_hline(data = hormesis_lines, aes(yintercept = max_response_fitted), color = "red", linetype = "solid", size = 0.75, alpha = 0.5)
		#geom_vline(data = hormesis_lines, aes(xintercept = max_dose_pred), color = "blue", linetype = "dashed", size = 0.25)
	}
	
	plot.graphic <- plot.graphic +
		scale_fill_manual(values = colors) +
		scale_color_manual(values = colors) +
		ylab("Fitness (Log2)") +
		doc_theme +
		theme(legend.position = "none") +
		facet_grid(cols = vars(Condition), rows = vars(Gene)) +
		guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
		ggtitle(paste(conditions, collapse = ", "))
	
	print(plot.graphic)
}

							 


gene_colors <- c(
	"nuoB" = "#6A3D9A",
	"lpxC" = "#33A02C",
	"lpxA" = "#33A02C",
	"glnS" = "#1F78B4",
	"murA" = "#FF7F00",
	"rpmB" = "#E31A1C")

# beautiful 
plot_gene_dose_effect(
	reduced_results,
	c("glnS", "murA", "nuoB", "lpxA", "aroC", "rpmB", "GO593_00515"), 
	c("None_0_T1 - None_0_T0", "None_0_T2 - None_0_T0"), 
	gene_colors)

# beautiful 
plot_gene_dose_effect(
	reduced_results,
	c("accD"), 
	c("Colistin_0.44_T2 - None_0_T0"), 
	gene_colors)

plot_gene_dose_effect(
	reduced_results,
	c("glnS", "murA", "nuoB", "lpxA"), 
	c("None_0_T2 - None_0_T0"), 
	gene_colors)

# beautiful 
plot_gene_dose_effect(
	reduced_results,
	c("glnS", "murA", "nuoB", "lpxA"), 
	c("None_0_T1 - None_0_T0"), 
	gene_colors)

plot_gene_dose_effect(
	reduced_results,
	c("murA", "nuoB", "lpxA"), 
	c("None_0_T1 - None_0_T0", "None_0_T2 - None_0_T0"), 
	gene_colors)


# not justifiable use, check model comparison results
# plot_gene_dose_effect(
# 	full_results,
# 	c("rpmB", "murA", "nuoB", "lpxC"), 
# 	c("None_0_T1 - None_0_T0"), 
# 	gene_colors)


plot_gene_dose_effect(
	reduced_results,
	consistent_genes %>% inner_join(curated_names_operons_pathways) %>% filter(Pathways %like% "tRNA") %>% pull(unique_name), 
	c("None_0_T2 - None_0_T0"), 
	gene_colors,
	bands = FALSE)

plot_gene_dose_effect(
	reduced_results,
	consistent_genes %>% inner_join(curated_names_operons_pathways) %>% filter(Pathways %like% "tRNA") %>% pull(unique_name), 
	c("Imipenem_0.09_T2 - None_0_T0"), 
	gene_colors,
	bands = FALSE)

plot_gene_dose_effect(
	full_results,
	c("nuoB"), 
	c("Rifampicin_0.34_T2 - None_0_T0"), 
	gene_colors,
	bands = FALSE,
	hormesis = hormesis_results)

plot_gene_dose_effect(
	full_results,
	c("glnS"), 
	c("Imipenem_0.09_T2 - None_0_T0"), 
	gene_colors,
	bands = FALSE,
	hormesis = hormesis_results)

plot_gene_dose_effect(
	reduced_results,
	c("glnS"), 
	c("Imipenem_0.09_T2 - None_0_T0"), 
	gene_colors,
	bands = FALSE)


plot_gene_dose_effect(
	full_results,
	c("sdhD"), 
	c("Rifampicin_0.34_T2 - None_0_T0"), 
	gene_colors,
	bands = FALSE,
	hormesis = hormesis_results)


plot_gene_dose_effect(
	full_results,
	c("rplR"), 
	c("Imipenem_0.09_T2 - None_0_T0"), 
	gene_colors,
	bands = FALSE)


# beautiful
plot_gene_dose_effect(
	full_results,
	c("glnS"), 
	c("Imipenem_0.06_T2 - None_0_T0"), 
	gene_colors)


 plot_gene_dose_effect(
 	reduced_results,
 	c("lysC", "glnS", "trpS", "tyrS"), 
	c("Imipenem_0.09_T2 - None_0_T0"), 
 	gene_colors)
 
 plot_gene_dose_effect(
 	full_results,
 	c("lysC", "glnS", "trpS", "tyrS"), 
 	c("Imipenem_0.09_T2 - None_0_T0"), 
 	gene_colors)
 
 plot_gene_dose_effect(
 	reduced_results,
 	c("lysC", "glnS", "trpS", "tyrS"), 
 	c("Meropenem_0.17_T2 - None_0_T0"), 
 	gene_colors)

 plot_gene_dose_effect(
 	full_results,
 	c("glnS", "lptB", "ispB", "rho"), 
 	c("Imipenem_0.09_T2 - None_0_T0"), 
 	gene_colors)
 
 plot_gene_dose_effect(
 	reduced_results,
 	c("rpmF", "glnS", "trpS", "acpT"), 
 	c("None_0_T2 - None_0_T0"), 
 	gene_colors)
 
 plot_gene_dose_effect(
 	full_results,
 	c("GO593_13220", "glnS", "lysC"), 
 	c("Meropenem_0.11_T2 - None_0_T0"), 
 	gene_colors)
 
 
 plot_gene_dose_effect(
 	full_results,
 	c("glnS", "valS"), 
 	c("Imipenem_0.09_T2 - None_0_T0"), 
 	gene_colors,
 	bands = FALSE,
 	hormesis = hormesis_results)
 
 plot_gene_dose_effect(
 	full_results,
 	c("glnS", "leuS"), 
 	c("Imipenem_0.09_T2 - None_0_T0"), 
 	gene_colors,
 	bands = FALSE,
 	hormesis = hormesis_results)
 
 plot_gene_dose_effect(
 	full_results,
 	c("topA"), 
 	c("None_0_T2 - None_0_T0"), 
 	gene_colors,
 	bands = FALSE,
 	hormesis = hormesis_results)
 
 plot_gene_dose_effect(
 	full_results,
 	c("nuoB"), 
 	c("Rifampicin_0.34_T2 - None_0_T0"), 
 	gene_colors,
 	bands = FALSE,
 	hormesis = hormesis_results)
 
 plot_gene_dose_effect(
 	reduced_results,
 	c("nuoB"), 
 	c("Rifampicin_0.34_T2 - None_0_T0"), 
 	gene_colors,
 	bands = FALSE)
 