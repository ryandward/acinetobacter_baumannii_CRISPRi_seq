# Load several packages from CRAN and Bioconductor
require('pacman')
p_load(
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
	magrittr
)

p_load(
	"data.table",
	"tidyverse",
	"broom",
	"modelr")

p_load_current_gh(
	"DoseResponse/drcData",
	"ryandward/drc",
	"hrbrmstr/hrbrthemes")

conflicted::conflicts_prefer(
	gtools::permute,
	dplyr::filter,
	dplyr::select,
	drc::gaussian)

doc_theme <- theme_ipsum(
	base_family = "Arial", 
	caption_margin = 12,
	axis_title_size = 12,
	axis_col = "black")


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

plot_gene_dose_effect <- function(results_data, unique_names, conditions, colors) {
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
	
	plot.graphic <- plot.fit_predictions %>% 
		filter(
			Gene %in% plot.genes & Condition %in% plot.conditions) %>%
		mutate(label = gsub("", "", Condition)) %>%
		ggplot() +
		geom_hline(
			yintercept = 0, 
			linetype = "solid", 
			color = "black", 
			linewidth = 0.75) +
		geom_vline(
			xintercept = 0, 
			linetype = "solid", 
			color = "black", 
			linewidth = 0.75) +
		geom_vline(
			data = plot.parameters,
			aes(xintercept = as.numeric(vuln.kd_50)),
			linetype = "dashed", 
			colour = "black",
			linewidth = 0.25) +
		geom_hline(
			data = plot.parameters,
			aes(yintercept = as.numeric(vuln.est)),
			colour = "black",
			linetype = "dashed", 
			linewidth = 0.25) +
		geom_line(
			alpha = 0.50, 
			size = 3, 
			aes(
				x = y_pred, 
				y = .fitted, 
				color = Gene)) +
		geom_point(
			data = plot.fit_points %>%
				filter(
					
					Gene %in% plot.genes &
						Condition %in% plot.conditions) %>%
				filter(
					Gene %in% plot.genes & 
						Condition %in% plot.conditions),
			shape = 20, 
			size = 3.5,
			aes(
				x = y_pred, 
				y = LFC.adj, 
				color = Gene)) + 
		# geom_ribbon(
		# 	data = plot.fit_predictions %>%
		# 		filter(
		# 			Gene %in% plot.genes &
		# 				Condition %in% plot.conditions) %>%
		# 		filter(
		# 			Gene %in% plot.genes &
		# 				Condition %in% plot.conditions),
		# 	alpha = 0.25,
		# aes(
		# 	x = y_pred,
		# 	y = .fitted,
		# 	ymin = .lower,
		# 	ymax = .upper,
		# 	fill = Gene)) +
		scale_fill_manual(values = colors) +
		scale_color_manual(values = colors) +
		ylab("Fitness (Log2)") +
		doc_theme +
		theme(legend.position = "none") +
	facet_wrap(~factor(Gene, levels = unique_names)) + 
		guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
		ggtitle(conditions)
	
	print(plot.graphic)
}

gene_colors <- c(
	"nuoB" = "#6A3D9A",
	"lpxC" = "#33A02C",
	"glnS" = "#1F78B4",
	"murA" = "#FF7F00",
	"rpmB" = "#E31A1C")

# beautiful 
plot_gene_dose_effect(
	reduced_results,
	c("rpmB", "murA", "nuoB", "lpxC"), 
	c("None_0_T1 - None_0_T0"), 
	gene_colors)

# not justifiable use, check model comparison results
# plot_gene_dose_effect(
# 	full_results,
# 	c("rpmB", "murA", "nuoB", "lpxC"), 
# 	c("None_0_T1 - None_0_T0"), 
# 	gene_colors)



#lrt p-value = 1.281345e-05
plot_gene_dose_effect(
	full_results,
	c("nuoB"), 
	c("Rifampicin_0.34_T2 - None_0_T0"), 
	gene_colors)

#lrt p-value = 1.281345e-05
plot_gene_dose_effect(
	reduced_results,
	c("nuoB"), 
	c("None_0_T2 - None_0_T0"), 
	gene_colors)



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
 	c("lysC", "glnS", "trpS", "tyrS"), 
 	c("Meropenem_0.17_T2 - None_0_T0"), 
 	gene_colors)
 