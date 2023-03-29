fit_predictions <- fread("Results/hormetic_fit_predictions.tsv.gz")
fit_points <- fread("Results/hormetic_fit_points.tsv.gz")
vulnerability <- fread("Results/hormetic_vulnerability_summary.tsv.gz")
model_performance <- fread("Results/hormetic_performance.tsv.gz")

process_data <- function(data) {
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

plot.fit_predictions <- process_data(fit_predictions)
plot.fit_points <- process_data(fit_points)


plot_gene_dose_effect <- function(unique_names, conditions, colors) {
	
	plot.genes <- unique_names
	plot.conditions <- conditions
	
	plot.fit_predictions <- plot.fit_predictions %>%
		filter(Gene %in% unique_names) %>%
		mutate(Gene = factor(Gene, levels = unique(unique_names)))
	
	plot.fit_points <- plot.fit_points %>% 
		filter(Gene %in% unique_names) %>%
		mutate(Gene = factor(Gene, levels = unique(unique_names)))
	
	plot.parameters <- vulnerability %>%
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
		guides(fill = guide_legend(nrow = 1, byrow = TRUE))
	
	print(plot.graphic)
}

gene_colors <- c(
	"nuoB" = "#6A3D9A",
	"lpxC" = "#33A02C",
	"glnS" = "#1F78B4",
	"murA" = "#FF7F00",
	"rpmB" = "#E31A1C")

plot_gene_dose_effect(
	c("rpmB", "murA", "nuoB", "lpxC"), 
	c("None_0_T1 - None_0_T0"), 
	gene_colors)

plot_gene_dose_effect(
	c("nuoB"), 
	c("Rifampicin_0.34_T2 - None_0_T0"), 
	gene_colors)

plot_gene_dose_effect(
	c("glnS"), 
	c("Imipenem_0.09_T2 - None_0_T0"), 
	gene_colors)
