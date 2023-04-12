require(conflicted)
require(pacman)

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

fit_predictions <- fread("Results/hormetic_fit_predictions.tsv.gz")
fit_points <- fread("Results/hormetic_fit_points.tsv.gz")

plot.genes <- c("lpxC","nuoB", "glnS", "murA")
# plot.genes <- c("lpxC", "nuoB",)
# plot.genes <- c("lpxC", "murA", "rpmB", "aroC", "GO593_00515")

doc_theme <- theme_ipsum(
	base_family = "Arial", 
	caption_margin = 12,
	axis_title_size = 12,
	axis_col = "black")

plot.conditions <- c(
	"None_0_T1 - None_0_T0",
	"None_0_T2 - None_0_T0",
	"Rifampicin_0.34_T1 - None_0_T0",
	"Rifampicin_0.34_T2 - None_0_T0",
	"Colistin_0.44_T1 - None_0_T0",
	"Colistin_0.44_T2 - None_0_T0",
	"Imipenem_0.09_T1 - None_0_T0",
	"Imipenem_0.09_T2 - None_0_T0")

# plot.conditions <- c(
# 	"None_0_T1 - None_0_T0",
# 	"None_0_T2 - None_0_T0",
# 	"Rifampicin_0.34_T1 - None_0_T0",
# 	"Rifampicin_0.34_T2 - None_0_T0",
# 	"Colistin_0.44_T1 - None_0_T0",
# 	"Colistin_0.44_T2 - None_0_T0")

# plot.conditions <- c(
# 	"None_0_T1 - None_0_T0",
# 	"None_0_T2 - None_0_T0")


plot.fit_predictions <-
	fit_predictions %>% 
	mutate(
		Timing = case_when(
			Condition %like% "T1" ~ "T1",
			Condition %like% "T2" ~ "T2"),
		Drug = case_when(
			Condition %like% "Rifampicin" ~ "Rifampicin",
			Condition %like% "Colistin" ~ "Colistin",
			Condition %like% "Meropenem" ~ "Meropenem",
			Condition %like% "Imipenem" ~ "Imipenem",
			Condition %like% "^None_0" ~ "No drug"),
		Drug = factor(Drug, levels = c("No drug", "Colistin", "Rifampicin", "Meropenem", "Imipenem")))

plot.fit_points <-
	fit_points %>% 
	mutate(
		Timing = case_when(
			Condition %like% "T1" ~ "T1",
			Condition %like% "T2" ~ "T2"),
		Drug = case_when(
			Condition %like% "Rifampicin" ~ "Rifampicin",
			Condition %like% "Colistin" ~ "Colistin",
			Condition %like% "Meropenem" ~ "Meropenem",
			Condition %like% "Imipenem" ~ "Imipenem",
			Condition %like% "^None_0" ~ "No drug"),
		Drug = factor(Drug, levels = c("No drug", "Colistin", "Rifampicin", "Meropenem", "Imipenem")))

# plot.labeller <- as_labeller(
# 	c(
# 		`None_0_T1 - None_0_T0` = "Induction Only (T1)",
# 		`None_0_T2 - None_0_T0` = "Induction Only (T2)",
# 		`Rifampicin_0.34_T1 - None_0_T0` = "Rifampicin (T1)",
# 		`Rifampicin_0.34_T2 - None_0_T0` = "Rifampicin (T2)",
# 		`Colistin_0.44_T1 - None_0_T0` = "Colistin (T1)",
# 		`Colistin_0.44_T2 - None_0_T0` = "Colistin (T2)"))

# plot.title <- bquote(bold("Gene dose effect on drug activity" ~ at ~ hour[36] ~ '(Confidence = 0.90)'))

plot.graphic <- plot.fit_predictions %>% 
	filter(
		Gene %in% plot.genes &
			Condition %in% plot.conditions) %>%
	mutate(label = gsub("", "", Condition)) %>%
	ggplot() +
	# geom_rect(
	# 	data = . %>% select(Drug) %>% unique,
	# 	aes(fill = Drug),
	# 	xmin = -Inf,
	# 	xmax = Inf,
	# 	ymin = -Inf,
	# 	ymax = Inf,
	# 	alpha = 0.15) +
	geom_hline(
		yintercept = 0, 
		linetype = "dashed", 
		color = "black", 
		linewidth = 0.5) +
	geom_line(
		alpha = 1, 
		linewidth = 2, 
		aes(
			x = y_pred, 
			y = .fitted, 
			color = interaction(Timing, Gene))) +
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
			color = interaction(Timing, Gene))) + 
	# geom_ribbon(
	# 	data = plot.fit_predictions %>%
	# 		filter(
	# 			Gene %in% plot.genes &
	# 				Condition %in% plot.conditions) %>%
	# 		filter(
	# 			Gene %in% plot.genes &
	# 				Condition %in% plot.conditions),
	# 	alpha = 0.25,
	# 	aes(
	# 		x = y_pred,
	# 		y = .fitted,
	# 		ymin = .lower,
	# 		ymax = .upper,
	# 		fill = interaction(Timing, Gene))) +
	scale_fill_manual(
		values = c(
			"T1.nuoB" = "#CAB2D6",
			"T2.nuoB" = "#6A3D9A",
			"T1.lpxC" = "#B2DF8A",
			"T2.lpxC" = "#33A02C",
			"T1.glnS" = "#A6CEE3",
			"T2.glnS" = "#1F78B4",
			"T1.murA" = "#FDBF6F",
			"T2.murA" = "#FF7F00",
			"T1.aroC" = "light grey",
			"T2.aroC" = "dark grey",
			"T1.rpmB" = "light grey",
			"T2.rpmB" = "dark grey",
			"T1.GO593_00515" = "light grey",
			"T2.GO593_00515" = "dark grey")) +
	scale_color_manual(
		values = c(
			"T1.nuoB" = "#CAB2D6",
			"T2.nuoB" = "#6A3D9A",
			"T1.lpxC" = "#B2DF8A",
			"T2.lpxC" = "#33A02C",
			"T1.glnS" = "#A6CEE3",
			"T2.glnS" = "#1F78B4",
			"T1.murA" = "#FDBF6F",
			"T2.murA" = "#FF7F00",
			"T1.aroC" = "light grey",
			"T2.aroC" = "dark grey",
			"T1.rpmB" = "light grey",
			"T2.rpmB" = "dark grey",
			"T1.GO593_00515" = "light grey",
			"T2.GO593_00515" = "dark grey")) +
	xlab("Knockdown") +
	ylab("Fitness (Log2)") +
	doc_theme +
	theme(legend.position = "bottom") +
	facet_grid(scales = "free_y",
		facets = c("Gene", "Drug")) + 
	guides(fill = guide_legend(nrow = 2, byrow = TRUE))

print(plot.graphic)

