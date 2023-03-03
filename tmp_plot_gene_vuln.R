plot.genes <- c("lpxC","nuoB")

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
	"Colistin_0.44_T2 - None_0_T0")


plot.fit_predictions <-
	fit_predictions %>% 
	mutate(
		Timing = case_when(
			Condition %like% "T1" ~ "T1",
			Condition %like% "T2" ~ "T2"),
		Drug = case_when(
			Condition %like% "^None_0" ~ "No drug",
			Condition %like% "Rifampicin" ~ "Rifampicin",
			Condition %like% "Colistin" ~ "Colistin"),
		Drug = factor(Drug, levels = c("No drug", "Colistin", "Rifampicin")))

plot.fit_points <-
	fit_points %>% 
	mutate(
		Timing = case_when(
			Condition %like% "T1" ~ "T1",
			Condition %like% "T2" ~ "T2"),
		Drug = case_when(
			Condition %like% "^None_0" ~ "No drug",
			Condition %like% "Rifampicin" ~ "Rifampicin",
			Condition %like% "Colistin" ~ "Colistin"),
		Drug = factor(Drug, levels = c("No drug", "Colistin", "Rifampicin")))

plot.labeller <- as_labeller(
	c(
		`None_0_T1 - None_0_T0` = "Induction Only (T1)",
		`None_0_T2 - None_0_T0` = "Induction Only (T2)",
		`Rifampicin_0.34_T1 - None_0_T0` = "Rifampicin (T1)",
		`Rifampicin_0.34_T2 - None_0_T0` = "Rifampicin (T2)",
		`Colistin_0.44_T1 - None_0_T0` = "Colistin (T1)",
		`Colistin_0.44_T2 - None_0_T0` = "Colistin (T2)"))

plot.title <- bquote(bold("Gene dose effect on drug activity" ~ at ~ hour[36] ~ '(Confidence = 0.90)'))

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
		size = 0.5) +
	geom_line(
		alpha = 1, 
		size = 2, 
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
			"No drug" = "white",
			"Colistin" = "#A6CEE3",
			"Rifampicin" = "#FB9A99",
			"T1.nuoB" = "#CAB2D6",
			"T2.nuoB" = "#6A3D9A",
			"T1.lpxC" = "#B2DF8A",
			"T2.lpxC" = "#33A02C")) +
	scale_color_manual(
		values = c(
			"T1.nuoB" = "#CAB2D6",
			"T2.nuoB" = "#6A3D9A",
			"T1.lpxC" = "#B2DF8A",
			"T2.lpxC" = "#33A02C")) +
	xlab("Knockdown") +
	ylab("Fitness (Log2)") +
	doc_theme +
	theme(legend.position = "bottom") +
	facet_grid(
		facets = c("Gene", "Drug")) + 
	guides(fill = guide_legend(nrow = 2, byrow = TRUE))

print(plot.graphic)

