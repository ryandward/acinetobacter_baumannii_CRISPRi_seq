# Load necessary packages
source("packages.R")
source("drc_logistic_functions.R")

curated_names_operons_pathways <- fread("curated_names_operons_pathways.tsv")

reduced_results <- tryCatch(
	read_results(file_names_reduced, output_dir = "Results", exclude_drc_fits = TRUE),
	error = function(e) {
		message("Failed to load reduced results: ", e$message)
		return(NULL)
	}
)

full_results <- tryCatch(
	read_results(file_names_full, output_dir = "Results", exclude_drc_fits = TRUE),
	error = function(e) {
		message("Failed to load full results: ", e$message)
		return(NULL)
	}
)



melted_results <- fread("Results/melted_results.tsv.gz", sep = "\t")
median_melted_results <- fread("Results/median_melted_results.tsv.gz", sep = "\t")
interest <- fread("interest.tsv", sep = "\t")
curated_names <- fread("curated_names.tsv", sep = "\t")
curated_names_operons_pathways <- fread("curated_names_operons_pathways.tsv")
consistent_genes <- fread("consistent_genes.tsv")


vuln.summary <- reduced_results$vuln.summary
fit_predictions <- reduced_results$fit_predictions
fit_points <- reduced_results$fit_points

annotated_fit_predictions <- 
	fit_predictions %>% 
	filter(Gene %in% consistent_genes$unique_name) %>%
	inner_join(curated_names_operons_pathways %>% rename(Gene = unique_name)) %>%
	rename(Pathway = Pathways) %>%
	mutate(Pathway = case_when(
		Pathway %like% "Ribosome" ~ "Ribosome",
		Pathway %like% "LOS" ~ "LOS",
		Pathway %like% "tRNA" ~ "tRNA Ligase",
		Pathway %like% "PG" ~ "PG/Division",
		Pathway %like% "Ox Phos" ~ "Ox Phos")) %>%
	mutate(Pathway = factor(Pathway, levels = c("Ribosome", "PG/Division", "Ox Phos", "tRNA Ligase", "LOS"))) %>%
	filter(!is.na(Pathway))

rough_fits <- annotated_fit_predictions %>% 
	filter(Condition %in% c(
		"None_0_T1 - None_0_T0")) %>%
	mutate(binned_y_pred = round(y_pred, 6)) %>% 
	group_by(Condition, Pathway, binned_y_pred) %>%
	summarise(binned_fit = mean(.fitted)) 

rough_fits <- rough_fits %>%
	nest(data = c(-Condition, -Pathway))



rough_fits <- rough_fits %>% 
	mutate(
		fit = map(
			data, 
			~ drm.try(
				data = .x, 
				binned_fit ~ binned_y_pred, 
				fct = BC.5(names = c("hill", "min_value", "max_value", "kd_50", "hormesis")))))

rough_fits <- rough_fits %>% 
	mutate(results = map(fit, glance)) %>%
	mutate(p.vals = map(fit, tidy)) %>%
	mutate(results = map(results, ~mutate(.x, logLik = c(logLik)))) %>%
	unnest(results)

rough_fits <- rough_fits %>%
	mutate(
		kd_50.tibble = map(
			p.vals, 
			~filter(.x, term == 'kd_50') %>%
				select(p.value) %>%
				rename (., vuln.p = p.value))) %>%
	mutate(
		vuln.tibble = map2(
			fit, 
			p.vals, 
			~augment.try(
				.x,
				newdata = .y %>% filter (term == "kd_50") %>% select (estimate)) %>% 
				rename (., vuln.est = .fitted, vuln.kd_50 = estimate))) %>%
	mutate(
		hill.tibble = map(
			p.vals, 
			~filter(.x, term == 'hill') %>%
				select(estimate, p.value) %>%
				rename (., hill.est = estimate, hill.p = p.value)))

rough_fits <- rough_fits %>%
	unnest(vuln.tibble) %>%
	unnest(kd_50.tibble) %>%
	unnest(hill.tibble) %>%
	select(-c(p.vals))

rough_fits <- rough_fits %>% 
	mutate(predictions = map2(fit, data, ~augment.try(
		.x,
		newdata = expand.grid(
			y_pred = seq(
				min(.y$binned_y_pred), 
				max(.y$binned_y_pred), 
				length = 250)),
		conf.int = T,
		conf.level = 0.90)))

rough_fit_predictions <- rough_fits %>% 
	select(Pathway, Condition, predictions) %>% 
	unnest(predictions)

plot.graphic <- annotated_fit_predictions %>% 
	filter(Condition %in% c(
		"None_0_T1 - None_0_T0")) %>%
	filter(Pathway != "Other") %>%
	filter(!is.na(Pathway)) %>%
	ggplot() +
	geom_hline(
		yintercept = 0, 
		linetype = "dashed", 
		color = "black", 
		size = 0.5) +
	geom_line(
		size = 1.5, 
		aes(
			x = y_pred, 
			y = .fitted, 
			group = Gene,
			colour = Pathway),
		alpha = 0.5) +
	geom_line(
		data = rough_fit_predictions,
		aes(x = y_pred, y = .fitted),
		colour = "black",
		lwd = 3
	) +
	xlab("Predicted Knockdown") +
	ylab("Log-2 Fitness Foldchange") +
	doc_theme +
	theme(legend.position = "none") +
	facet_wrap(facets = c("Pathway"), nrow = 2) +
	scale_colour_manual(
		values = c(
			"Ribosome" = "#FB9A99", 
			"Ox Phos" = "#CAB2D6",
			"LOS" = "#B2DF8A",
			"PG/Division" = "#FDBF6F"))


print(plot.graphic)
