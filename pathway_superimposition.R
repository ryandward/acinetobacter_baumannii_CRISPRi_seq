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

gene_pathway_mapping <- annotated_fit_predictions %>%
	select(Gene, Pathway) %>%
	distinct()

annotated_vuln.summary <- vuln.summary %>%
	rename(Gene = unique_name, Condition = condition) %>%
	inner_join(annotated_fit_predictions %>% select(Condition, Gene) %>% unique) %>%
	left_join(gene_pathway_mapping)

interaction_colors <- c(
	"T1.Ribosome" = "#FB9A99",
	"T2.Ribosome" = "#E31A1C",
	"T1.Ox Phos" = "#CAB2D6",
	"T2.Ox Phos" = "#6A3D9A",
	"T1.LOS" = "#B2DF8A",
	"T2.LOS" = "#33A02C",
	"T1.PG/Division" = "#FDBF6F",
	"T2.PG/Division" = "#FF7F00",
	"T1.tRNA Ligase" = "#A6CEE3",
	"T2.tRNA Ligase" = "#1F78B4"
)

pathway_colors <- c(
	"Ribosome" = "#E31A1C",
	"Ox Phos" = "#6A3D9A",
	"LOS" = "#33A02C",
	"PG/Division" = "#FF7F00",
	"tRNA Ligase" = "#1F78B4"
)

medians <- annotated_fit_predictions %>%
	filter(Condition %in% c("None_0_T1 - None_0_T0", "None_0_T2 - None_0_T0")) %>%
	mutate(Condition = recode(
		Condition,
		"None_0_T1 - None_0_T0" = "T1",
		"None_0_T2 - None_0_T0" = "T2")) %>%
	group_by(Pathway, y_pred) %>%
	summarize(median_fitted = median(.fitted, na.rm = TRUE))


annotated_fit_predictions %>%
	filter(Condition %in% c("None_0_T1 - None_0_T0", "None_0_T2 - None_0_T0")) %>%
	mutate(Condition = recode(
		Condition,
		"None_0_T1 - None_0_T0" = "T1",
		"None_0_T2 - None_0_T0" = "T2"),
		Combined = interaction(Condition, Pathway, sep = ".", lex.order = TRUE)) %>%
	ggplot(aes(x = y_pred, y = .fitted, group = AB19606, color = Combined)) +
	geom_line(size = 0.75, alpha = 0.75) +
		scale_color_manual(values = c(pathway_colors, interaction_colors)) +
	labs(x = "Predicted Knockdown", y = "Log Logistic Fit") +
	doc_theme +
	theme(legend.position = "none",
				axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
	geom_point(
		data = (
			annotated_vuln.summary %>%
				filter(Condition %in% c("None_0_T1 - None_0_T0", "None_0_T2 - None_0_T0")) %>%
				mutate(Condition = recode(
					Condition,
					"None_0_T1 - None_0_T0" = "T1",
					"None_0_T2 - None_0_T0" = "T2")) ), 
		aes(x = vuln.kd_50, y = vuln.est, group = Gene), colour = "black", size = 2, alpha = 1, shape = 20) +
	facet_grid(Pathway ~ Condition) 
	# stat_smooth(aes(group = Pathway), colour = alpha("black", 0.75), method = "loess", span = 0.25, linetype = "twodash", size = 1.75, se = FALSE)
	# stat_summary(aes(group = interaction(y_pred, Pathway), color = NULL), fun = median, geom = "point", size = 1.75, color = alpha("black", 0.65))

	


