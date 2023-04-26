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

pathway_colors <- c(
	"Ox Phos" = "#6A3D9A",
	"LOS" = "#33A02C",
	"tRNA Ligase" = "#1F78B4",
	"PG/Division" = "#FF7F00",
	"Ribosome" = "#E31A1C")

annotated_fit_predictions %>%
	filter(Condition %in% c("None_0_T1 - None_0_T0", "None_0_T2 - None_0_T0")) %>%
	filter(Pathway %in% c("Ox Phos", "LOS")) %>%
	ggplot(aes(x = y_pred, y = .fitted, group = AB19606, color = Pathway)) +
	# Plot individual fitted lines for each gene (locus tag) within a pathway
	geom_line(linewidth = 1.0, alpha = 0.25) +
	# Calculate and plot the smoothed central tendency line for each pathway on the fly
	stat_smooth(aes(group = Pathway), method = "loess", linetype = "dotdash", linewidth = 1.5, se = FALSE) +
	# Apply the specified colors to the pathways
	scale_color_manual(values = pathway_colors) +
	labs(x = "y_pred", y = "Fitted values") +
	doc_theme +
	facet_grid(Condition ~ Pathway, scales = "free_y")

