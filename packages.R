# Load required packages
require('pacman')
require(conflicted)
require(progress)

# Load packages using pacman
p_load(
	broom,
	colourpicker,
	data.table,
	dplyr,
	drc,
	edgeR,
	ggplot2,
	ggrepel,
	gtools,
	Hmisc,
	lmtest,
	magrittr,
	modelr,
	poolr,
#	pracma,
	purrr,
	R.utils,
	RColorBrewer,
	Rtsne,
	scales,
	statmod,
	svglite,
	tidyr,
	tidyverse,
	vegan
)

# Load packages from GitHub
p_load_current_gh(
	"DoseResponse/drcData",
	"hrbrmstr/hrbrthemes",
	"ryandward/drc",
	"jokergoo/ComplexHeatmap"
)

# Resolve conflicts using the conflicted package
conflicted::conflicts_prefer(
	gtools::permute,
	dplyr::filter,
	dplyr::select,
	drc::gaussian
)

# Add definitions

doc_theme <- theme_ipsum(
	base_family = "Arial", 
	caption_margin = 12,
	axis_title_size = 12,
	axis_col = "black")

# Update file_names lists with the new DRC objects entry
message("Updating file_names lists...\n")
file_names_full <- list(
	vuln_summary = "hormetic_vulnerability_summary_full.tsv.gz",
	fit_predictions = "hormetic_fit_predictions_full.tsv.gz",
	fit_points = "hormetic_fit_points_full.tsv.gz",
	model_performance = "hormetic_performance_full.tsv.gz",
	model_parameters = "hormetic_parameters_full.tsv.gz",
	drc_fits = "drc_fits_full.RDS"
)

file_names_reduced <- list(
	vuln_summary = "hormetic_vulnerability_summary_reduced.tsv.gz",
	fit_predictions = "hormetic_fit_predictions_reduced.tsv.gz",
	fit_points = "hormetic_fit_points_reduced.tsv.gz",
	model_performance = "hormetic_performance_reduced.tsv.gz",
	model_parameters = "hormetic_parameters_reduced.tsv.gz",
	drc_fits = "drc_fits_reduced.RDS"
)
message("File_names lists updated.\n")

safe_fread <- purrr::possibly(.f = fread, otherwise = data.frame())

model_comparisons <- safe_fread("model_comparisons.tsv")
hormesis_results <- safe_fread("model_comparisons_hormesis.tsv")
aba_key <- safe_fread("aba_key.tsv")
interested.genes <- safe_fread("curated_names.tsv") %>% pull(unique_name)
curated_names <- safe_fread("curated_names.tsv")
melted_results <- safe_fread("Results/melted_results.tsv.gz")
median_melted_results <- safe_fread("Results/median_melted_results.tsv.gz")