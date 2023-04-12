# Load necessary packages
require(conflicted)
require(pacman)
source("drc_logistic_functions.R")

# Load data
aba_key <- fread("aba_key.tsv")
interested.genes <- fread("curated_names.tsv") %>% pull(unique_name)
interested.conditions <- c( # List of conditions
	"None_0_T1 - None_0_T0",
	"None_0_T2 - None_0_T0",
	"Colistin_0.44_T1 - None_0_T0",
	"Colistin_0.44_T2 - None_0_T0",
	"Rifampicin_0.34_T1 - None_0_T0",
	"Rifampicin_0.34_T2 - None_0_T0",
	"Meropenem_0.11_T1 - None_0_T0",
	"Meropenem_0.17_T1 - None_0_T0",
	"Meropenem_0.11_T2 - None_0_T0",
	"Meropenem_0.17_T2 - None_0_T0",
	"Imipenem_0.06_T1 - None_0_T0",
	"Imipenem_0.09_T1 - None_0_T0",
	"Imipenem_0.06_T2 - None_0_T0",
	"Imipenem_0.09_T2 - None_0_T0") 

# Read results
melted_results <- fread("Results/melted_results.tsv.gz", sep = "\t")
median_melted_results <- fread("Results/median_melted_results.tsv.gz", sep = "\t")
curated_names <- fread("curated_names.tsv", sep = "\t")
aba_key <- fread("aba_key.tsv", sep = "\t")

# Pre-process data
melted_results$y_pred <- NULL
melted_results[aba_key, on = .(spacer), y_pred := y_pred]

# Define error handling functions
drm.try <- possibly(drm, otherwise = NA)
augment.try <- possibly(augment, otherwise = NA)
glance.try <- possibly(glance, otherwise = NA)
tidy.try <- possibly(tidy, otherwise = NA)

# Estimate dose-response models
mismatches <- melted_results %>%
	filter(y_pred > 0) %>%
	# mutate(y_pred = case_when(y_pred < 0 ~ 0, TRUE ~ y_pred)) %>%
	filter(condition %in% interested.conditions) %>%
	filter(unique_name %in% interested.genes) %>%
	filter(type == "mismatch") %>%
	inner_join(
		melted_results %>%
			filter(condition %in% interested.conditions) %>%
			filter(unique_name %in% interested.genes) %>%
			filter(type == "perfect") %>%
			group_by(unique_name, condition) %>%
			select(unique_name, condition, LFC.adj) %>%
			filter(abs(LFC.adj) == max(abs(LFC.adj))) %>%
			rename(response.max = LFC.adj)) %>%
	nest(data = c(-condition, -unique_name, -response.max))

# Test genes and parameters
test.genes <- c("murA", "rpmB", "aroC", "GO593_00515", "glnS", "nuoB", "lpxC")
L.4.parameters <- c("shape", "min_value", "max_value", "kd_50")
BC.5.parameters <- c("shape", "min_value", "max_value", "kd_50", "hormesis")

BC.5_model <- mismatches %>%
	# filter(unique_name %in% test.genes) %>%
	mutate(fit = case_when(
		response.max > 0 ~ map2(data, response.max, ~drm.try(
				data = .x, 
				LFC.adj ~ y_pred, 
				control = drmc(method = "Nelder-Mead", noMessage = TRUE, maxIt = 50000), 
				fct = BC.5.logistic(fixed = c(NA, 0, .y, NA, NA), names = BC.5.parameters))),
		response.max < 0 ~ map2(data, response.max, ~drm.try(
				data = .x, 
				LFC.adj ~ y_pred, 
				control = drmc(method = "Nelder-Mead", noMessage = TRUE, maxIt = 50000), 
				fct = BC.5.logistic(fixed = c(NA, .y, 0, NA, NA), names = BC.5.parameters)))))

# Define BC.5_model and BC.5_reduced_model
BC.5_reduced_model <- mismatches %>%
	# filter(unique_name %in% test.genes) %>%
	mutate(fit = case_when(
		response.max > 0 ~ map2(data, response.max, ~drm.try(
				data = .x, 
				LFC.adj ~ y_pred, 
				control = drmc(method = "Nelder-Mead", noMessage = TRUE, maxIt = 50000), 
				fct = BC.5.logistic(fixed = c(NA, 0, .y, NA, 0), names = BC.5.parameters))),
		response.max < 0 ~ map2(data, response.max, ~drm.try(
				data = .x, 
				LFC.adj ~ y_pred, 
				control = drmc(method = "Nelder-Mead", noMessage = TRUE, maxIt = 50000), 
				fct = BC.5.logistic(fixed = c(NA, .y, 0, NA, 0), names = BC.5.parameters)))))


##########################################################################################
# Main script for full model
BC.5_model_processed <- process_mismatches(BC.5_model)
BC.5_model_summary <- compute_vuln_summary(BC.5_model_processed)
BC.5_model_predictions <- compute_predictions(BC.5_model_processed)
BC.5_model_performance <- compute_model_performance(BC.5_model_processed)
BC.5_model_parameters <- compute_model_parameters(BC.5_model_processed)
BC.5_model_fit_points <- extract_fit_points(BC.5_model_processed)

# Main script for reduced model
BC.5_reduced_processed <- process_mismatches(BC.5_reduced_model)
BC.5_reduced_summary <- compute_vuln_summary(BC.5_reduced_processed)
BC.5_reduced_predictions <- compute_predictions(BC.5_reduced_processed)
BC.5_reduced_performance <- compute_model_performance(BC.5_reduced_processed)
BC.5_reduced_parameters <- compute_model_parameters(BC.5_reduced_processed)
BC.5_reduced_fit_points <- extract_fit_points(BC.5_reduced_processed)

# Save results for full and reduced models
file_names_full <- list(
	vuln_summary = "hormetic_vulnerability_summary_full.tsv.gz",
	fit_predictions = "hormetic_fit_predictions_full.tsv.gz",
	fit_points = "hormetic_fit_points_full.tsv.gz",
	model_performance = "hormetic_performance_full.tsv.gz",
	model_parameters = "hormetic_parameters_full.tsv.gz")

file_names_reduced <- list(
	vuln_summary = "hormetic_vulnerability_summary_reduced.tsv.gz",
	fit_predictions = "hormetic_fit_predictions_reduced.tsv.gz",
	fit_points = "hormetic_fit_points_reduced.tsv.gz",
	model_performance = "hormetic_performance_reduced.tsv.gz",
	model_parameters = "hormetic_parameters_reduced.tsv.gz")

# Save results for full and reduced models
save_results(BC.5_model_summary, BC.5_model_predictions, BC.5_model_fit_points, 
	BC.5_model_performance, BC.5_model_parameters, file_names_full)

save_results(BC.5_reduced_summary, BC.5_reduced_predictions, BC.5_reduced_fit_points, 
	BC.5_reduced_performance, BC.5_reduced_parameters, file_names_reduced)

# Compare the models
model_comparisons <- compare_models(
	BC.5_model,
	BC.5_reduced_model)
