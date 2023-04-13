# Load necessary packages
source("drc_logistic_functions.R")

require(conflicted)
require(pacman)
require(progress)

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

# Test genes and parameters
test.genes <- c("murA", "rpmB", "aroC", "GO593_00515", "glnS", "nuoB", "lpxC")
test.genes <- c(test.genes, c("rpmF", "glnS", "trpS", "acpT"))

L.4.parameters <- c("shape", "min_value", "max_value", "kd_50")
BC.5.parameters <- c("shape", "min_value", "max_value", "kd_50", "hormesis")

# Estimate dose-response models
mismatches <- melted_results %>%
	# filter(unique_name %in% test.genes) %>%
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

##########################################################################################

# Full Model
message("Checking if full model results are in memory...\n")
if (exists("full_results")) {
	message("Full model results found in memory.\n")
} else {
	message("Checking if full model results files exist...\n")
	if (check_files_exist(file_names_full)) {
		message("Full model results files found. Loading results...\n")
		full_results <- read_results(file_names_full)
		message("Full model results loaded successfully.\n")
	} else {
		message("Full model results files not found. Proceeding with model fitting...\n")
		total <- length(mismatches$data)
		count <- 1
		
		BC.5_model <- mismatches %>%
			mutate(fit = map2(data, response.max, ~{
				message(sprintf("\rFitting full model: %d/%d", count, total))
				result <- if (.y > 0) {
					drm.try(
						data = .x, 
						LFC.adj ~ y_pred, 
						# control = drmc(method = "Nelder-Mead", maxIt = 1e7, relTol = 1e-25),
						control = drmc(method = "L-BFGS-B", maxIt = 1e7, relTol = 1e-25),
						lowerl = c(-50, 0, -20), upperl = c(50, 5, 20),
						start = c(1, 0.5, 0),
						fct = BC.5(fixed = c(NA, 0, .y, NA, NA), names = BC.5.parameters))
				} else {
					drm.try(
						data = .x, 
						LFC.adj ~ y_pred, 
						# control = drmc(method = "Nelder-Mead", maxIt = 1e7, relTol = 1e-25),
						control = drmc(method = "L-BFGS-B", maxIt = 1e7, relTol = 1e-25),
						lowerl = c(-50, 0, -20), upperl = c(50, 5, 20),
						start = c(1, 0.5, 0),
						fct = BC.5(fixed = c(NA, .y, 0, NA, NA), names = BC.5.parameters))
				}
				count <<- count + 1
				result
			}))
		message("\nModel fitting complete.\n")
	}
}

# Reduced Model
message("Checking if reduced model results are in memory...\n")
if (exists("reduced_results")) {
	message("Reduced model results found in memory.\n")
} else {
	message("Checking if reduced model results files exist...\n")
	if (check_files_exist(file_names_reduced)) {
		message("Reduced model results files found. Loading results...\n")
		reduced_results <- read_results(file_names_reduced)
		message("Reduced model results loaded successfully.\n")
	} else {
		message("Reduced model results files not found. Proceeding with model fitting...\n")
		total <- length(mismatches$data)
		count <- 1
		
		BC.5_reduced_model <- mismatches %>%
			mutate(fit = map2(data, response.max, ~{
				message(sprintf("\rFitting reduced model: %d/%d", count, total))
				result <- if (.y > 0) {
					drm.try(
						data = .x, 
						LFC.adj ~ y_pred, 
						# control = drmc(method = "Nelder-Mead", maxIt = 1e7, relTol = 1e-25),
						control = drmc(method = "L-BFGS-B", maxIt = 1e7, relTol = 1e-25),
						lowerl = c(-50, 0), upperl = c(50, 5),
						start = c(1, 0.5),
						fct = BC.5(fixed = c(NA, 0, .y, NA, 0), names = BC.5.parameters))
				} else {
					drm.try(
						data = .x, 
						LFC.adj ~ y_pred, 
						# control = drmc(method = "Nelder-Mead", maxIt = 1e7, relTol = 1e-25),
						control = drmc(method = "L-BFGS-B", maxIt = 1e7, relTol = 1e-25),
						lowerl = c(-50, 0), upperl = c(50, 5),
						start = c(1, 0.5),
						fct = BC.5(fixed = c(NA, .y, 0, NA, 0), names = BC.5.parameters))
				}
				count <<- count + 1
				result
			}))
		message("\nModel fitting complete.\n")
	}
}

						

##########################################################################################
message("Beginning to process full model...\n")

if (exists("full_results")) {
	message("Full model results already loaded in memory.\n")
} else if (check_files_exist(file_names_full)) {
	message("Loading full model results from disk...\n")
	full_results <- read_results(file_names_full)
	message("Full model results loaded successfully.\n")
} else {
	message("Processing BC.5_model...")
	BC.5_model_processed <- process_mismatches(BC.5_model)
	drc_fits <- BC.5_model_processed %>% select(unique_name, condition, fit)
	
	message("Computing results for full model...\n")
	full_results <- list(
		vuln.summary = compute_vuln_summary(BC.5_model_processed),
		fit_predictions = compute_predictions(BC.5_model_processed),
		fit_points = extract_fit_points(BC.5_model_processed),
		model_performance = compute_model_performance(BC.5_model_processed),
		model_parameters = compute_model_parameters(BC.5_model_processed),
		drc_fits = drc_fits
	)
	save_results(full_results, file_names_full)
}

message("Full model processed.\n")

message("Beginning to process reduced model...\n")
if (exists("reduced_results")) {
	message("Reduced model results already loaded in memory.\n")
} else if (check_files_exist(file_names_reduced)) {
	message("Loading reduced model results from disk...\n")
	reduced_results <- read_results(file_names_reduced)
	message("Reduced model results loaded successfully.\n")
} else {
	message("Processing BC.5_reduced_model...\n")
	BC.5_reduced_processed <- process_mismatches(BC.5_reduced_model)
	drc_fits <- BC.5_reduced_processed %>% select(unique_name, condition, fit)
	
	message("Computing results for reduced model...\n")
	reduced_results <- list(
		vuln.summary = compute_vuln_summary(BC.5_reduced_processed),
		fit_predictions = compute_predictions(BC.5_reduced_processed),
		fit_points = extract_fit_points(BC.5_reduced_processed),
		model_performance = compute_model_performance(BC.5_reduced_processed),
		model_parameters = compute_model_parameters(BC.5_reduced_processed),
		drc_fits = drc_fits
	)
	save_results(reduced_results, file_names_reduced)
}

message("Reduced model processed.\n")


################################################################################
# Compare the models
model_comparisons_file_name <- "hormetic_model_comparisons.tsv"
check_and_load_model_comparisons(model_comparisons_file_name)


if (is.null(model_comparisons)) {
	message("Calculating model comparisons...\n")
	model_comparisons <- compare_models(
		full_results$drc_fits,
		reduced_results$drc_fits
	)
	
	closeAllConnections()  # turn off sink FIX LATER
	
	model_comparisons <- inner_join(
		full_results$model_performance %>% select(unique_name, condition, logLik) %>% rename(hormetic_logLik = logLik),
		reduced_results$model_performance %>% select(unique_name, condition, logLik) %>% rename(reduced_logLik = logLik)
	) %>% inner_join(model_comparisons)
	
	fwrite(model_comparisons, file.path("Results", model_comparisons_file_name), sep = "\t")
}
