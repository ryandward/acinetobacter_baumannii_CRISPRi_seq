# Load necessary packages
source("packages.R")
source("drc_logistic_functions.R")


# Load data in
aba_key <- fread("aba_key.tsv")

interested.genes <- fread("curated_names.tsv") %>% pull(unique_name)

interested.conditions <- c( 
	# List of conditions that make sense to analyze, i.e., all versus T0 without induction.
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

# The maximum knockdown dose provided by the system
max_y_pred <- aba_key %>% select(y_pred) %>% summarize(max(y_pred, na.rm = TRUE)) %>% as.numeric

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
glance.try <- possibly(glance, otherwise = tibble(logLik = NA_real_))
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

# L.4.parameters <- c("hill", "min_value", "max_value", "kd_50")
# Define the parameter names for the BC.5 model
BC.5.parameters <- c("shape", "min_value", "max_value", "kd_50", "hormesis")

# This function creates constraints for the model parameters.
# It takes four lists as input: fixed_params, lowerl_params, upperl_params, and start_params.
# fixed_params: Parameters with fixed values.
# lowerl_params: Parameters with lower limits.
# upperl_params: Parameters with upper limits.
# start_params: Parameters with starting values for optimization.
create_constraints <- function(fixed_params, lowerl_params, upperl_params, start_params) {
	param_names <- BC.5.parameters
	
	# Initialize the constraints and start values with default values
	fixed_constraints <- setNames(rep(NA, length(param_names)), param_names)
	lowerl_constraints <- setNames(rep(-Inf, length(param_names)), param_names)
	upperl_constraints <- setNames(rep(Inf, length(param_names)), param_names)
	start_values <- setNames(rep(NA, length(param_names)), param_names)
	
	# Update the constraints and start values with the input values
	fixed_constraints[names(fixed_params)] <- fixed_params
	lowerl_constraints[names(lowerl_params)] <- lowerl_params
	upperl_constraints[names(upperl_params)] <- upperl_params
	start_values[names(start_params)] <- start_params
	
	# Remove fixed parameters from the other constraint lists
	non_fixed_names <- setdiff(param_names, names(fixed_params))
	
	# Return the final constraints and start values
	return(list(fixed = unlist(fixed_constraints),
							lowerl = unlist(lowerl_constraints[non_fixed_names]),
							upperl = unlist(upperl_constraints[non_fixed_names]),
							start = unlist(start_values[non_fixed_names])))
}

# This function returns the model constraints based on the model number and response_max.
# model_number: An integer representing the model type (1, 2, 3, or 4).
# response_max: The maximum response value for the data.
get_model_constraints <- function(model_number, response_max) {
	if (model_number == 1) {
		# Model 1: Positive response with a lower limit of 0 for the response.
		# min_value is fixed at 0, which implies that the response starts at 0 and increases.
		fixed_params = list(min_value = 0)
		lowerl_params = list(max_value = 0, kd_50 = 0)
		upperl_params = list(shape = 0, max_value = response_max, kd_50 = max_y_pred)
		start_params = list(shape = 0, max_value = response_max, kd_50 = 0.5, hormesis = 0)
	} else if (model_number == 2) {
		# Model 2: Negative response with an upper limit of 0 for the response.
		# max_value is fixed at 0, which implies that the response starts at a high level and decreases.
		fixed_params = list(max_value = 0)
		lowerl_params = list(shape = 0, min_value = response_max, kd_50 = 0)
		upperl_params = list(min_value = 0, kd_50 = max_y_pred)
		start_params = list(shape = 0, min_value = response_max, kd_50 = 0.5, hormesis = 0)
	} else if (model_number == 3) {
		# Model 3: Positive response without hormesis.
		# min_value and hormesis are fixed at 0, which implies that the response starts at 0 and increases,
		# and there is no hormesis effect.
		fixed_params = list(min_value = 0, hormesis = 0)
		lowerl_params = list(max_value = 0, kd_50 = 0)
		upperl_params = list(shape = 0, max_value = response_max, kd_50 = max_y_pred)
		start_params = list(shape = 0, max_value = response_max, kd_50 = 0.5)
	} else {
		# Model 4: Negative response without hormesis.
		# max_value and hormesis are fixed at 0, which implies that the response starts at a high level and decreases,
		# and there is no hormesis effect.
		fixed_params = list(max_value = 0, hormesis = 0)
		lowerl_params = list(shape = 0, min_value = response_max, kd_50 = 0)
		upperl_params = list(min_value = 0, kd_50 = max_y_pred)
		start_params = list(shape = 0, min_value = response_max, kd_50 = 0.5)
	}
	return(create_constraints(fixed_params, lowerl_params, upperl_params, start_params))
}

# Check, load, or calculate the models

message("Checking if BC.5_model results are in memory...\n")
if (exists("BC.5_model")) {
	message("BC.5_model results found in memory.\n")
} else {
	message("Checking if BC.5_model results files exist...\n")
	if (check_files_exist(file_names_full)) {
		message("BC.5_model results files found. Loading results...\n")
		BC.5_model <- read_results(file_names_full)
		message("BC.5_model results loaded successfully.\n")
	} else {
		message("BC.5_model results files not found. Proceeding with model fitting...\n")
		total <- length(mismatches$data) * 2
		count <- 1
		
		BC.5_model <- mismatches %>%
			mutate(
				fit_model_1 = map2(data, response.max, ~{
					cat(sprintf("\rFitting Model 1 and 2: %d/%d... ", count, total))
					model_constraints_1 <- get_model_constraints(1, .y)
					result <- drm.try(
						data = .x,
						LFC.adj ~ y_pred,
						control = drmc(method = "L-BFGS-B", maxIt = 1e8, relTol = 1e-25),
						lowerl = unlist(model_constraints_1$lowerl),
						upperl = unlist(model_constraints_1$upperl),
						start = unlist(model_constraints_1$start),
						fct = BC.5(fixed = unlist(model_constraints_1$fixed), names = BC.5.parameters))
					count <<- count + 1
					result
				}),
				fit_model_2 = map2(data, response.max, ~{
					cat(sprintf("\rFitting Model 1 and 2: %d/%d... ", count, total))
					model_constraints_2 <- get_model_constraints(2, .y)
					result <- drm.try(
						data = .x,
						LFC.adj ~ y_pred,
						control = drmc(method = "L-BFGS-B", maxIt = 1e8, relTol = 1e-25),
						lowerl = unlist(model_constraints_2$lowerl),
						upperl = unlist(model_constraints_2$upperl),
						start = unlist(model_constraints_2$start),
						fct = BC.5(fixed = unlist(model_constraints_2$fixed), names = BC.5.parameters))
					count <<- count + 1
					result
				})
			)
		
		
		
		BC.5_model_comparison <- BC.5_model %>%
			mutate(
				logLik_model_1 = map_dbl(fit_model_1, ~glance.try(.x)$logLik),
				logLik_model_2 = map_dbl(fit_model_2, ~glance.try(.x)$logLik),
				better_model = case_when(
					is.na(logLik_model_1) & is.na(logLik_model_2) ~ NA_character_,
					is.na(logLik_model_1) ~ "model_2",
					is.na(logLik_model_2) ~ "model_1",
					logLik_model_1 > logLik_model_2 ~ "model_1",
					TRUE ~ "model_2"
				)
			) %>%
			select(unique_name, condition, logLik_model_1, logLik_model_2, better_model)
		
		BC.5_model <- BC.5_model %>%
			mutate(
				fit = if_else(
					BC.5_model_comparison$better_model == "model_1",
					fit_model_1,
					fit_model_2
				)
			) %>% select(-fit_model_1, -fit_model_2)
		message("...Done!")
	}
}

message("Checking if BC.5_reduced_model results are in memory...\n")
if (exists("BC.5_reduced_model")) {
	message("BC.5_reduced_model results found in memory.\n")
} else {
	message("Checking if BC.5_reduced_model results files exist...\n")
	if (check_files_exist(file_names_reduced)) {
		message("BC.5_reduced_model results files found. Loading results...\n")
		BC.5_reduced_model <- read_results(file_names_reduced)
		message("BC.5_reduced_model results loaded successfully.\n")
	} else {
		message("BC.5_reduced_model results files not found. Proceeding with model fitting...\n")
		total <- length(mismatches$data) * 2
		count <- 1
		
		BC.5_reduced_model <- mismatches %>%
			mutate(
				fit_model_3 = map2(data, response.max, ~{
					cat(sprintf("\rFitting Model 3 and 4: %d/%d... ", count, total))
					model_constraints_3 <- get_model_constraints(3, .y)
					result <- drm.try(
						data = .x,
						LFC.adj ~ y_pred,
						control = drmc(method = "L-BFGS-B", maxIt = 1e8, relTol = 1e-25),
						lowerl = unlist(model_constraints_3$lowerl),
						upperl = unlist(model_constraints_3$upperl),
						start = unlist(model_constraints_3$start),
						fct = BC.5(fixed = unlist(model_constraints_3$fixed), names = BC.5.parameters))
					count <<- count + 1
					result
				}),
				fit_model_4 = map2(data, response.max, ~{
					cat(sprintf("\rFitting Model 3 and 4: %d/%d... ", count, total))
					model_constraints_4 <- get_model_constraints(4, .y)
					result <- drm.try(
						data = .x,
						LFC.adj ~ y_pred,
						control = drmc(method = "L-BFGS-B", maxIt = 1e8, relTol = 1e-25),
						lowerl = unlist(model_constraints_4$lowerl),
						upperl = unlist(model_constraints_4$upperl),
						start = unlist(model_constraints_4$start),
						fct = BC.5(fixed = unlist(model_constraints_4$fixed), names = BC.5.parameters))
					count <<- count + 1
					result
				})
			)
		
		BC.5_reduced_model_comparison <- BC.5_reduced_model %>%
			mutate(
				logLik_model_3 = map_dbl(fit_model_3, ~glance.try(.x)$logLik),
				logLik_model_4 = map_dbl(fit_model_4, ~glance.try(.x)$logLik),
				better_model = case_when(
					is.na(logLik_model_3) & is.na(logLik_model_4) ~ NA_character_,
					is.na(logLik_model_3) ~ "model_4",
					is.na(logLik_model_4) ~ "model_3",
					logLik_model_3 > logLik_model_4 ~ "model_3",
					TRUE ~ "model_4"
				)
			) %>%
			select(unique_name, condition, logLik_model_3, logLik_model_4, better_model)
		
		
		BC.5_reduced_model <- BC.5_reduced_model %>%
			mutate(
				fit = if_else(
					BC.5_reduced_model_comparison$better_model == "model_3",
					fit_model_3,
					fit_model_4
				)
			) %>% select(-fit_model_3, -fit_model_4)
		message("...Done!")
		
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
file_names <- list(
	drc_fits = "hormetic_model_comparisons.tsv"
)

object_name <- "model_comparisons"
check_and_load_model_comparisons(file_names, object_name)

#################################################################################

# Create responses data frame
response_data <- full_results$fit_predictions %>%
	select(Gene, Condition, y_pred, .fitted) %>%
	rename(unique_name = Gene, condition = Condition, fitted = .fitted) %>%
	filter(y_pred != 0) %>%
	group_by(unique_name, condition) %>%
	summarize(
		highest_fitted = max(fitted),
		lowest_fitted = min(fitted),
		highest_y_pred = y_pred[which.max(fitted)],
		lowest_y_pred = y_pred[which.min(fitted)],
		.groups = 'drop'
	)

# Create intermediate_phenotypes data frame
intermediate_phenotypes <- response_data %>%
	inner_join(mismatches %>% select(-data)) %>%
	mutate(opposite_direction = sign(highest_fitted) != sign(response.max) | sign(lowest_fitted) != sign(response.max)) %>%
	arrange(highest_y_pred)

# Create filtered_data data frame
filtered_results <- model_comparisons %>%
	full_join(full_results$model_parameters %>% filter(term == "hormesis")) %>%
	filter(p.value <= 0.05 & lrt_p_value <= 0.05) %>%
	inner_join(intermediate_phenotypes) %>%
	filter(opposite_direction == TRUE)

# Create full_estimates data frame
full_estimates <- full_results$model_parameters %>% data.table %>% dcast(unique_name + condition ~ term, value.var = "estimate")

# Create hormesis_with_parameters data frame
hormesis_results <- filtered_results %>%
	select(-term, -statistic, -std.error, -estimate) %>%
	rename(hormesis_p_value = p.value) %>%
	inner_join(full_estimates) %>% tibble

##########################################################################################

model_comparisons %>% fwrite("Results/model_comparisons.tsv", sep = "\t")
hormesis_results %>% fwrite("Results/model_comparisons_hormesis.tsv", sep = "\t")

