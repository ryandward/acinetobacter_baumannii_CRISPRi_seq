# Brain-Cousens Hormetic DRC Logistic Model (Not LLogistic)
# Ryan Ward
# Tuesday, April 11, 2023

require(pacman)

pacman::p_load(purrr, lmtest, tidyverse, data.table, data.table, tidyverse, broom, modelr, Hmisc, R.utils)
# p_load_current_gh("DoseResponse/drcData", "hrbrmstr/hrbrthemes")
# p_load_current_gh("ryandward/drc", "jokergoo/ComplexHeatmap")


# Load packages
library(drcData)
library(drc)
library(hrbrthemes)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::summarize)


################################################################################
# BC.5 Logistic function model using L.4.parameters and BC.5.parameters
###############################################################################
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



################################################################################
# Functions to extract and predict based on models
################################################################################

# Processes mismatches in the model and calculates p-values, tibble and other related variables
process_mismatches <- function(mismatches) {
	mismatches %>%
		filter(!is.na(fit)) %>%
		mutate(p.vals = map(fit, tidy)) %>%
		mutate(kd_50.tibble = map(p.vals, ~filter(.x, term == 'kd_50') %>%
																select(p.value) %>%
																rename(., vuln.p = p.value))) %>%
		mutate(vuln.tibble = map2(fit, p.vals, ~augment.try(.x, newdata = .y %>% filter(term == "kd_50") %>% select(estimate)) %>%
																rename(., vuln.est = .fitted, vuln.kd_50 = estimate))) %>%
		unnest(vuln.tibble) %>%
		unnest(kd_50.tibble) %>%
		select(-c(p.vals)) %>%
		mutate(Gene = factor(unique_name, levels = unique(interested.genes))) %>%
		mutate(Condition = factor(condition, levels = unique(interested.conditions)))
}

# Computes a summary of vulnerability results
compute_vuln_summary <- function(processed_results) {
	vuln.summary <- processed_results %>% select(unique_name, condition, vuln.est, vuln.kd_50, vuln.p)
	vuln.summary <- processed_results %>% select(condition, unique_name, response.max) %>% inner_join(vuln.summary)
	vuln.summary
}

# Computes predictions for the fitted model
compute_predictions <- function(processed_results) {
	processed_results %>%
		mutate(predictions = map2(fit, data, ~augment.try(.x, newdata = expand.grid(y_pred = seq(0, max(1, max(.y$y_pred)), length = 100)), conf.int = T, conf.level = 0.90))) %>%
		select(Gene, Condition, predictions) %>%
		unnest(predictions)
}

# Extracts data points from the fitted model
extract_fit_points <- function(processed_results) {
	processed_results %>%
		select(Gene, Condition, data) %>%
		unnest(data)
}

# Computes model performance metrics such as AIC, BIC, logLik, and df.residual
compute_model_performance <- function(processed_results) {
	processed_results %>%
		mutate(bayes = map(fit, glance),
					 bayes = map(bayes, ~mutate(.x, logLik = c(logLik)))) %>%
		unnest(bayes) %>%
		select(unique_name, condition, AIC, BIC, logLik, df.residual)
}

# Computes model parameters and their statistics
compute_model_parameters <- function(processed_results) {
	processed_results %>%
		mutate(perf = map(fit, tidy)) %>%
		unnest(perf) %>%
		select(unique_name, condition, term, estimate, std.error, statistic, p.value) %>%
		mutate_if(is.numeric, round, 3)
}

################################################################################
# Functions to compare between models
################################################################################

# Modify calculate_lrt function to accept any two models
calculate_lrt <- function(this.gene, this.condition, filtered_HA, filtered_H0) {
	
	# Get the model fits for the given gene and condition
	this.HA <- filtered_HA %>%
		filter(unique_name == this.gene, condition == this.condition) %>%
		pull(fit)
	
	this.H0 <- filtered_H0 %>%
		filter(unique_name == this.gene, condition == this.condition) %>%
		pull(fit)
	
	# Perform the likelihood ratio test
	lrt.result <- lrtest(this.HA[[1]], this.H0[[1]])
	
	# Return the likelihood ratio test p-value
	return(lrt.result$'Pr(>Chisq)'[2])
}

# Compare between two DRC models
compare_models <- function(HA_model, H0_model) {
	
	# Get all unique gene and condition combinations from both models
	HA_gene_condition_combos <- HA_model %>%
		select(unique_name, condition) %>%
		distinct()
	
	H0_gene_condition_combos <- H0_model %>%
		select(unique_name, condition) %>%
		distinct()
	
	# Find gene-condition combinations present in both models
	common_gene_condition_combos <- inner_join(HA_gene_condition_combos, H0_gene_condition_combos, by = c("unique_name", "condition"))
	
	# Filter both data frames once
	filtered_HA <- HA_model %>% select(unique_name, condition, fit)
	filtered_H0 <- H0_model %>% select(unique_name, condition, fit)
	
	# Initialize an empty data frame to store the results
	results <- tibble(unique_name = character(), condition = character(), lrt_p_value = numeric())
	
	# Iterate through each common gene and condition combination
	for (i in seq_along(common_gene_condition_combos$unique_name)) {
		this.gene <- common_gene_condition_combos$unique_name[i]
		this.condition <- common_gene_condition_combos$condition[i]
		
		# Calculate LRT p-values with error handling
		lrt_p_value <- tryCatch({
			calculate_lrt(this.gene, this.condition, filtered_HA, filtered_H0)
		}, error = function(e) {
			NA
		})
		

		# Store the p-values
		results <- results %>%
		add_row(unique_name = this.gene, condition = this.condition, lrt_p_value = lrt_p_value)

	}
	
	results <- results %>% mutate_if(is.numeric, round, 3)
	
	return(results)
}

################################################################################
# Functions to check, save, and read results
################################################################################

# Check if all files in the list exist
check_files_exist <- function(file_names, output_dir = "Results") {
	all_files_exist <- TRUE
	for (file_name in file_names) {
		if (!file.exists(file.path(output_dir, file_name))) {
			all_files_exist <- FALSE
			break
		}
	}
	return(all_files_exist)
}


# Function to check if all files in the list exist, excluding the last one (drc_fits)
check_processed_files_exist <- function(file_names, output_dir = "Results") {
	existing_files <- list()
	for (i in 1:(length(file_names) - 1)) {
		if (file.exists(file.path(output_dir, file_names[[i]]))) {
			existing_files[[i]] <- file_names[[i]]
		}
	}
	return(existing_files)
}


# Save results function including DRC objects
save_results <- function(results, file_names, output_dir = "Results") {
	if (!dir.exists(output_dir)) {dir.create(output_dir)}
	
	fwrite(results$vuln.summary, file.path(output_dir, file_names$vuln_summary), sep = "\t")
	fwrite(results$fit_predictions, file.path(output_dir, file_names$fit_predictions))
	fwrite(results$fit_points, file.path(output_dir, file_names$fit_points))
	fwrite(results$model_performance, file.path(output_dir, file_names$model_performance), sep = "\t")
	fwrite(results$model_parameters, file.path(output_dir, file_names$model_parameters), sep = "\t")
	# saveRDS(results$drc_fits, file.path(output_dir, file_names$drc_fits))
}

# Read results from the saved files
read_results <- function(file_names, output_dir = "Results", exclude_drc_fits = FALSE) {
	vuln.summary <- fread(file.path(output_dir, file_names$vuln_summary), sep = "\t")
	fit_predictions <- fread(file.path(output_dir, file_names$fit_predictions))
	fit_points <- fread(file.path(output_dir, file_names$fit_points))
	model_performance <- fread(file.path(output_dir, file_names$model_performance), sep = "\t")
	model_parameters <- fread(file.path(output_dir, file_names$model_parameters), sep = "\t")
	
	drc_fits <- NULL
	if (!exclude_drc_fits && file.exists(file.path(output_dir, file_names$drc_fits))) {
		drc_fits <- readRDS(file.path(output_dir, file_names$drc_fits))
	}
	
	return(list(
		vuln.summary = vuln.summary,
		fit_predictions = fit_predictions,
		fit_points = fit_points,
		model_performance = model_performance,
		model_parameters = model_parameters,
		drc_fits = drc_fits
	))
}


check_file_exist <- function(file_path) {
	return(file.exists(file_path))
}


check_and_load_model_comparisons <- function(HA, H0, file_name, directory = "Results") {
	file_path <- file.path(directory, file_name)
	
	if (check_file_exist(file_path)) {
		message(paste("File found. Loading from disk..."))
		loaded_object <- fread(file_path, sep = "\t")
		message("Loaded successfully.")
	} else {
		message("File not found on disk. Calculating...")
		loaded_object <- calculate_model_comparisons(HA, H0, directory)
	}
	
	return(loaded_object)
}


calculate_model_comparisons <- function(full_results, reduced_results, object_name, directory = "Results") {
	model_comparisons <- compare_models(
		full_results$drc_fits,
		reduced_results$drc_fits
	)
	
	# closeAllConnections()  # turn off sink FIX LATER
	
	model_comparisons <- inner_join(
		full_results$model_performance %>% select(unique_name, condition, logLik) %>% rename(hormetic_logLik = logLik),
		reduced_results$model_performance %>% select(unique_name, condition, logLik) %>% rename(reduced_logLik = logLik)
	) %>% inner_join(model_comparisons)
	
	assign(object_name, model_comparisons, envir = .GlobalEnv)
	# fwrite(model_comparisons, paste0(directory, "/", file_names$drc_fits), sep = "\t")
}
