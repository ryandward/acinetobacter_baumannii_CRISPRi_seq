# Brain-Cousens Hormetic DRC Logistic Model (Not LLogistic)
# Ryan Ward
# Tuesday, April 11, 2023

require(pacman)

pacman::p_load(lmtest, tidyverse, data.table)

p_load(data.table, tidyverse, broom, modelr, Hmisc)
# p_load_current_gh("DoseResponse/drcData", "ryandward/drc", "hrbrmstr/hrbrthemes")

# Load packages
library(drcData)
library(drc)
library(hrbrthemes)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::summarize)

L.4.parameters <- c("hill", "min_value", "max_value", "kd_50")
BC.5.parameters <- c("shape", "min_value", "max_value", "kd_50", "hormesis")

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
compute_vuln_summary <- function(mismatches) {
	vuln.summary <- mismatches %>% select(unique_name, condition, vuln.est, vuln.kd_50, vuln.p)
	vuln.summary <- mismatches %>% select(condition, unique_name, response.max) %>% inner_join(vuln.summary)
	vuln.summary
}

# Computes predictions for the fitted model
compute_predictions <- function(mismatches) {
	mismatches %>%
		mutate(predictions = map2(fit, data, ~augment.try(.x, newdata = expand.grid(y_pred = seq(0, max(1, max(.y$y_pred)), length = 250)), conf.int = T, conf.level = 0.90))) %>%
		select(Gene, Condition, predictions) %>%
		unnest(predictions)
}

# Extracts data points from the fitted model
extract_fit_points <- function(mismatches) {
	mismatches %>%
		select(Gene, Condition, data) %>%
		unnest(data)
}

# Computes model performance metrics such as AIC, BIC, logLik, and df.residual
compute_model_performance <- function(mismatches) {
	mismatches %>%
		mutate(bayes = map(fit, glance),
					 bayes = map(bayes, ~mutate(.x, logLik = c(logLik)))) %>%
		unnest(bayes) %>%
		select(unique_name, condition, AIC, BIC, logLik, df.residual)
}

# Computes model parameters and their statistics
compute_model_parameters <- function(mismatches) {
	mismatches %>%
		mutate(perf = map(fit, tidy)) %>%
		unnest(perf) %>%
		select(unique_name, condition, term, estimate, std.error, statistic, p.value)
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

# Modify calculate_anova function to accept any two models
calculate_anova <- function(this.gene, this.condition, filtered_HA, filtered_H0) {
	
	# Get the model fits for the given gene and condition
	this.HA <- filtered_HA %>%
		filter(unique_name == this.gene, condition == this.condition) %>%
		pull(fit)
	
	this.H0 <- filtered_H0 %>%
		filter(unique_name == this.gene, condition == this.condition) %>%
		pull(fit)
	
	# Perform the likelihood ratio test and suppress the output
	sink(tempfile())
	anova.result <- anova(this.HA[[1]], this.H0[[1]])
	sink()
	
	# Restore console output and return the ANOVA p-value
	return(anova.result$'p value'[2])
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
	results <- tibble(unique_name = character(), condition = character(), lrt_p_value = numeric(), anova_p_value = numeric())
	
	# Iterate through each common gene and condition combination
	for (i in seq_along(common_gene_condition_combos$unique_name)) {
		this.gene <- common_gene_condition_combos$unique_name[i]
		this.condition <- common_gene_condition_combos$condition[i]
		
		# Calculate LRT and ANOVA p-values with error handling
		lrt_p_value <- tryCatch({
			calculate_lrt(this.gene, this.condition, filtered_HA, filtered_H0)
		}, error = function(e) {
			NA
		})
		
		anova_p_value <- tryCatch({
			calculate_anova(this.gene, this.condition, filtered_HA, filtered_H0)
		}, error = function(e) {
			NA
		})
		
		# Store the p-values
		results <- results %>%
			add_row(unique_name = this.gene, condition = this.condition, lrt_p_value = lrt_p_value, anova_p_value = anova_p_value)
	}
	
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

# Save results function including DRC objects
save_results <- function(results, file_names, output_dir = "Results") {
	if (!dir.exists(output_dir)) {dir.create(output_dir)}
	
	fwrite(results$vuln.summary, file.path(output_dir, file_names$vuln_summary), sep = "\t")
	fwrite(results$fit_predictions, file.path(output_dir, file_names$fit_predictions))
	fwrite(results$fit_points, file.path(output_dir, file_names$fit_points))
	fwrite(results$model_performance, file.path(output_dir, file_names$model_performance), sep = "\t")
	fwrite(results$model_parameters, file.path(output_dir, file_names$model_parameters), sep = "\t")
	saveRDS(results$drc_fits, file.path(output_dir, file_names$drc_fits))
}

# Read results from the saved files
read_results <- function(file_names, output_dir = "Results") {
	vuln.summary <- fread(file.path(output_dir, file_names$vuln_summary), sep = "\t")
	fit_predictions <- fread(file.path(output_dir, file_names$fit_predictions))
	fit_points <- fread(file.path(output_dir, file_names$fit_points))
	model_performance <- fread(file.path(output_dir, file_names$model_performance), sep = "\t")
	model_parameters <- fread(file.path(output_dir, file_names$model_parameters), sep = "\t")
	
	drc_fits <- NULL
	if (file.exists(file.path(output_dir, file_names$drc_fits))) {
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


check_and_load_model_comparisons <- function(file_names, object_name, directory = "Results") {
	file_path <- file.path(directory, file_names$drc_fits)
	
	if (exists(object_name)) {
		message(paste(object_name, "found in memory."))
	} else if (check_file_exist(file_path)) {
		message(paste(object_name, "file found. Loading from disk..."))
		assign(object_name, fread(file_path, sep = "\t"), envir = .GlobalEnv)
		message(paste(object_name, "loaded successfully."))
	} else {
		message(paste(object_name, "not found in memory or on disk. Calculating..."))
		calculate_model_comparisons(full_results, reduced_results, object_name)
	}
}


calculate_model_comparisons <- function(full_results, reduced_results, object_name) {
	model_comparisons <- compare_models(
		full_results$drc_fits,
		reduced_results$drc_fits
	)
	
	closeAllConnections()  # turn off sink FIX LATER
	
	model_comparisons <- inner_join(
		full_results$model_performance %>% select(unique_name, condition, logLik) %>% rename(hormetic_logLik = logLik),
		reduced_results$model_performance %>% select(unique_name, condition, logLik) %>% rename(reduced_logLik = logLik)
	) %>% inner_join(model_comparisons)
	
	assign(object_name, model_comparisons, envir = .GlobalEnv)
	fwrite(model_comparisons, file_names$drc_fits, sep = "\t")
}
