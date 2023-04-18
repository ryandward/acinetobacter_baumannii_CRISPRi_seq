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

linear_BC.5 <- function(
		fixed = c(NA, NA, NA, NA, NA), 
		names = c("b", "c", "d", "e", "f")){
	return(linear_braincousens(fixed = fixed, names = names))
}

linear_braincousens <- function (fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", 
																																					"d", "e", "f"), method = c("1", "2", "3", "4"), ssfct = NULL, 
																 fctName, fctText) 
{
	numParm <- 5
	if (!is.character(names) | !(length(names) == numParm)) {
		stop("Not correct 'names' argument")
	}
	if (!(length(fixed) == numParm)) {
		stop("Not correct 'fixed' argument")
	}
	notFixed <- is.na(fixed)
	parmVec <- rep(0, numParm)
	parmVec[!notFixed] <- fixed[!notFixed]
	parmVec1 <- parmVec
	parmVec2 <- parmVec
	fct <- function(dose, parm) {
		parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
		parmMat[, notFixed] <- parm
		parmMat[, 2] + (parmMat[, 3] + parmMat[, 5] * dose - 
											parmMat[, 2])/(1 + exp(parmMat[, 1] * (dose - 
																														 	parmMat[, 4])))
	}
	if (FALSE) {
		ssfct <- function(dataFra) {
			dose2 <- dataFra[, 1]
			resp3 <- dataFra[, 2]
			startVal <- rep(0, numParm)
			startVal[3] <- max(resp3) + 0.001
			startVal[2] <- min(resp3) - 0.001
			startVal[5] <- 0
			if (length(unique(dose2)) == 1) {
				return((c(NA, NA, startVal[3], NA, NA))[notFixed])
			}
			indexT2 <- (dose2 > 0)
			if (!any(indexT2)) {
				return((rep(NA, numParm))[notFixed])
			}
			dose3 <- dose2[indexT2]
			resp3 <- resp3[indexT2]
			logitTrans <- log((startVal[3] - resp3)/(resp3 - 
																							 	startVal[2] + 0.001))
			logitFit <- lm(logitTrans ~ dose3)
			startVal[4] <- (-coef(logitFit)[1]/coef(logitFit)[2])
			startVal[1] <- coef(logitFit)[2]
			return(startVal[notFixed])
		}
	}
	if (!is.null(ssfct)) {
		ssfct <- ssfct
	}
	else {
		ssfct <- function(dframe) {
			initval <- llogistic()$ssfct(dframe)
			initval[5] <- 0
			return(initval[notFixed])
		}
	}
	names <- names[notFixed]
	deriv1 <- function(dose, parm) {
		parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
		parmMat[, notFixed] <- parm
		t1 <- parmMat[, 3] - parmMat[, 2] + parmMat[, 5] * dose
		t2 <- exp(parmMat[, 1] * (dose - parmMat[, 4]))
		t3 <- 1 + t2
		t4 <- (1 + t2)^(-2)
		cbind(-t1 * t2 * t4, 1 - 1/t3, 1/t3, t1 * t2 * parmMat[, 1] * t4, dose/t3)[, notFixed]
	}
	deriv2 <- NULL
	edfct <- function(parm, respl, reference, type, lower = 0.001, 
										upper = 1000, ...) {
		interval <- c(lower, upper)
		parmVec[notFixed] <- parm
		p <- EDhelper(parmVec, respl, reference, type)
		tempVal <- (100 - p)/100
		helpEqn <- function(dose) {
			expVal <- exp(parmVec[1] * (dose - parmVec[4]))
			parmVec[5] * (1 + expVal * (1 - parmVec[1])) - (parmVec[3] - 
																												parmVec[2]) * expVal * parmVec[1]/dose
		}
		maxAt <- uniroot(helpEqn, interval)$root
		eqn <- function(dose) {
			tempVal * (1 + exp(parmVec[1] * (dose - parmVec[4]))) - 
				(1 + parmVec[5] * dose/(parmVec[3] - parmVec[2]))
		}
		EDp <- uniroot(eqn, lower = maxAt, upper = upper)$root
		EDdose <- EDp
		tempVal1 <- exp(parmVec[1] * (EDdose - parmVec[4]))
		tempVal2 <- parmVec[3] - parmVec[2]
		derParm <- c(tempVal * tempVal1 * (EDdose - parmVec[4]), 
								 -parmVec[5] * EDdose/((tempVal2)^2), parmVec[5] * 
								 	EDdose/((tempVal2)^2), -tempVal * tempVal1 * 
								 	parmVec[1], -EDdose/tempVal2)
		derDose <- tempVal * tempVal1 * parmVec[1] - 
			parmVec[5]/tempVal2
		EDder <- derParm/derDose
		return(list(EDp, EDder[notFixed]))
	}
	maxfct <- function(parm, lower = 0.001, upper = 1000) {
		parmVec[notFixed] <- parm
		if (parmVec[1] < 1) {
			stop("Brain-Cousens model with b<1 not meaningful")
		}
		if (parmVec[5] < 0) {
			stop("Brain-Cousens model with f<0 not meaningful")
		}
		optfct <- function(t) {
			expTerm1 <- parmVec[5] * t
			expTerm2 <- exp(parmVec[1] * (t - parmVec[4]))
			return(parmVec[5] * (1 + expTerm2) - (parmVec[3] - 
																							parmVec[2] + expTerm1) * expTerm2 * parmVec[1]/t)
		}
		ED1 <- edfct(parm, 1, lower, upper)[[1]]
		doseVec <- exp(seq(log(1e-06), log(ED1), length = 100))
		maxDose <- uniroot(optfct, c((doseVec[optfct(doseVec) > 
																						0])[1], ED1))$root
		return(c(maxDose, fct(maxDose, matrix(parm, 1, length(names)))))
	}
	returnList <- list(fct = fct, ssfct = ssfct, names = names, 
										 deriv1 = deriv1, deriv2 = deriv2, edfct = edfct, maxfct = maxfct, 
										 name = ifelse(missing(fctName), as.character(match.call()[[1]]), 
										 							fctName), text = ifelse(missing(fctText), "Brain-Cousens (hormesis)", 
										 																			fctText), noParm = sum(is.na(fixed)))
	class(returnList) <- "linear_braincousens"
	invisible(returnList)
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
