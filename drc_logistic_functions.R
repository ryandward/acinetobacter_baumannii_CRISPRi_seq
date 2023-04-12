# Brain-Cousens Hormetic DRC Logistic Model (Not LLogistic)
# Ryan Ward
# Tuesday, April 11, 2023

pacman::p_load(lmtest)

L.4.parameters <- c("hill", "min_value", "max_value", "kd_50")
BC.5.parameters <- c("shape", "min_value", "max_value", "kd_50", "hormesis")

# BC.5 Logistic function model using L.4.parameters and BC.5.parameters
BC.5.logistic <-
	function(fixed = c(NA, NA, NA, NA, NA),
						names = c("b", "c", "d", "e", "f"),
						...)
	{
		numParm <- 5
		if (!is.character(names) | !(length(names) == numParm)) {
			stop("Not correct 'names' argument")
		}
		if (!(length(fixed) == numParm)) {
			stop("Not correct length of 'fixed' argument")
		}
		return(braincousens_logistic(
			names = names,
			fixed = fixed,
			fctName = as.character(match.call()[[1]]),
			...
		))
	}

# Brain-Cousens logistic model for hormesis with a linear x
braincousens_logistic <- function(
		fixed = c(NA, NA, NA, NA, NA),
		names = c("b", "c", "d", "e", "f"),
		method = c("1", "2", "3", "4"),
		ssfct = NULL,
		fctName,
		fctText) {
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
		parmMat[, 2] + (parmMat[, 3] + parmMat[, 5] * dose - parmMat[, 2]) / (1 + exp(parmMat[, 1] * (dose - parmMat[, 4])))
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
			logitTrans <- log((startVal[3] - resp3) / (resp3 -
																								 	startVal[2] + 0.001))
			logitFit <- lm(logitTrans ~ log(dose3))
			startVal[4] <- exp((-coef(logitFit)[1] / coef(logitFit)[2]))
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
		t4 <- (1 + t2) ^ (-2)
		cbind(-t1 * t2 * t4, 1 - 1 / t3, 1 / t3, t1 * t2 * parmMat[, 1] * t4, dose /
						t3)[, notFixed]
	}
	deriv2 <- NULL
	edfct <- function(
		parm,
		respl,
		reference,
		type,
		lower = 0.001,
		upper = 1000,
		...) {
		interval <- c(lower, upper)
		parmVec[notFixed] <- parm
		p <- EDhelper(parmVec, respl, reference, type)
		tempVal <- (100 - p) / 100
		helpEqn <- function(dose) {
			expVal <- exp(parmVec[1] * (log(dose) - log(parmVec[4])))
			parmVec[5] * (1 + expVal * (1 - parmVec[1])) - 
				(parmVec[3] - parmVec[2]) * expVal * parmVec[1] / dose
		}
		maxAt <- uniroot(helpEqn, interval)$root
		eqn <- function(dose) {
			tempVal * (1 + exp(parmVec[1] * (log(dose) - log(parmVec[4])))) -
				(1 + parmVec[5] * dose / (parmVec[3] - parmVec[2]))
		}
		EDp <- uniroot(eqn, lower = maxAt, upper = upper)$root
		EDdose <- EDp
		tempVal1 <- exp(parmVec[1] * (log(EDdose) - log(parmVec[4])))
		tempVal2 <- parmVec[3] - parmVec[2]
		derParm <-
			c(
				tempVal * tempVal1 * (log(EDdose) - log(parmVec[4])),
				-parmVec[5] * EDdose / ((tempVal2) ^ 2),
				parmVec[5] *
					EDdose / ((tempVal2) ^ 2),
				-tempVal * tempVal1 *
					parmVec[1] / parmVec[4],
				-EDdose / tempVal2
			)
		derDose <- tempVal * tempVal1 * parmVec[1] / EDdose -
			parmVec[5] / tempVal2
		EDder <- derParm / derDose
		return(list(EDp, EDder[notFixed]))
	}
	maxfct <- function(parm,
										 lower = 0.001,
										 upper = 1000) {
		parmVec[notFixed] <- parm
		if (parmVec[1] < 1) {
			stop("Brain-Cousens model with b<1 not meaningful")
		}
		if (parmVec[5] < 0) {
			stop("Brain-Cousens model with f<0 not meaningful")
		}
		optfct <- function(t) {
			expTerm1 <- parmVec[5] * t
			expTerm2 <- exp(parmVec[1] * (log(t) - log(parmVec[4])))
			return(parmVec[5] * (1 + expTerm2) - (parmVec[3] -
																							parmVec[2] + expTerm1) * expTerm2 * parmVec[1] /
						 	t)
		}
		ED1 <- edfct(parm, 1, lower, upper)[[1]]
		doseVec <- exp(seq(log(1e-06), log(ED1), length = 100))
		maxDose <- uniroot(optfct, c((doseVec[optfct(doseVec) >
																						0])[1], ED1))$root
		return(c(maxDose, fct(maxDose, matrix(
			parm, 1, length(names)
		))))
	}
	returnList <- list(
		fct = fct,
		ssfct = ssfct,
		names = names,
		deriv1 = deriv1,
		deriv2 = deriv2,
		edfct = edfct,
		maxfct = maxfct,
		name = ifelse(missing(fctName), as.character(match.call()[[1]]),
									fctName),
		text = ifelse(missing(fctText), "Brain-Cousens (hormesis)",
									fctText),
		noParm = sum(is.na(fixed))
	)
	class(returnList) <- "braincousens_logistic"
	invisible(returnList)
}

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

# Saves the results of the analysis, including vulnerability summary, predictions, data points, model performance, and model parameters
save_results <- function(
		vuln.summary, fit_predictions, fit_points, model_performance, 
		model_parameters, file_names, output_dir = "Results") {
	if (!dir.exists(output_dir)) {dir.create(output_dir)}
	
	fwrite(vuln.summary, file.path(output_dir, file_names$vuln_summary), sep = "\t")
	fwrite(fit_predictions, file.path(output_dir, file_names$fit_predictions))
	fwrite(fit_points, file.path(output_dir, file_names$fit_points))
	fwrite(model_performance, file.path(output_dir, file_names$model_performance), sep = "\t")
	fwrite(model_parameters, file.path(output_dir, file_names$model_parameters), sep = "\t")
}

compare_models <- function(full_model, reduced_model) {
	
	# Check if both models have the same unique combinations of gene and condition
	full_model_combinations <- full_model %>% select(unique_name, condition) %>% distinct()
	reduced_model_combinations <- reduced_model %>% select(unique_name, condition) %>% distinct()
	
	if (!identical(full_model_combinations, reduced_model_combinations)) {
		stop("The full and reduced models do not have the same unique combinations of genes and conditions.")
	}
	
	# Create a nested data frame for each unique combination of gene and condition
	all_gene_conditions <- full_model_combinations %>%
		nest_by(unique_name, condition)
	
	# Calculate ANOVA and LRT test results for all genes and conditions
	results <- all_gene_conditions %>%
		mutate(
			HA_fit = list(collect(select(full_model, unique_name, condition, fit))),
			H0_fit = list(collect(select(reduced_model, unique_name, condition, fit))),
			anova_result = purrr::map2(HA_fit, H0_fit, ~anova(.x[[1]]$fit[[1]], .y[[1]]$fit[[1]])),
			lrt_result = purrr::map2(HA_fit, H0_fit, ~lrtest(.x[[1]]$fit[[1]], .y[[1]]$fit[[1]]))
		) %>%
		select(-HA_fit, -H0_fit) %>%
		unnest(cols = c(anova_result, lrt_result))
	
	# Extract p-values from the test results and return them in a tidy format
	p_values <- results %>%
		mutate(
			anova_p_value = purrr::map_dbl(anova_result, ~.x[["Pr(>F)"]][2]),
			lrt_p_value = purrr::map_dbl(lrt_result, ~.x[["Pr(>Chisq)"]][1])
		) %>%
		select(unique_name, condition, anova_p_value, lrt_p_value)
	
	return(p_values)
}



##########################################################################################
# scratch pad

# Define LRT and ANOVA calculation functions
calculate_lrt <- function(this.gene, this.condition, this.HA, this.H0) {
	
	# Filter the BC.5.model and LBC.5.reduced data frames for the given gene and condition
	this.HA <- BC.5_model %>%
		filter(unique_name == this.gene) %>%
		filter(condition == this.condition) %>%
		select(fit) %>%
		pull(fit)
	
	this.H0 <- BC.5_reduced_model %>%
		filter(unique_name == this.gene) %>%
		filter(condition == this.condition) %>%
		select(fit) %>%
		pull(fit)
	
	# Perform the likelihood ratio test
	lrt.result <- lrtest(this.HA[[1]], this.H0[[1]])
	
	# Return the likelihood ratio test result
	# return(lrt.result)
	return(lrt.result)
}

calculate_anova <- function(this.gene, this.condition, this.HA, this.H0) {
	
	# Filter the BC.5.model and LBC.5.reduced data frames for the given gene and condition
	this.HA <- BC.5_model %>%
		filter(unique_name == this.gene) %>%
		filter(condition == this.condition) %>%
		select(fit) %>%
		pull(fit)
	
	this.H0 <- BC.5_reduced_model %>%
		filter(unique_name == this.gene) %>%
		filter(condition == this.condition) %>%
		select(fit) %>%
		pull(fit)
	
	# Perform the likelihood ratio test
	anova.result <- anova(this.HA[[1]], this.H0[[1]])
	
	# Return the likelihood ratio test result
	# return(lrt.result) 
	return(anova.result)
}