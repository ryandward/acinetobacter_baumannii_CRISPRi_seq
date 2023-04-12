require(conflicted)
require(pacman)
source("drc_logistic_functions.R")

aba_key <- fread("aba_key.tsv")

p_load(
	"data.table",
	"tidyverse",
	"broom",
	"modelr")

p_load_current_gh(
	"DoseResponse/drcData",
	"ryandward/drc",
	"hrbrmstr/hrbrthemes")

conflicted::conflicts_prefer(
	gtools::permute,
	dplyr::filter,
	dplyr::select,
	drc::gaussian)

# https://jgeb.springeropen.com/articles/10.1186/s43141-020-00048-4

# all genes
interested.genes <- fread("curated_names.tsv") %>% pull(unique_name)

interested.conditions <- c(
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

# read results

melted_results <- fread("Results/melted_results.tsv.gz", sep = "\t")

median_melted_results <- fread("Results/median_melted_results.tsv.gz", sep = "\t")

curated_names <- fread("curated_names.tsv", sep = "\t")

aba_key <- fread("aba_key.tsv", sep = "\t")


melted_results$y_pred <- NULL

melted_results[aba_key, on = .(spacer), y_pred := y_pred]

drm.try <- possibly(drm, otherwise = NA)
augment.try <- possibly(augment, otherwise = NA)
glance.try <- possibly(glance, otherwise = NA)
tidy.try <- possibly(tidy, otherwise = NA)

#estimate dose-response models

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

################################################################################

test.genes <- c("murA", "rpmB", "aroC", "GO593_00515", "glnS", "nuoB", "lpxC")

LL.4.parameters <- c("shape", "min_value", "max_value", "kd_50")
BC.5.parameters <- c("shape", "min_value", "max_value", "kd_50", "hormesis")

################################################################################



L.4.model <- mismatches %>%
	filter(unique_name %in% test.genes) %>%
	mutate(fit = case_when(
		response.max > 0 ~ map2(
			data,
			response.max,
			~drm.try(
				data = .x, 
				LFC.adj ~ y_pred, 
				control = drmc(method = "Nelder-Mead", noMessage = TRUE, maxIt = 50000), 
				fct = L.4(fixed = c(NA, 0, .y, NA), names = LL.4.parameters)
			)
		),
		response.max < 0 ~ map2(
			data,
			response.max,
			~drm.try(
				data = .x, 
				LFC.adj ~ y_pred, 
				control = drmc(method = "Nelder-Mead", noMessage = TRUE, maxIt = 50000), 
				fct = L.4(fixed = c(NA, .y, 0, NA), names = LL.4.parameters)
			)
		)
	))

BC.5.model <- mismatches %>%
	filter(unique_name %in% test.genes) %>%
	mutate(fit = case_when(
		response.max > 0 ~ map2(
			data,
			response.max,
			~drm.try(
				data = .x, 
				LFC.adj ~ y_pred, 
				control = drmc(method = "Nelder-Mead", noMessage = TRUE, maxIt = 50000), 
				fct = BC.5.logistic(fixed = c(NA, 0, .y, NA, NA), names = BC.5.parameters)
			)
		),
		response.max < 0 ~ map2(
			data,
			response.max,
			~drm.try(
				data = .x, 
				LFC.adj ~ y_pred, 
				control = drmc(method = "Nelder-Mead", noMessage = TRUE, maxIt = 50000), 
				fct = BC.5.logistic(fixed = c(NA, .y, 0, NA, NA), names = BC.5.parameters)
			)
		)
	))

################################################################################

calculate_lrt <- function(this.gene, this.condition, this.HA, this.H0) {
	# Filter the BC.5.model and LL.4.model data frames for the given gene and condition
	this.HA <- BC.5.model %>%
		filter(unique_name == this.gene) %>%
		filter(condition == this.condition) %>%
		select(fit) %>%
		pull(fit)
	
	this.H0 <- LL.4.model %>%
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
	# Filter the BC.5.model and LL.4.model data frames for the given gene and condition
	this.HA <- BC.5.model %>%
		filter(unique_name == this.gene) %>%
		filter(condition == this.condition) %>%
		select(fit) %>%
		pull(fit)
	
	this.H0 <- LL.4.model %>%
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


############################