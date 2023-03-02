require(conflicted)
require(pacman)

p_load(
	"data.table", 
	"tidyverse", 
	"broom", 
	"modelr")

p_load_current_gh(
	"DoseResponse/drcData",
	"ryandward/drc",
	"hrbrmstr/hrbrthemes")

conflict_prefer("gaussian", "drc")

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflicted::conflicts_prefer(gtools::permute)

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
	"Imipenem_0.09_T2 - None_0_T0",
	"Colistin_0.44_T1 - None_0_T1",
	"Colistin_0.44_T2 - None_0_T2",
	"Rifampicin_0.34_T1 - None_0_T1",
	"Rifampicin_0.34_T2 - None_0_T2",
	"Meropenem_0.11_T1 - None_0_T1",
	"Meropenem_0.17_T1 - None_0_T1",
	"Meropenem_0.11_T2 - None_0_T2",
	"Meropenem_0.17_T2 - None_0_T2",
	"Imipenem_0.06_T1 - None_0_T1",
	"Imipenem_0.09_T1 - None_0_T1",
	"Imipenem_0.06_T2 - None_0_T2",
	"Imipenem_0.09_T2 - None_0_T2"
)

# read results

melted_results <- fread("Results/melted_results.tsv.gz", sep = "\t")

median_melted_results <- fread("Results/median_melted_results.tsv.gz", sep = "\t")

curated_names <- fread("curated_names.tsv", sep = "\t")

aba_key <- fread("aba_key.tsv", sep = "\t")

# function to rescale between 0 and 1
# range01 <- function(x){(x - min(x))/(max(x) - min(x))}

# melted_results[!is.na(y_pred), y_pred := range01(y_pred)]

# melted_results <- melted_results %>% group_by(target) %>% mutate(y_pred = range01(y_pred))

# melted_results <- tibble(melted_results)

melted_results$y_pred <- NULL

melted_results <- melted_results %>% inner_join(aba_key %>% select(y_pred, spacer)) 

# define parameters of interest for dose-response plots

drm.try <- possibly(drm, otherwise = NA)
augment.try <- possibly(augment, otherwise = NA)
glance.try <- possibly(glance, otherwise = NA)
tidy.try <- possibly(tidy, otherwise = NA)

#estimate dose-response models

mismatches <- melted_results %>%
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
			rename(response_max = LFC.adj)) %>%
	nest(data = c(-condition, -unique_name, -response_max))

################################################################################

L4.parameters <- c("hill", "min_value", "max_value", "kd_50")
BC5.parameters <- c("shape", "min_value", "max_value", "kd_50", "hormesis")

# In CRISPRi, either the min or max value should be constrained to 0.
zero_min <- c(NA, 0, NA, NA, NA)
zero_max <- c(NA, NA, 0, NA, NA)

# names(zero_min) <- names(zero_max) <- BC5.parameters
# 
# mismatches <- mismatches %>%
# 	mutate(fit = map(data, ~ drm.try(data = .x, LFC.adj ~ exp(y_pred), fct = BC.5(names = BC5.parameters))))


# # look for each gene whether max or min is closer to zero
# mismatches <- mismatches %>% 
# 	mutate(fit = map(
# 		data, ~ drm.try(
# 			data = .x, 
# 			LFC.adj ~ y_pred, 
# 			fct = L.4(names = L4.parameters))))
# 
# # which end varies the most away from zero?
# away_from_zero <- mismatches %>%
# 	mutate(
# 		results_look = map(fit, tidy)) %>%
# 	select(unique_name, condition, results_look) %>%
# 	unnest(results_look) %>%
# 	group_by(unique_name, condition) %>%
# 	filter(term %in% c("min_value", "max_value")) %>%
# 	filter(abs(estimate) == max(abs(estimate))) %>%
# 	select(term) %>%
# 	rename("variable_term" = "term")
# 
# mismatches <- mismatches %>% inner_join(away_from_zero)

mismatches <- mismatches %>% 
	# filter(unique_name == "lpxC" | unique_name == "nuoB") %>%
	mutate(fit = case_when(
		response_max > 0 ~ map2(data, response_max, ~ drm.try(data = .x, LFC.adj ~ y_pred, fct = BC.5(fixed = c(NA, 0, .y, NA, NA), names = BC5.parameters))),
		response_max < 0 ~ map2(data, response_max, ~ drm.try(data = .x, LFC.adj ~ y_pred, fct = BC.5(fixed = c(NA, .y, 0, NA, NA), names = BC5.parameters)))))


# mismatches <- mismatches %>% filter(unique_name == "lpxC" | unique_name == "nuoB") %>%
# 	mutate(fit = case_when(
# 		response_max > 0 ~ map2(data, response_max, ~ drm.try(data = .x, LFC.adj ~ y_pred, fct = L.4(fixed = c(NA, 0, .y, NA), names = L4.parameters))),
# 		response_max < 0 ~ map2(data, response_max, ~ drm.try(data = .x, LFC.adj ~ y_pred, fct = L.4(fixed = c(NA, .y, 0, NA), names = L4.parameters)))))

################################################################################

mismatches <- mismatches %>% 
	filter(!is.na(fit)) %>%
	mutate(p.vals = map(fit, tidy)) 

mismatches <- mismatches %>% 
	mutate(
		kd_50.tibble = map(
			p.vals,
			~filter(.x, term == 'kd_50') %>%
				select(p.value) %>%
				rename (., vuln.p = p.value))) %>%
	mutate(
		vuln.tibble = map2(
			fit,
			p.vals,
			~augment.try(
				.x,
				newdata = .y %>% filter (term == "kd_50") %>% select (estimate)) %>%
				rename (., vuln.est = .fitted, vuln.kd_50 = estimate))) %>%
	mutate(
		hill.tibble = map(
			p.vals,
			~filter(.x, term == 'hill') %>%
				select(estimate, p.value) %>%
				rename (., hill.est = estimate, hill.p = p.value)))
# 
mismatches <- mismatches %>%
	unnest(vuln.tibble) %>%
	unnest(kd_50.tibble) %>%
	select(-c(p.vals))

mismatches <- mismatches %>% 
	mutate(Gene = factor(unique_name, levels = unique(interested.genes))) %>%
	mutate(Condition = factor(condition, levels = unique(interested.conditions)))

# write and save a table that contains parameter estimates related to dose-response curves

vuln.summary <- mismatches %>% select(
	unique_name,
	condition,
	vuln.est,
	vuln.kd_50,
	vuln.p)

vuln.summary %>% fwrite("Results/hormetic_vulnerability_summary.tsv.gz", sep = "\t")

# add 90% confidence interval predictions to the dose-response curves

mismatches <- mismatches %>% 
	mutate(predictions = map2(fit, data, ~augment.try(
		.x,
		newdata = expand.grid(
			y_pred = seq(
				min(.y$y_pred),
				# 0,
				max(.y$y_pred),
				# 1,
				length = 250)),
		conf.int = T,
		conf.level = 0.90)))

# write and save fitted points, and fitted predictions

fit_predictions <- mismatches %>% 
	select(Gene, Condition, predictions) %>% 
	unnest(predictions)

fit_predictions %>% fwrite("Results/hormetic_fit_predictions.tsv.gz")

fit_points <- mismatches %>% 
	select(Gene, Condition, data) %>% 
	unnest(data)

fit_points %>% fwrite("Results/hormetic_fit_points.tsv.gz")

