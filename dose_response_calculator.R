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

# https://jgeb.springeropen.com/articles/10.1186/s43141-020-00048-4

# interested.genes <- fread("../../genes_in_pathways.tsv") %>% pull(unique_name)
# interested.genes <- interested.genes %>% c("GO593_00515")

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
	"Imipenem_0.09_T2 - None_0_T2")

# read results

melted_results <- fread(
	"Results/melted_results.tsv.gz", sep = "\t")

curated_names <- fread(
	"curated_names.tsv", sep = "\t")

aba_key <- fread(
	"aba_key.tsv", sep = "\t")

# function to rescale between 0 and 1
# range01 <- function(x){(x - min(x))/(max(x) - min(x))}
# 
# melted_results[!is.na(y_pred), y_pred := range01(y_pred)]

melted_results$y_pred <- NULL

melted_results <- melted_results %>% inner_join(aba_key %>% select(y_pred, spacer)) 

melted_results <- tibble(melted_results)

# define parameters of interest for dose-response plots

L.4.parameters <- c("hill", "min_value", "max_value", "kd_50")


drm.try <- possibly(drm, otherwise = NA)

augment.try <- possibly(augment, otherwise = NA)

#estimate dose-response models

mismatches <- melted_results %>%
	filter(condition %in% interested.conditions) %>%
	filter(unique_name %in% interested.genes) %>%
	filter(type == "mismatch") %>%
	nest(data = c(-condition, -unique_name))

mismatches <- mismatches %>% mutate(fit = map(data, ~ drm.try(data = .x, LFC.adj ~y_pred, fct = L.4(names = L.4.parameters))))

mismatches <- mismatches %>% mutate(parameters = map(fit, tidy)) 

mismatches <- mismatches %>% 
	inner_join(
		mismatches %>% 
			select(unique_name, condition, parameters) %>% 
			unnest(parameters) %>% filter(term == "kd_50") %>% 
			rename(vuln.est = estimate, vuln.p = p.value) %>% 
			select(unique_name, condition, vuln.est, vuln.p))

mismatches <- mismatches %>% 
	mutate(vulnerability = map2(fit, vuln.est, ~augment(.x, newdata = data.frame(vuln.est = .y)))) %>%
	mutate(vuln.kd_50 = map_dbl(vulnerability, ".fitted")) %>% 
	select(unique_name, condition, vuln.est, vuln.kd_50, vuln.p, fit, data)

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

vuln.summary %>% fwrite("Results/vulnerability_summary.tsv.gz", sep = "\t")

# add 90% confidence interval predictions to the dose-response curves

mismatches <- mismatches %>% 
	mutate(predictions = map2(fit, data, ~augment.try(
		.x,
		newdata = expand.grid(
			y_pred = seq(
				min(.y$y_pred), 
				max(.y$y_pred), 
				length = 250)),
		conf.int = T,
		conf.level = 0.90)))

# write and save fitted points, and fitted predictions

fit_predictions <- mismatches %>% 
	select(Gene, Condition, predictions) %>% 
	unnest(predictions)

fit_predictions %>% fwrite("Results/fit_predictions.tsv.gz")

fit_points <- mismatches %>% 
	select(Gene, Condition, data) %>% 
	unnest(data)

fit_points %>% fwrite("Results/fit_points.tsv.gz")



