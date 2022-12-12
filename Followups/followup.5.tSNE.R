source("Followups/followup.5.depletion.R")
source("Followups/followup.5.depletion.colistin.R")
source("Followups/followup.5.depletion.rifampicin.R")

depletion_experiments <-
	followup.5.depletion %>% mutate(drug = "no drug") %>%
	rbind(followup.5.depletion.rifampicin %>% mutate(drug = "rifampicin"), fill = T) %>% 
	rbind(followup.5.depletion.colistin %>% mutate(drug = "colistin"), fill = T) %>%
	filter(time <= 15 * time * 60 * 60) %>%
	select(drug, well, OD600, time) %>%
	pivot_wider(id_cols = c(drug, well), names_from = time, values_from = OD600)

depletion_experiments_tSNE <- depletion_experiments %>% 
	unite(drug_well, c("drug", "well")) %>% 
	column_to_rownames("drug_well") %>%
	Rtsne(perplexity = 50)

depletion_experiments_tSNE <-	depletion_experiments_tSNE$Y %>% 
	data.table

depletion_experiments_tSNE <- depletion_experiments %>% 
	select(drug, well) %>% 
	cbind(depletion_experiments_tSNE) %>% 
	inner_join(
		followup.5.depletion.colistin.wells %>% 
			rbind(followup.5.depletion.rifampicin.wells) %>% 
			rbind(followup.5.depletion.wells)) %>%
	mutate(
		induced = case_when(
			induced == "off" ~ "0",
			induced == "on" ~ "0.5",
			TRUE ~ induced)) %>%
	mutate(induced = as.numeric(induced)) %>%
	mutate(rep = as.character(rep)) %>%
	mutate(induced = paste(induced, "IPTG"))



depletion_experiments_tSNE %>% 
	filter(drug == "colistin") %>% 
	ggplot(aes(x = V1, y = V2)) + 
	geom_point(aes(fill = strain, shape = rep), size = 5, alpha = 0.5) + 
	facet_grid(induced ~ dose + drug) + 
	scale_shape_manual(values = c(21, 22, 24)) +
	guides(fill = guide_legend(
		override.aes = list(shape = 23))) +
	doc_theme

depletion_experiments_tSNE %>% 
	filter(drug == "rifampicin") %>% 
	ggplot(aes(x = V1, y = V2)) + 
	geom_point(aes(fill = strain, shape = rep), size = 5, alpha = 0.5) + 
	facet_grid(induced ~ dose + drug) + 
	scale_shape_manual(values = c(21, 22, 24)) +
	guides(fill = guide_legend(
		override.aes = list(shape = 23))) +
	doc_theme

depletion_experiments_tSNE %>% 
	filter(drug == "no drug") %>% 
	ggplot(aes(x = V1, y = V2)) + 
	geom_point(aes(fill = strain, shape = rep), size = 5, alpha = 0.5) + 
	facet_grid(induced ~ dose + drug) + 
	scale_shape_manual(values = c(21, 22, 24)) +
	guides(fill = guide_legend(
		override.aes = list(shape = 23))) +
	doc_theme

