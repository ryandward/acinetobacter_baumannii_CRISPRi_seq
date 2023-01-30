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

doc_theme <- theme_ipsum(
	base_family = "Arial", 
	caption_margin = 12,
	axis_title_size = 12,
	axis_col = "black")


conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

competition_followups <- fread("competition_followups.tsv")

drug_sets <- competition_followups %>% select(set, antibiotic) %>% filter(antibiotic != "none") %>% unique

competition_deltas <- competition_followups %>% 
	mutate(
		treatment = case_when(antibiotic == "none" ~ "none", TRUE ~ "drug"),
		treatment = factor(treatment, levels = c("none", "drug")),
		`unique name` = case_when(`unique name` =="" ~ "control", TRUE ~ `unique name`)) %>%
	dcast(set + strains + `unique name` ~ timing + treatment + replicate, value.var = "proportion") %>%
	mutate(
		T0_none_A_d = (T0_none_A/T0_none_A),
		T0_none_B_d = (T0_none_B/T0_none_B),
		T1_none_A_d = (T1_none_A/T0_none_A),
		T1_none_B_d = (T1_none_B/T0_none_B),
		T1_drug_A_d = (T1_drug_A/T0_none_A),
		T1_drug_B_d = (T1_drug_B/T0_none_B),
		T2_none_A_d = (T2_none_A/T0_none_A),
		T2_none_B_d = (T2_none_B/T0_none_B),
		T2_drug_A_d = (T2_drug_A/T0_none_A),
		T2_drug_B_d = (T2_drug_B/T0_none_B))

delta_columns <- competition_deltas %>% 
	colnames %>% `[`(competition_deltas %>% colnames %>% grepl("_d$", .))

competition_deltas <- competition_deltas %>% 
	select(set, strains, `unique name`, all_of(delta_columns)) %>% 
	melt(id.vars = c("set", "strains", "unique name")) %>% 
	mutate(variable = gsub("_d$", "", variable)) %>%
	separate(variable, c("timing", "treatment", "replicate")) %>%
	mutate(treatment = factor(treatment, levels = c("none", "drug"))) 

competition_deltas %>% inner_join(drug_sets) %>% 
	mutate(
		time_drug = paste(timing, treatment, sep = "\n"),
		time_drug = case_when(time_drug %like% "T0" ~ "T0", TRUE ~ time_drug),
		time_drug = factor(time_drug, levels = c("T0", "T1\nnone", "T1\ndrug", "T2\nnone", "T2\ndrug"))) %>% 
	ggplot(aes(y = value, x = time_drug, fill = `unique name`)) + 
	geom_bar(position = "dodge", stat = "identity") + 
	facet_grid( ~ set + antibiotic ,  scales = "free", space = "free") + doc_theme +
	scale_fill_manual(
		values = c(
			"control" = "grey", 
			"nuoB" = "#6A3D9A",
			"lpxC" = "#33A02C",
			"murA" = "#FF7F00",
			"glnS" = "#1F78B4"))

