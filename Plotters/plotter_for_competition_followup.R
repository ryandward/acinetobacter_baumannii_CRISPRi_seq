# Load necessary libraries
library(conflicted)
library(pacman)

# Use pacman to load the following packages
p_load(data.table, tidyverse, broom, modelr, Hmisc)
p_load_current_gh("DoseResponse/drcData", "ryandward/drc", "hrbrmstr/hrbrthemes")

# Conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("summarize", "dplyr")

# Define a custom theme for the ggplot visualization
doc_theme <- theme_ipsum(base_family = "Arial", caption_margin = 12, axis_title_size = 12, axis_col = "black")

# Read the data from "competition_followups.tsv"
competition_followups <- fread("competition_followups.tsv")

# Create a unique set of antibiotics from the competition follow-up data
drug_sets <- competition_followups %>% 
  select(set, antibiotic) %>% 
  filter(antibiotic != "none") %>% 
  unique

# Create additional binary variable to indicate whether a strain has been treated with any drug, 
# turn empty named strains into "control", then dcast into wide-form data
competition_raw <- competition_followups %>% 
	mutate(
		treatment = case_when(antibiotic == "none" ~ "none", TRUE ~ "drug"),
		treatment = factor(treatment, levels = c("none", "drug")),
		`unique name` = case_when(`unique name` == "" ~ "control", TRUE ~ `unique name`)) %>%
	dcast(set + strains + `unique name` ~ timing + treatment + replicate, value.var = "proportion")

# Define list of data-filled columns, which in our case begin with T for T0/T1/T2
data_columns <- competition_raw %>% 
	colnames %>% `[`(competition_raw %>% colnames %>% grepl("^T", .))

# Beautify data a little more 
competition_raw <- competition_raw %>% 
	select(set, strains, `unique name`, all_of(data_columns)) %>% 
	melt(id.vars = c("set", "strains", "unique name")) %>% 
	mutate(variable = gsub("_d$", "", variable)) %>%
	separate(variable, c("timing", "treatment", "replicate")) %>%
	mutate(treatment = factor(treatment, levels = c("none", "drug"))) 

# Filter the raw competition data for "timing == "T0"", add a column "treatment" with the value "drug", 
# bind the filtered data to the original competition data with fill, join with the drug_sets, 
# and mutate the columns "timing" and "treatment" to factor variables
competition_raw %>%
	filter(timing == "T0") %>%
	mutate(treatment = "drug") %>%
	rbind(competition_raw) %>%
	inner_join(drug_sets) %>%
	mutate(
		timing = factor(timing, levels = c("T0", "T1", "T2")),
		treatment = factor(treatment, c("none", "drug"))
	) %>%
	filter(`unique name` != "control") %>%
	ggplot(aes(x = timing, y = value, group = interaction(`unique name`, treatment))) +
	stat_summary(
		aes(linetype = treatment, alpha = treatment),
		fun.data = "mean_sdl",
		geom = "line",
		lwd = 0.75
	) +
	geom_point(
		data = . %>% filter(timing != "T0" | treatment != "drug"),
		aes(colour = `unique name`, alpha = treatment),
		size = 2.5,
		position = position_dodge(width = 0.15)
	) +
	# stat_summary(
	# 	data = . %>% filter(timing != "T0" | treatment != "drug"),
	# 	fun.data = "mean_sdl",
	# 	geom = "errorbar",
	# 	position = position_dodge2(width = 1),
	# 	lwd = 1,
	# 	mapping = aes(colour = `unique name`, alpha = treatment)
	# ) +
	scale_colour_manual(
		values = c(
			"control" = "dark grey",
			"nuoB" = "#6A3D9A",
			"lpxC" = "#33A02C",
			"murA" = "#FF7F00",
			"glnS" = "#1F78B4"
		)
	) +
	scale_linetype_manual(
		values = c(
			"none" = "solid",
			"drug" = "twodash"
		)
	) +
	scale_alpha_manual(
		values = c(
			"drug" = 1,
			"none" = 0.25
		)
	) +
	doc_theme +
	facet_grid(~set + antibiotic)

