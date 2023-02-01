# Load necessary libraries
library(conflicted)
library(pacman)

# Use pacman to load the following packages
p_load(data.table, tidyverse, broom, modelr)
p_load_current_gh("DoseResponse/drcData", "ryandward/drc", "hrbrmstr/hrbrthemes")

# Define a custom theme for the ggplot visualization
doc_theme <- theme_ipsum(base_family = "Arial", caption_margin = 12, axis_title_size = 12, axis_col = "black")

# Prefer the dplyr version of select and filter functions in case of conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Read the data from "competition_followups.tsv"
competition_followups <- fread("competition_followups.tsv")

# Create a unique set of antibiotics from the competition follow-up data
drug_sets <- competition_followups %>% 
  select(set, antibiotic) %>% 
  filter(antibiotic != "none") %>% 
  unique

# Filter the raw competition data for "timing == "T0"", add a column "treatment" with the value "drug", 
# bind the filtered data to the original competition data with fill, join with the drug_sets, 
# and mutate the columns "timing" and "treatment" to factor variables
competition_raw %>%
	filter(timing == "T0") %>%
	mutate(treatment = "drug") %>%
	rbind(competition_raw, fill = TRUE) %>%
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
	stat_summary(
		data = . %>% filter(timing != "T0" | treatment != "drug"),
		fun.data = "mean_sdl",
		geom = "errorbar",
		position = position_dodge2(width = 1),
		lwd = 1,
		mapping = aes(colour = `unique name`, alpha = treatment)
	) +
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
			"drug" = 0.5,
			"none" = 1
		)
	) +
	doc_theme +
	facet_grid(~set + antibiotic)

