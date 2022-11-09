library(pacman)

p_load(data.table, tidyverse, hrbrthemes, growthcurver)


doc_theme <- theme_ipsum(
	base_family = "Arial", 
	caption_margin = 12,
	axis_title_size = 12,
	axis_col = "black")

followup.7.OD600 <- fread(
	"Followups/followup.7.OD600.tsv",
	header = TRUE)

followup.7.ThT <- fread(
	"Followups/followup.7.ThT.tsv",
	header = TRUE)

##########################################################################################
# wrangle the data from the sunrise machine... 

followup.7.OD600 <- followup.7.OD600[ -c(1,2) ,] %>% 
	rename(well = "Cycle Nr.") %>%
	melt(id.vars = "well", variable.name = "time", value.name = "OD600") %>% 
	mutate(time = (as.numeric(time) - 1) * 300)

followup.7.OD600[, OD600 := as.numeric(OD600)]
followup.7.OD600[, time := as.numeric(time)]
followup.7.OD600 <- followup.7.OD600[OD600 != "NA"]

######################################################

# followup.7.ThT <- followup.7.ThT[ -c(1,2) ,] %>% 
# 	rename(well = "Cycle Nr.") %>%
# 	melt(id.vars = "well", variable.name = "time", value.name = "ThT") %>% 
# 	mutate(time = (as.numeric(time) - 1) * 300)
# 
# min_over <- followup.7.ThT %>% filter(ThT == "OVER") %>% select(time) %>% min
# 
# followup.7.ThT <- followup.7.ThT %>% filter(time < min_over)
# 
# followup.7.ThT[, ThT := as.numeric(ThT)]
# followup.7.ThT[, time := as.numeric(time)]
# followup.7.ThT <- followup.7.ThT[ThT != "NA"]
# 
# followup.7.all <- followup.7.OD600 %>% 
# 	inner_join(followup.7.ThT)

followup.7.all <- followup.7.OD600

##########################################################################################
# declare a pattern that corresponds to the wells that actually have cells in them

wells.filled <-
	CJ(well.row = toupper(letters[1:8]),
		 well.col = c(1:12))[
		 	, .(well = paste0(well.row, well.col))]

followup.7.all.strains <- c("nuoB", "nuoF", "lpxC", "control.2") %>% rep(each = 12) %>% rep(96/length(.)) %>% 
	data.table(strain = .)

followup.7.all.reps <- 1:3 %>% rep(96/length(.)) %>%
	data.table(rep = .)

followup.7.all.induced <- c("on", "on") %>% rep(each = 96/length(.)) %>%
	data.table(induced = .)

followup.7.all.CCCP <- c("no CCCP", "CCCP") %>% rep(each = 96/length(.)) %>%
	data.table(CCCP = .)

followup.7.all.drug <- c("colistin") %>% rep(96) %>%
	data.table(drug = .)

followup.7.all.dose <- c(0, 6, 4, 2) %>% rep(each = 3) %>% rep(96/length(.)) %>% 
	data.table(dose = .)

followup.7.all.wells <- cbind(
	wells.filled, 
	followup.7.all.strains, 
	followup.7.all.reps, 
	followup.7.all.CCCP, 
	followup.7.all.drug, 
	followup.7.all.dose)

##########################################################################################
# combine the data about the identity of each well with the plate reader data

followup.7.all <- followup.7.all %>% full_join(followup.7.all.wells)
followup.7.all[, dose := factor(as.character(dose))]
##########################################################################################

followup.7.all <- followup.7.all[, SummarizeGrowth(
	time/60/60,
	OD600)$vals,
	by = .(well)] %>%
	inner_join(followup.7.all)

##########################################################################################
# 
followup.7.all.plot <-
	followup.7.all %>%
	mutate(hour = time/60/60) %>%
	filter(!is.na(strain)) %>%
	filter(hour <= 15) %>%
	ggplot(aes(x = hour, y = OD600, fill = strain, colour = strain)) +
	stat_smooth(
		fullrange = TRUE,
		level = 0.99999,
		method = "gam",
		formula = y ~ s(x, bs = "cs")) +
	# geom_hline(yintercept = 0.5, linetype="dashed", color = "red") +
	facet_wrap(facets = c("drug", "CCCP", "dose"), ncol = 4) +
	scale_colour_brewer(palette = "Dark2") +
	scale_fill_brewer(palette = "Dark2") +
	ggtitle(
		bquote(bold("Growth Phenotypes for")~bolditalic("Acinetobacter baumannii.")),
		subtitle = "Induced 18 hours before exposure to antibiotics.") +
	doc_theme

print(followup.7.all.plot)
# 
# 
# #########################################################################
# 
# followup.7.all.plot <- 
# 	followup.7.all %>%
# 	mutate(hour = time/60/60) %>%
# 	filter(!is.na(strain)) %>%
# 	filter(CCCP == "yes") %>%
# 	ggplot(aes(x = hour, y = ThT, fill = strain, colour = strain)) + 
# 	stat_smooth(
# 		fullrange = TRUE, 
# 		level = 0.99999,
# 		method = "gam",
# 		formula = y ~ s(x, bs = "cs")) +
# 	# geom_hline(yintercept = 0.5, linetype="dashed", color = "red") +
# 	facet_wrap(facets = c("drug", "dose"), ncol = 4) +
# 	scale_colour_brewer(palette = "Dark2") +
# 	scale_fill_brewer(palette = "Dark2") +
# 	ggtitle(
# 		bquote(bold("Raw ThT Fluorescence for")~bolditalic("Acinetobacter baumannii.")),
# 		subtitle = "Induced 18 hours before exposure to antibiotics.") +
# 	doc_theme
# 
# print(followup.7.all.plot)
# 
# #########################################################################
# 
# followup.7.all.plot <- 
# 	followup.7.all %>%
# 	mutate(hour = time/60/60) %>%
# 	filter(!is.na(strain)) %>%
# 	filter(CCCP == "yes") %>%
# 	ggplot(aes(x = hour, y = ThT/OD600, fill = strain, colour = strain)) + 
# 	stat_smooth(
# 		fullrange = TRUE, 
# 		level = 0.99999,
# 		method = "gam",
# 		formula = y ~ s(x, bs = "cs")) +
# 	# geom_hline(yintercept = 0.5, linetype="dashed", color = "red") +
# 	facet_wrap(facets = c("drug", "dose"), ncol = 4) +
# 	scale_colour_brewer(palette = "Dark2") +
# 	scale_fill_brewer(palette = "Dark2") +
# 	ggtitle(
# 		bquote(bold("ThT Fluorescence/OD600 for")~bolditalic("Acinetobacter baumannii.")),
# 		subtitle = "Induced 18 hours before exposure to antibiotics.") +
# 	doc_theme
# 
# print(followup.7.all.plot)


#########################################################################
# 
# followup.7.all %>% 
# 	ggplot(aes(x = dose, y = t_mid, fill = strain)) + 
# 	geom_boxplot(position = "dodge") +
# 	ggtitle(
# 		bquote(bold("Hours to reach half capacity")~bolditalic("Acinetobacter baumannii."))) +
# 	doc_theme +
# 	scale_fill_brewer(palette = "Dark2") +
# 	ylim(0, NA) +
# 	facet_wrap(facets = c("drug", "induced"), ncol = 2) %>%
# 	print
# 
