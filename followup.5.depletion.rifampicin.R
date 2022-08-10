library(pacman)

p_load(data.table, tidyverse, hrbrthemes)


doc_theme <- theme_ipsum(
	base_family = "Arial", 
	caption_margin = 12,
	axis_title_size = 12,
	axis_col = "black")

followup.5.depletion.rifampicin <- fread(
	"Followups/followup.5.depletion.rifampicin.tsv",
	header = TRUE)

##########################################################################################
# wrangle the data from the sunrise machine... 


followup.5.depletion.rifampicin <- followup.5.depletion.rifampicin[-1:-2,]
followup.5.depletion.rifampicin <- followup.5.depletion.rifampicin %>% rename(well = `Cycle Nr.`)


followup.5.depletion.rifampicin <- 
	melt(
		followup.5.depletion.rifampicin, 
		id.vars = "well", 
		variable.name = "time", 
		value.name = "OD600", 
		na.rm = TRUE)


followup.5.depletion.rifampicin[, OD600 := as.numeric(OD600)]
followup.5.depletion.rifampicin[, time := as.numeric(levels(time))[time]]

followup.5.depletion.rifampicin[, time := (time - 1) * 300]

followup.5.depletion.rifampicin <- followup.5.depletion.rifampicin[OD600 != "NA"]

##########################################################################################
# declare a pattern that corresponds to the wells that actually have cells in them

wells.filled <-
	CJ(well.row = toupper(letters[1:8]),
		 well.col = c(1:12))[
		 	, .(well = paste0(well.row, well.col))]

followup.5.depletion.rifampicin.strains <- c("nuoB", "nuoF", "lpxC", "control.2") %>% rep(each = 12) %>% rep(96/length(.)) %>% 
	data.table(strain = .)

followup.5.depletion.rifampicin.reps <- 4:6 %>% rep(96/length(.)) %>%
	data.table(rep = .)

followup.5.depletion.rifampicin.induced <- c("off", "on") %>% rep(each = 96/length(.)) %>%
	data.table(induced = .)

followup.5.depletion.rifampicin.drug <- c("rifampicin") %>% rep(96) %>%
	data.table(drug = .)

followup.5.depletion.rifampicin.dose <- c(0, 0.96, 0.48, 0.24) %>% rep(each = 3) %>% rep(96/length(.)) %>% 
	data.table(dose = .)

followup.5.depletion.rifampicin.wells <- cbind(
	wells.filled, 
	followup.5.depletion.rifampicin.strains, 
	followup.5.depletion.rifampicin.reps, 
	followup.5.depletion.rifampicin.induced, 
	followup.5.depletion.rifampicin.drug, 
	followup.5.depletion.rifampicin.dose)

##########################################################################################
# combine the data about the identity of each well with the plate reader data
followup.5.depletion.rifampicin <- followup.5.depletion.rifampicin %>% full_join(followup.5.depletion.rifampicin.wells)
followup.5.depletion.rifampicin[, dose := factor(as.character(dose))]

followup.5.depletion.rifampicin <- followup.5.depletion.rifampicin %>% filter(time <= (18*60*60))

followup.5.depletion.rifampicin <- followup.5.depletion.rifampicin[,
																															 SummarizeGrowth(
																															 	time/60/60,
																															 	OD600)$vals,
																															 by = .(well)] %>% 
	inner_join(followup.5.depletion.rifampicin)

followup.5.depletion.rifampicin.plot <- 
	followup.5.depletion.rifampicin %>%
	filter(dose %in% c(0, 0.48, 0.96)) %>%
	mutate(dose = paste(dose, "ng/uL")) %>%
	filter(strain %in% c("control.2", "lpxC", "nuoB")) %>%
	mutate(hour = time/60/60) %>%
	filter(!is.na(strain)) %>%
	filter(note != "cannot fit data") %>%
	filter(induced == "on") %>%
	ggplot(aes(x = hour, y = OD600, fill = strain, colour = strain)) + 
	stat_smooth(
		alpha = 1,
		fullrange = TRUE, 
		level = 0.95,
		method = "gam",
		formula = y ~ s(x, bs = "cs")) +
	# geom_hline(yintercept = 0.5, linetype="dashed", color = "red") +
	facet_wrap(facets = c("dose"), ncol = 4) +
	scale_fill_ipsum()+
	scale_colour_ipsum() +
	# ggtitle(
	# 	bquote(bold("Growth Phenotypes for")~bolditalic("Acinetobacter baumannii.")),
	# 	subtitle = ~"Induced 18 hours before exposure to antibiotics.") +
	doc_theme

print(followup.5.depletion.rifampicin.plot)

#########################################################################

followup.5.depletion.rifampicin %>% 
	filter(dose != 0.96) %>% 
	ggplot(aes(x = dose, y = t_mid, fill = strain)) + 
	geom_boxplot(position = "dodge") +
	ggtitle(
		bquote(bold("Hours to reach half capacity")~bolditalic("Acinetobacter baumannii."))) +
	doc_theme +
	scale_fill_brewer(palette = "Dark2") +
	ylim(0, NA) +
	facet_wrap(facets = c("drug", "induced"), ncol = 2)
