library(pacman)

p_load(data.table, tidyverse, hrbrthemes)

followup.5.depletion.colistin <- fread(
	"Followups/followup.5.depletion.colistin.tsv",
	header = TRUE)

##########################################################################################
# wrangle the data from the sunrise machine... 

followup.5.depletion.colistin <- followup.5.depletion.colistin[, -c(2,3)] %>% 
	melt(id.vars = "Cycle Nr.", variable.name = "well", value.name = "OD600") %>% 
	rename(time = "Cycle Nr.") %>%
	mutate(time = (time - 1) * 300)

followup.5.depletion.colistin[, OD600 := as.numeric(OD600)]
followup.5.depletion.colistin[, time := as.numeric(time)]
followup.5.depletion.colistin <- followup.5.depletion.colistin[OD600 != "NA"]

##########################################################################################
# declare a pattern that corresponds to the wells that actually have cells in them

wells.filled <-
	CJ(well.row = toupper(letters[1:8]),
		 well.col = c(1:12))[
		 	, .(well = paste0(well.row, well.col))]

followup.5.depletion.colistin.strains <- c("nuoB", "nuoF", "lpxC", "control.2") %>% rep(each = 12) %>% rep(96/length(.)) %>% 
	data.table(strain = .)

followup.5.depletion.colistin.reps <- 4:6 %>% rep(96/length(.)) %>%
	data.table(rep = .)

followup.5.depletion.colistin.induced <- c("off", "on") %>% rep(each = 96/length(.)) %>%
	data.table(induced = .)

followup.5.depletion.colistin.drug <- c("colistin") %>% rep(96) %>%
	data.table(drug = .)

followup.5.depletion.colistin.dose <- c(0, 6, 4, 2) %>% rep(each = 3) %>% rep(96/length(.)) %>% 
	data.table(dose = .)

followup.5.depletion.colistin.wells <- cbind(
	wells.filled, 
	followup.5.depletion.colistin.strains, 
	followup.5.depletion.colistin.reps, 
	followup.5.depletion.colistin.induced, 
	followup.5.depletion.colistin.drug, 
	followup.5.depletion.colistin.dose)

##########################################################################################
# combine the data about the identity of each well with the plate reader data
followup.5.depletion.colistin <- followup.5.depletion.colistin %>% full_join(followup.5.depletion.colistin.wells)
followup.5.depletion.colistin[, dose := factor(as.character(dose))]

followup.5.depletion.colistin <- followup.5.depletion.colistin %>% filter(time <= (18*60*60))
followup.5.depletion.colistin <- followup.5.depletion.colistin[,
	SummarizeGrowth(
		time/60/60,
		OD600)$vals,
	by = .(well)] %>% 
	inner_join(followup.5.depletion.colistin)

followup.5.depletion.colistin.plot <- 
	followup.5.depletion.colistin %>%
	mutate(hour = time/60/60) %>%
	filter(!is.na(strain)) %>%
	filter(note != "cannot fit data") %>%
	ggplot(aes(x = hour, y = OD600, fill = strain, colour = strain)) + 
	stat_smooth(
		fullrange = TRUE, 
		level = 0.99999,
		method = "gam",
		formula = y ~ s(x, bs = "cs")) +
	geom_hline(yintercept = 0.5, linetype="dashed", color = "red") +
	facet_wrap(facets = c("drug", "dose", "induced"), ncol = 4) +
	scale_colour_brewer(palette = "Dark2") +
	scale_fill_brewer(palette = "Dark2") +
	ggtitle(
		bquote(bold("Growth Phenotypes for")~bolditalic("Acinetobacter baumannii.")),
		subtitle = ~"Induced 18 hours before exposure to antibiotics.") +
	theme_ipsum()

print(followup.5.depletion.colistin.plot)

#########################################################################
# 
# followup.5.depletion.colistin %>% 
# 	ggplot(aes(x = dose, y = t_mid, fill = strain)) + 
# 	geom_boxplot(position = "dodge") +
# 	ggtitle(
# 		bquote(bold("Hours to reach half capacity")~bolditalic("Acinetobacter baumannii."))) +
# 	theme_ipsum()+
# 	scale_fill_brewer(palette = "Dark2") +
# 	ylim(0, NA) +
# 	facet_wrap(facets = c("drug", "induced"), ncol = 2) %>%
# 	print
# 
