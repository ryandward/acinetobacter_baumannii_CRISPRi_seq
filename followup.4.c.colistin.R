library(pacman)

p_load(data.table, tidyverse, hrbrthemes)

p_load_current_gh("briandconnelly/growthcurve")

followup.colistin <- fread(
	"Followups/followup.colistin.tsv",
	header = TRUE)

##########################################################################################
# wrangle the data from the sunrise machine... 

followup.colistin <- followup.colistin %>% 
	melt(id.vars = c("time", "row"), variable.name = "column", value.name = "OD600") %>% 
	arrange(time, row, column) %>%
	unite("well", c(row, column), sep = "")

followup.colistin[, OD600 := as.numeric(OD600)]
followup.colistin[, time := as.numeric(time)]

##########################################################################################
# declare a pattern that corresponds to the wells that actually have cells in them

wells.filled <-
	CJ(well.row = toupper(letters[2:7]),
		 well.col = c(1:9))[
		 	, .(well = paste0(well.row, well.col))]

followup.colistin.strains <- c("control.2", "lpxC", "nuoH") %>% rep(each = 9) %>% rep(2) %>% 
	data.table(strain = .)

followup.colistin.reps <- rep(1, 54) %>%
	data.table(rep = .)

followup.colistin.induced <- c("off", "off") %>% rep(each = 3*9) %>%
	data.table(induced = .)

followup.colistin.drug <- rep("colistin", 54) %>%
	data.table(drug = .)

followup.colistin.dose <- (10000/10/2)/(2^(0:8)) %>% rep(6) %>% 
	data.table(dose = .)

followup.colistin.wells <- cbind(
	wells.filled, 
	followup.colistin.strains, 
	followup.colistin.reps, 
	followup.colistin.induced, 
	followup.colistin.drug, 
	followup.colistin.dose)

##########################################################################################
# combine the data about the identity of each well with the plate reader data
followup.colistin <- followup.colistin %>% full_join(followup.colistin.wells)
followup.colistin[, dose := factor(as.character(dose))]

# followup.colistin <- followup.colistin %>% filter(time <= (18*60*60))

followup.colistin.plot <- followup.colistin %>% 
	mutate(hour = time/60/60) %>%
	mutate(dose = as.numeric(levels(dose))[dose]) %>% 
	arrange(desc(dose)) %>% 
	mutate(dose = round(dose, 3)) %>%
	mutate(dose = factor(dose, levels = unique(dose))) %>%
	filter(!is.na(strain)) %>%
	ggplot(aes(x = hour, y = OD600, fill = strain, colour = strain)) + 
	stat_smooth(
		fullrange = TRUE, 
		level = 0.99999,
		method = "gam",
		formula = y ~ s(x, bs = "cs")) +
	facet_wrap(facets = c("drug", "dose"), ncol = 3) +
	scale_colour_brewer(palette = "Dark2") +
	scale_fill_brewer(palette = "Dark2") +
	ggtitle(
		bquote(bold("Colistin phenotypes for")~bolditalic("Acinetobacter baumannii."))) +
	theme_ipsum()

print(followup.colistin.plot)