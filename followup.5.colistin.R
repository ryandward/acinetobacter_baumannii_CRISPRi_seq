library(pacman)

p_load(data.table, tidyverse, hrbrthemes)

p_load_current_gh("briandconnelly/growthcurve")

followup.colistin.2 <- fread(
	"Followups/followup.colistin.2.tsv",
	header = TRUE)

##########################################################################################
# wrangle the data from the sunrise machine... 

followup.colistin.2 <- followup.colistin.2[, -c(1,3)] %>% 
	melt(id.vars = "Time [s]", variable.name = "well", value.name = "OD600") %>% 
	rename(time = "Time [s]")

followup.colistin.2[, OD600 := as.numeric(OD600)]
followup.colistin.2[, time := as.numeric(time)]

##########################################################################################
# declare a pattern that corresponds to the wells that actually have cells in them

wells.filled <-
	CJ(well.row = toupper(letters[3:5]),
		 well.col = c(2:11))[
		 	, .(well = paste0(well.row, well.col))]

followup.colistin.2.strains <- c("control.2", "lpxC", "nuoH") %>% rep(each = 10) %>% rep(1) %>% 
	data.table(strain = .)

followup.colistin.2.reps <- rep(1, 30) %>%
	data.table(rep = .)

followup.colistin.2.induced <- c("off") %>% rep(each = 3*10) %>%
	data.table(induced = .)

followup.colistin.2.drug <- rep("colistin", 30) %>%
	data.table(drug = .)

followup.colistin.2.dose <- (100)/(2^(0:9)) %>% rep(3) %>% 
	data.table(dose = .)

followup.colistin.2.wells <- cbind(
	wells.filled, 
	followup.colistin.2.strains, 
	followup.colistin.2.reps, 
	followup.colistin.2.induced, 
	followup.colistin.2.drug, 
	followup.colistin.2.dose)

##########################################################################################
# combine the data about the identity of each well with the plate reader data
followup.colistin.2 <- followup.colistin.2 %>% full_join(followup.colistin.2.wells)
followup.colistin.2[, dose := factor(as.character(dose))]

# followup.colistin.2 <- followup.colistin.2 %>% filter(time <= (18*60*60))

followup.colistin.2.plot <- followup.colistin.2 %>% 
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
	facet_wrap(facets = c("drug", "dose"), ncol = 5) +
	ylim(0, NA) +
	scale_colour_brewer(palette = "Dark2") +
	scale_fill_brewer(palette = "Dark2") +
	ggtitle(
		bquote(bold("Colistin phenotypes for")~bolditalic("Acinetobacter baumannii."))) +
	theme_ipsum()

print(followup.colistin.2.plot)
