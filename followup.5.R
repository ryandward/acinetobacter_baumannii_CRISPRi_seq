library(pacman)

p_load(data.table, tidyverse, hrbrthemes)

followup.5 <- fread(
	"Followups/followup.5.tsv",
	header = TRUE)

##########################################################################################
# wrangle the data from the sunrise machine... 

followup.5 <- followup.5[, -c(1,3)] %>% 
	melt(id.vars = "Time [s]", variable.name = "well", value.name = "OD600") %>% 
	rename(time = "Time [s]")

followup.5[, OD600 := as.numeric(OD600)]
followup.5[, time := as.numeric(time)]

##########################################################################################
# declare a pattern that corresponds to the wells that actually have cells in them

wells.filled <-
	CJ(well.row = toupper(letters[2:7]),
		 well.col = c(1:12))[
		 	, .(well = paste0(well.row, well.col))]

followup.5.strains <- c("nuoH", "nuoB", "nuoF", "lpxC", "lpxA", "control.2") %>% rep(2) %>% rep(6) %>% 
	data.table(strain = .)

followup.5.reps <- 4:6 %>% rep(each = 12) %>% rep(2) %>%
	data.table(rep = .)

followup.5.induced <- c(
	c(0.0, 0.1) %>% rep(each = 6) %>% rep(3), 
	c(0.5, 1.0) %>% rep(each = 6) %>% rep(3)) %>%
	data.table(induced = .)

followup.5.drug <- c("no drug", "no drug") %>% rep(each = 6) %>% rep(6) %>%
	data.table(drug = .)

followup.5.dose <- c(0, 0, 0, 0) %>% rep(each = 3) %>% rep(6) %>% 
	data.table(dose = .)

followup.5.wells <- cbind(
	wells.filled, 
	followup.5.strains, 
	followup.5.reps, 
	followup.5.induced, 
	followup.5.drug, 
	followup.5.dose)

##########################################################################################
# combine the data about the identity of each well with the plate reader data
followup.5 <- followup.5 %>% full_join(followup.5.wells)
followup.5[, dose := factor(as.character(dose))]

followup.5 <- followup.5 %>% filter(time <= (18*60*60))

followup.5.plot <- 
	followup.5 %>%
	mutate(hour = time/60/60) %>%
	filter(!is.na(strain)) %>%
	ggplot(aes(x = hour, y = OD600, fill = strain, colour = strain)) + 
	stat_smooth(
		fullrange = TRUE, 
		level = 0.99999,
		method = "gam",
		formula = y ~ s(x, bs = "cs")) +
	facet_wrap(facets = c("drug", "dose", "induced"), ncol = 4) +
	scale_colour_brewer(palette = "Dark2") +
	scale_fill_brewer(palette = "Dark2") +
	ggtitle(bquote(bold("Growth Phenotypes for")~bolditalic("Acinetobacter baumannii"))) +
	theme_ipsum()

print(followup.5.plot)

#########################################################################
followup.5.stat.plot <-
	followup.5[
		!is.na(strain), 
		SummarizeGrowth(
			time/60/60, 
			OD600,
			blank = followup.5[
				is.na(strain), 
				median(OD600), 
				by = .(time)]$V1)$vals, 
		by = .(well)] %>% 
	inner_join(followup.5) %>% 
	mutate(dose = case_when(dose == 0 ~ "no drug", dose != 0 ~ "with drug")) %>%
	ggplot(aes(x = dose, y = t_mid, fill = strain)) + 
	geom_boxplot(position = "dodge") +
	ggtitle(
		bquote(bold("Hours to reach half capacity")~bolditalic("Acinetobacter baumannii.")),
		subtitle = ~"Induced 18 hours before exposure to antibiotics.") +
	theme_ipsum()+
	scale_fill_brewer(palette = "Dark2") +
	ylim(0, NA) +
	facet_wrap(facets = c("drug", "induced"), ncol = 2)

print(followup.5.stat.plot)
