library(pacman)

p_load(data.table, tidyverse, hrbrthemes)

followup.2 <- fread(
	"Followups/followup2.tsv",
	header = TRUE)

##########################################################################################
# wrangle the data from the sunrise machine... 

setnames(followup.2, t(followup.2)[,1])
followup.2 <- followup.2[-1:-2,]
followup.2 <- followup.2 %>% rename(well = `Time [s]`)


followup.2 <- 
	melt(
		followup.2, 
		id.vars = "well", 
		variable.name = "time", 
		value.name = "OD600", 
		na.rm = TRUE)

followup.2[, OD600 := as.numeric(OD600)]
followup.2[, time := as.numeric(levels(time))[time]]

##########################################################################################
# declare a pattern that corresponds to the wells that actually have cells in them

wells.filled <-
	CJ(well.row = toupper(letters[2:7]),
		 well.col = c(1:12))[
		 	, .(well = paste0(well.row, well.col))]

followup.2.strains <- c("control.2", "nuoF", "nuoB") %>% rep(each = 12) %>% rep(2) %>% 
	data.table(strain = .)

followup.2.reps <- 1:3 %>% rep(4*6) %>%
	data.table(rep = .)

followup.2.induced <- c("off","on") %>% rep(each = 3*12) %>%
	data.table(induced = .)

followup.2.drug <- c("colistin", "rifampicin") %>% rep(each = 6) %>% rep(6) %>%
	data.table(drug = .)

followup.2.dose <- c(0,0.44, 0.34, 0) %>% rep(each = 3) %>% rep(6) %>% 
	data.table(dose = .)

followup.2.wells <- cbind(
	wells.filled, 
	followup.2.strains, 
	followup.2.reps, 
	followup.2.induced, 
	followup.2.drug, 
	followup.2.dose)

##########################################################################################
# combine the data about the identity of each well with the plate reader data
followup.2 <- followup.2 %>% full_join(followup.2.wells)
followup.2[, dose := factor(as.character(dose))]

followup.2 <- followup.2 %>% filter(time <= (18*60*60))

followup.2.plot <- 
	followup.2 %>% 
	mutate(hour = time/60/60) %>%
	filter(!is.na(strain)) %>%
	ggplot(aes(x = hour, y = OD600, fill = strain, colour = strain)) + 
	stat_smooth(
		fullrange = TRUE, 
		level = 0.999,
		method = "gam",
		formula = y ~ s(x, bs = "cs")) +
	facet_wrap(facets = c("drug", "dose", "induced"), ncol = 4) +
	scale_colour_brewer(palette = "Dark2") +
	scale_fill_brewer(palette = "Dark2") +
	ggtitle(bquote(bold("Growth Phenotypes for")~bolditalic("Acinetobacter baumannii"))) +
	theme_ipsum()

print(followup.2.plot)

################################################################################