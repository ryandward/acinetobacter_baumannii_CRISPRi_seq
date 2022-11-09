library(pacman)

p_load(data.table, tidyverse, hrbrthemes)

followup.3 <- fread(
	"Followups/followup3.tsv",
	header = TRUE)

################################################################################
# wrangle the data from the sunrise machine... 

setnames(followup.3, t(followup.3)[,1])
followup.3 <- followup.3[-1:-2,]
followup.3 <- followup.3 %>% rename(well = `Time [s]`)


followup.3 <- 
	melt(
		followup.3, 
		id.vars = "well", 
		variable.name = "time", 
		value.name = "OD600", 
		na.rm = TRUE)

followup.3[, OD600 := as.numeric(OD600)]
followup.3[, time := as.numeric(levels(time))[time]]

################################################################################
# declare a pattern that corresponds to the wells that actually have cells in them

wells.filled <-
	CJ(well.row = toupper(letters[2:7]),
		 well.col = c(1:12))[
		 	, .(well = paste0(well.row, well.col))]

followup.3.strains <- c("control.2", "nuoH", "lpxC") %>% rep(each = 12) %>% rep(2) %>% 
	data.table(strain = .)

followup.3.reps <- 1:3 %>% rep(4*6) %>%
	data.table(rep = .)

followup.3.induced <- c("off","on") %>% rep(each = 3*12) %>%
	data.table(induced = .)

followup.3.drug <- c("colistin", "rifampicin") %>% rep(each = 6) %>% rep(6) %>%
	data.table(drug = .)

followup.3.dose <- c(0,0.80, 0.80, 0) %>% rep(each = 3) %>% rep(6) %>% 
	data.table(dose = .)

followup.3.wells <- cbind(
	wells.filled, 
	followup.3.strains, 
	followup.3.reps, 
	followup.3.induced, 
	followup.3.drug, 
	followup.3.dose)

################################################################################
# combine the data about the identity of each well with the plate reader data
followup.3 <- followup.3 %>% full_join(followup.3.wells)
followup.3[, dose := factor(as.character(dose))]

followup.3 <- followup.3 %>% filter(time <= (18*60*60))

followup.3.plot <- 
	followup.3 %>% 
	filter(!is.na(strain)) %>%
	mutate(hour = time/60/60) %>%
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

print(followup.3.plot)

################################################################################