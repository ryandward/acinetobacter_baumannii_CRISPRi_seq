library(pacman)

p_load(data.table, tidyverse)

followup.1 <- fread(
	"Followups/20220520 abau followup.tsv",
	header = TRUE)

##########################################################################################
# wrangle the data from the sunrise machine... 

followup.1 <- followup.1[`Raw data` %like% "[0-9]+s"]
followup.1 <- followup.1[, -2]
followup.1[, `Raw data` := gsub("s", "", `Raw data`)]
followup.1 <- t(followup.1)
colnames(followup.1) <- followup.1[1, ]
followup.1 <- followup.1[-1, ]

transposed_wells <- 
	CJ(well.row = toupper(letters[1:8]),
		 well.col = c(1:12))[
		 	, .(well = paste0(well.row, well.col))]

rownames(followup.1) <- transposed_wells$well
followup.1 <- data.table(followup.1, keep.rownames = "well")

followup.1 <- 
	melt(
		followup.1, 
		id.vars = "well", 
		variable.name = "time", 
		value.name = "OD600", 
		na.rm = TRUE)

followup.1[, OD600 := as.numeric(OD600)]
followup.1[, time := as.numeric(levels(time))[time]]

##########################################################################################
# declare a pattern that corresponds to the wells that actually have cells in them

wells.filled <-
	CJ(well.row = toupper(letters[2:7]),
		 well.col = c(2:10))[
		 	, .(well = paste0(well.row, well.col))]

followup.1.strains <- c("control.2", "nuoH", "nuoB") %>% 
	rep(18) %>% 
	data.table(strain = .)

followup.1.reps <- rep(c(1:3), each = 9) %>% 
	rep(2) %>%
	data.table(rep = .)

followup.1.induced <- rep(c("off","on"), each = 27) %>%
	data.table(induced = .)

followup.1.drug <- rep("colistin", 54) %>% 
	data.table(drug = .)

followup.1.dose <- rep(c(0, 0.44, 0.22), each = 3) %>% 
	rep(6) %>% 
	data.table(dose = .)

followup.1.wells <- cbind(
	wells.filled, 
	followup.1.strains, 
	followup.1.reps, 
	followup.1.induced, 
	followup.1.drug, 
	followup.1.dose)

##########################################################################################
# combine the data about the identity of each well with the plate reader data

followup.1 <- followup.1 %>% full_join(followup.1.wells)
followup.1[, dose := factor(as.character(dose))]

followup.1.plot <- 
	followup.1 %>% 
	filter(!is.na(strain)) %>%
 	ggplot(aes(x = time, y = OD600, fill = dose, colour = dose)) + 
	stat_smooth(fullrange = TRUE, level = 0.999) +
	facet_wrap(facets = c("drug", "strain", "induced"), ncol = 2) +
	scale_colour_brewer(palette = "Set1") +
	scale_fill_brewer(palette = "Set1") +
	theme_ipsum()


print(followup.1.plot)


followup.1.plot <- 
	followup.1 %>% 
	filter(!is.na(strain)) %>%
	ggplot(aes(x = time, y = OD600, fill = dose, colour = strain)) + 
	stat_smooth(fullrange = TRUE, level = 0.999) +
	facet_wrap(facets = c("drug", "dose", "induced"), ncol = 2) +
	scale_colour_brewer(palette = "Set1") +
	scale_fill_brewer(palette = "Set1") +
	theme_ipsum()


print(followup.1.plot)

################################################################################
