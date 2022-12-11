library(pacman)

p_load(data.table, tidyverse, hrbrthemes)

followup.5.depletion <- fread(
	"Followups/followup.5.depletion.tsv",
	header = TRUE)

##########################################################################################
# wrangle the data from the infinite machine... 

setnames(followup.5.depletion, t(followup.5.depletion)[,1])
followup.5.depletion <- followup.5.depletion[-1:-2,]
followup.5.depletion <- followup.5.depletion %>% rename(well = `Time [s]`)


followup.5.depletion <- 
	melt(
		followup.5.depletion, 
		id.vars = "well", 
		variable.name = "time", 
		value.name = "OD600", 
		na.rm = TRUE)

followup.5.depletion[, OD600 := as.numeric(OD600)]
followup.5.depletion[, time := as.numeric(levels(time))[time]]

##########################################################################################
# declare a pattern that corresponds to the wells that actually have cells in them

wells.filled <-
	CJ(well.row = toupper(letters[2:7]),
		 well.col = c(1:12))[
		 	, .(well = paste0(well.row, well.col))]

followup.5.depletion.strains <- c("nuoH", "nuoB", "nuoF", "lpxC", "lpxA", "control.2") %>% rep(2) %>% rep(6) %>% 
	data.table(strain = .)

followup.5.depletion.reps <- 4:6 %>% rep(each = 12) %>% rep(2) %>%
	data.table(rep = .)

followup.5.depletion.induced <- c(
	c(0.0, 0.1) %>% rep(each = 6) %>% rep(3), 
	c(0.5, 1.0) %>% rep(each = 6) %>% rep(3)) %>%
	data.table(induced = .)

followup.5.depletion.drug <- c("no drug", "no drug") %>% rep(each = 6) %>% rep(6) %>%
	data.table(drug = .)

followup.5.depletion.dose <- c(0, 0, 0, 0) %>% rep(each = 3) %>% rep(6) %>% 
	data.table(dose = .)

followup.5.depletion.wells <- cbind(
	wells.filled, 
	followup.5.depletion.strains, 
	followup.5.depletion.reps, 
	followup.5.depletion.induced, 
	followup.5.depletion.drug, 
	followup.5.depletion.dose)

##########################################################################################
# combine the data about the identity of each well with the plate reader data
followup.5.depletion <- followup.5.depletion %>% full_join(followup.5.depletion.wells)
followup.5.depletion[, dose := factor(as.character(dose))]

followup.5.depletion <- followup.5.depletion %>% filter(time <= (18*60*60))

followup.5.depletion.plot <- 
	followup.5.depletion %>%
	mutate(hour = time/60/60) %>%
	filter(!is.na(strain)) %>%
	ggplot(aes(x = hour, y = OD600, fill = strain, colour = strain)) + 
	stat_smooth(
		fullrange = TRUE, 
		level = 0.99999,
		method = "gam",
		formula = y ~ s(x, bs = "cs")) +
	facet_wrap(facets = c("drug", "dose", "induced"), ncol = 2) +
	scale_colour_brewer(palette = "Dark2") +
	scale_fill_brewer(palette = "Dark2") +
	ggtitle(bquote(bold("IPTG Growth Phenotypes for")~bolditalic("Acinetobacter baumannii"))) +
	theme_ipsum()

print(followup.5.depletion.plot)

#########################################################################
followup.5.depletion.stat.plot <-
	followup.5.depletion[
		!is.na(strain), 
		SummarizeGrowth(
			time/60/60, 
			OD600,
			blank = followup.5.depletion[
				is.na(strain), 
				median(OD600), 
				by = .(time)]$V1)$vals, 
		by = .(well)] %>% 
	inner_join(followup.5.depletion) %>% 
	mutate(dose = case_when(dose == 0 ~ "no drug", dose != 0 ~ "with drug")) %>%
	ggplot(aes(x = dose, y = t_mid, fill = strain)) + 
	geom_boxplot(position = "dodge") +
	ggtitle(
		bquote(bold("Hours to reach half capacity")~bolditalic("Acinetobacter baumannii."))) +
	theme_ipsum()+
	scale_fill_brewer(palette = "Dark2") +
	ylim(0, NA) +
	facet_wrap(facets = c("drug", "induced"), ncol = 2)

print(followup.5.depletion.stat.plot)

################################################################################
# wells moving forward

CJ(well.row = toupper(letters[1:8]),
	 well.col = c(1:12)) %>% 
	mutate(well = paste0(well.row, well.col)) %>% 
	left_join(
		followup.5.depletion %>% 
			filter(strain %in% c("control.2", "nuoB", "nuoF", "lpxA", "lpxC") & 
						 	induced %in% c(0, 0.5)) %>% 
			select(well, strain, induced, rep) %>% 
			unique %>% 
			mutate(id = paste(strain, rep, induced, sep = ", "))) %>% 
	replace_na(list(id = "---")) %>% 
	pivot_wider(id_cols = well.row, names_from = well.col, values_from = id) %>% 
	print

################################################################################


followup.5.depletion.tSNE <- 
	followup.5.depletion %>%
	filter(!is.na(rep)) %>% 
	select(well, OD600, time) %>% 
	pivot_wider(id_cols = well, names_from = time, values_from = OD600) %>% 
	column_to_rownames("well") %>% 
	Rtsne(perplexity = 10)

followup.5.depletion.tSNE <-
	followup.5.depletion.tSNE$Y %>% 
	data.table

followup.5.depletion %>%
	filter(!is.na(rep)) %>% 
	select(well, OD600, time) %>% 
	pivot_wider(id_cols = well, names_from = time, values_from = OD600) %>%
	select(well) %>%
	inner_join(followup.5.depletion.wells) %>%
	cbind(followup.5.depletion.tSNE) %>%
	mutate(rep = as.character(rep)) %>%
	ggplot(aes(x = V1, y = V2, fill = strain, shape = rep)) +
	geom_point(size = 3, alpha = 0.5) +
	facet_grid(facets = c("induced")) +
	scale_shape_manual(values = c(21, 22, 24)) +
	guides(fill = guide_legend(
		override.aes = list(shape = 23)))

