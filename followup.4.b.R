library(pacman)

p_load(data.table, tidyverse, hrbrthemes)

p_load_current_gh("briandconnelly/growthcurve")

followup.4.b <- fread(
	"Followups/followup.4.b.tsv",
	header = TRUE)

##########################################################################################
# wrangle the data from the sunrise machine... 

followup.4.b <- followup.4.b[, -c(1,3)] %>% 
	melt(id.vars = "Time [s]", variable.name = "well", value.name = "OD600") %>% 
	rename(time = "Time [s]")

followup.4.b[, OD600 := as.numeric(OD600)]
followup.4.b[, time := as.numeric(time)]

##########################################################################################
# declare a pattern that corresponds to the wells that actually have cells in them

wells.filled <-
	CJ(well.row = toupper(letters[2:7]),
		 well.col = c(1:12))[
		 	, .(well = paste0(well.row, well.col))]

followup.4.b.strains <- c("control.2", "lpxC", "nuoH") %>% rep(each = 12) %>% rep(2) %>% 
	data.table(strain = .)

followup.4.b.reps <- 1:3 %>% rep(4*6) %>%
	data.table(rep = .)

followup.4.b.induced <- c("off","on") %>% rep(each = 3*12) %>%
	data.table(induced = .)

followup.4.b.drug <- c("colistin", "rifampicin") %>% rep(each = 6) %>% rep(6) %>%
	data.table(drug = .)

followup.4.b.dose <- c(0,0.80, 0.40, 0) %>% rep(each = 3) %>% rep(6) %>% 
	data.table(dose = .)

followup.4.b.wells <- cbind(
	wells.filled, 
	followup.4.b.strains, 
	followup.4.b.reps, 
	followup.4.b.induced, 
	followup.4.b.drug, 
	followup.4.b.dose)

##########################################################################################
# combine the data about the identity of each well with the plate reader data
followup.4.b <- followup.4.b %>% full_join(followup.4.b.wells)
followup.4.b[, dose := factor(as.character(dose))]

# followup.4.b <- followup.4.b %>% filter(time <= (18*60*60))

followup.4.b.plot <- 
	followup.4.b %>% 
	mutate(hour = time/60/60) %>%
	filter(!is.na(strain)) %>%
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

print(followup.4.b.plot)

################################################################################
followup.4.b.stat.plot <-
	followup.4.b[
	!is.na(strain), 
	SummarizeGrowth(
		time/60/60, 
		OD600,
		blank = followup.4.b[
			is.na(strain), 
			median(OD600), 
			by = .(time)]$V1)$vals, 
	by = .(well)] %>% 
	inner_join(followup.4.b) %>% 
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

print(followup.4.b.stat.plot)
