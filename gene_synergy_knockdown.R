library(tidyverse)
library(data.table)
library(scales)
library(ggallin)
library(hrbrthemes)
library(broom)
library(ggtext)
library(conflicted)

conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")

growth_means_inhibition <- 0.15
colistin.max <- 3
rifampicin.max <- 1.28


doc_theme <- theme_ipsum(
	base_family = "Arial", 
	caption_margin = 12,
	axis_title_size = 12,
	axis_col = "black")


growth.nt <- fread(
	"/Users/ryandward/R/acinetobacter_baumannii_CRISPRi_seq/synergy/growth.nt.tsv",
	header = T)
growth.lpxC <- fread(
	"/Users/ryandward/R/acinetobacter_baumannii_CRISPRi_seq/synergy/growth.lpxC.tsv",
	header = T)
growth.nuoB <- fread(
	"/Users/ryandward/R/acinetobacter_baumannii_CRISPRi_seq/synergy/growth.nuoB.tsv",
	header = T)

colistin.dose <- data.table(
	row = LETTERS[1:8], 
	colistin = c(colistin.max/2^(0:6), 0))

rifampicin.dose <- data.table(
	col = 3:12, 
	rifampicin = c(rifampicin.max/2^(0:8), 0))


growth <- rbind(
	growth.nt %>% 
		pivot_longer(cols = -`<>`, names_to = "col", values_to = "growth") %>% 
		mutate(gene = "nt"),
	growth.lpxC %>% 
		pivot_longer(cols = -`<>`, names_to = "col", values_to = "growth") %>% 
		mutate(gene = "lpxC"),
	growth.nuoB %>% 
		pivot_longer(cols = -`<>`, names_to = "col", values_to = "growth") %>% 
		mutate(gene = "nuoB")) %>%
	rename(row = `<>`) %>%
	group_by(gene) %>% 
	mutate(
		growth.adj = growth - min(growth), 
		growth.prop = growth.adj/max(growth)) %>%
	mutate(col = as.numeric(col)) %>%
	inner_join(rifampicin.dose) %>%
	inner_join(colistin.dose) %>%
	filter(colistin < 3 & rifampicin < 1.28) 


growth$gene <- factor(
	growth$gene,
	levels=c("nt","lpxC","nuoB"),
	labels=c("non-targeting", "bolditalic(lpxC)","bolditalic(nuoB)"))

# growth %>% filter(
#   growth.prop < 0.1 & rifampicin == 0) %>% 
#   summarise(mic.colistin = min(colistin))

# growth %>% arrange(desc(colistin), desc(rifampicin)) %>% filter(gene == "nuoB") %>% pivot_wider(id_cols = colistin, names_from = rifampicin, values_from = growth.prop)

mic <- growth %>% 
	filter(growth.prop < growth_means_inhibition) %>% {
		full_join(
			filter(., rifampicin == 0) %>% 
				summarise(mic.colistin = min(colistin)),
			filter(., colistin == 0) %>% 
				summarise(mic.rifampicin = min(rifampicin)))}

growth <- growth %>% 
	inner_join(mic) %>% 
	mutate(
		fic.colistin = colistin/mic.colistin,
		fic.rifampicin = rifampicin/mic.rifampicin)

synergy <- growth %>% 
	filter(fic.colistin + fic.rifampicin == 0.5) %>% 
	nest %>% 
	mutate(fit = map(data, ~ glm(colistin ~ rifampicin, data = .))) %>% 
	mutate(predictions = map2(
		fit, 
		data, 
		~augment(
			.x, 
			newdata = expand_grid(
				rifampicin = seq(
					min((.y %>% filter(fic.rifampicin <= 0.5) %>% select(rifampicin))), 
					max((.y %>% filter(fic.rifampicin <= 0.5) %>% select(rifampicin))), 
					length = 32))))) %>% 
	unnest(predictions) %>% 
	rename(colistin = .fitted) %>% 
	select(gene, rifampicin, colistin)

min.rifampicin <- growth %>% 
	ungroup %>% 
	filter(rifampicin != 0) %>% 
	summarise(min(rifampicin)) %>% 
	pull

min.colistin <- growth %>% 
	ungroup %>% 
	filter(colistin != 0) %>% 
	summarise(min(colistin)) %>% 
	pull

max.rifampicin <- growth %>% 
	ungroup %>% 
	filter(rifampicin != 0) %>% 
	summarise(max(rifampicin)) %>% 
	pull

max.colistin <- growth %>% 
	ungroup %>% 
	filter(colistin != 0) %>% 
	summarise(max(colistin)) %>% 
	pull


mic.rifampicin <- growth %>% 
	filter(colistin == 0 & growth.prop < growth_means_inhibition) %>% 
	summarise(mic.rifampicin = min(rifampicin))

mic.colistin <- growth %>% 
	filter(rifampicin == 0 & growth.prop < growth_means_inhibition) %>% 
	summarise(mic.colistin = min(colistin))

synergy_tally <-
	growth %>% 
	filter(fic.colistin + fic.rifampicin < 0.5 & growth.prop < growth_means_inhibition) %>% 
	tally %>%
	mutate(synergy_count = paste(n, "wells"))

growth <- growth %>%
	inner_join(synergy_tally)



growth %>% 
	left_join(mic.rifampicin) %>% 
	left_join(mic.colistin) %>%
	mutate(
		`Percent Growth` = (100 * (growth.prop)),
		Interaction = case_when(
			rifampicin > 0 & 
				colistin > 0 & 
				fic.colistin + fic.rifampicin < 0.5 & 
				growth.prop < growth_means_inhibition ~ "Synergy (No Growth)",
			growth.prop < growth_means_inhibition ~ "No Growth",
			growth.prop >= growth_means_inhibition ~ "Growth",
			TRUE ~ "Additive"),
		Interaction = factor(Interaction, levels = c("Growth", "No Growth", "Synergy (No Growth)")),
		rifampicin = case_when(
			rifampicin == 0 ~ min.rifampicin/2, TRUE ~ rifampicin),
		colistin = case_when(
			colistin == 0 ~ min.colistin/2, TRUE ~ colistin)) %>%
	ggplot(aes(x = rifampicin, y = colistin)) + 
	# ggtitle("Gene Knockdown Alters Synergy") +
	geom_point(
		pch = 21, 
		cex = 5,
		alpha = 0.75,
		aes(fill = Interaction)) +
	# geom_point(
	#   # alpha = 1,
	#   pch = 20,
	#   cex = 5,
	#   colour = "black",
	#   aes(alpha = Growth)) +
	geom_line(
		data = synergy %>% filter(
			rifampicin >= min.rifampicin & colistin >= min.colistin) %>% 
			rbind(growth %>% filter(fic.colistin + fic.rifampicin == 0.5)),
		lwd = 1, lty = "twodash", 
		colour = "red") + 
	doc_theme +
	scale_x_continuous(trans = 'log10') +
	scale_y_continuous(trans = 'log10') +
	
	# geom_richtext(
	#   data = . %>%
	#     filter( 
	#       (colistin == mic.colistin ) &
	#         (rifampicin == mic.rifampicin )),
	#   aes(
	#     x = rifampicin,
	#     y = min.colistin,
	#     label = paste("MIC =", rifampicin)),
	#   angle = 90) +
# geom_richtext(
#   data = . %>%
#     filter( 
#       (rifampicin == mic.rifampicin ) &
#         (colistin == mic.colistin )),
#   aes(
#     x = min.rifampicin,
#     y = colistin,
#     label = paste("MIC =", colistin))) +
# geom_line(
#   colour = "dark cyan",
#   lwd = 2,
#   alpha = 0.25,
#   data = . %>% 
#     group_by(gene) %>% {
#       full_join(
#         filter(., growth.prop < growth_means_inhibition) %>% 
#           group_by(colistin, fic.colistin, gene) %>% 
#           summarise(rifampicin = min(rifampicin), fic.rifampicin = min(fic.rifampicin)),
#         filter(., growth.prop < growth_means_inhibition) %>% 
#           group_by(rifampicin, fic.rifampicin, gene) %>% 
#           summarise(colistin = min(colistin), fic.colistin = min(fic.colistin)))
#       } %>%
#     arrange(
#       rifampicin, desc(colistin))) +
geom_richtext(
	data = . %>% filter(colistin == mic.colistin & rifampicin == min.rifampicin/2),
	aes(label = str_wrap(paste("MIC:",colistin), 12), y = colistin, hjust = 0),
	fill = "cornsilk",
	size = 3) +
	geom_richtext(
		data = . %>% filter(rifampicin == mic.rifampicin & colistin == min.colistin/2),
		aes(label = str_wrap(paste("MIC:",rifampicin), 12), x = rifampicin, hjust = 0), 
		angle = 90,
		fill = "cornsilk",
		size = 3) +
	theme(
		legend.position = 'bottom', legend.direction = "horizontal",
		plot.background = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
	scale_fill_manual(
		values = c("black", "white", "red"),
		na.value = "white") +
	scale_alpha_manual(
		values = c(100, 0)) +
	facet_wrap(
		facets = c("gene"),
		labeller = label_parsed) +
	coord_cartesian(
		xlim = c(min.rifampicin/2.5, max.rifampicin*2),
		ylim = c(min.colistin/2.5, max.colistin*2))





