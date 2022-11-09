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

growth_means_inhibition <- 0.1
fosfomycin.max <- 2500
imipenem.max <- 8


doc_theme <- theme_ipsum(
	base_family = "Arial", 
	caption_margin = 12,
	axis_title_size = 12,
	axis_col = "black")

growth.wt <- fread(
	"Synergy/amy/growth.fos.imi.med.tsv",
	header = T)

fosfomycin.dose <- data.table(
	row = LETTERS[1:8], 
	fosfomycin = c(fosfomycin.max/2^(0:6), 0))

imipenem.dose <- data.table(
	col = 1:8, 
	imipenem = c(imipenem.max/2^(0:6), 0))

growth <- rbind(
	growth.wt %>% 
		pivot_longer(cols = -`<>`, names_to = "col", values_to = "growth") %>% 
		mutate(gene = "wt")) %>%
	rename(row = `<>`) %>%
	group_by(gene) %>% 
	mutate(
		growth.adj = growth - min(growth), 
		growth.prop = growth.adj/max(growth)) %>%
	mutate(col = as.numeric(col)) %>%
	inner_join(imipenem.dose) %>%
	inner_join(fosfomycin.dose) 


growth$gene <- factor(
	growth$gene,
	levels=c("nt.1","nt.2"),
	labels=c("wildtype~1", "wildtype~2"))

mic <- growth %>% 
	filter(growth.prop < growth_means_inhibition) %>% {
		full_join(
			filter(., imipenem == 0) %>% 
				summarise(mic.fosfomycin = min(fosfomycin)),
			filter(., fosfomycin == 0) %>% 
				summarise(mic.imipenem = min(imipenem)))}

growth <- growth %>% 
	inner_join(mic) %>% 
	mutate(
		fic.fosfomycin = fosfomycin/mic.fosfomycin,
		fic.imipenem = imipenem/mic.imipenem)

synergy <- growth %>% 
	filter(fic.fosfomycin + fic.imipenem == 0.5) %>% 
	nest %>% 
	mutate(fit = map(data, ~ glm(fosfomycin ~ imipenem, data = .))) %>% 
	mutate(predictions = map2(
		fit, 
		data, 
		~augment(
			.x, 
			newdata = expand_grid(
				imipenem = seq(
					min((.y %>% filter(fic.imipenem <= 0.5) %>% select(imipenem))), 
					max((.y %>% filter(fic.imipenem <= 0.5) %>% select(imipenem))), 
					length = 32))))) %>% 
	unnest(predictions) %>% 
	rename(fosfomycin = .fitted) %>% 
	select(gene, imipenem, fosfomycin)

min.imipenem <- growth %>% 
	ungroup %>% 
	filter(imipenem != 0) %>% 
	summarise(min(imipenem)) %>% 
	pull

min.fosfomycin <- growth %>% 
	ungroup %>% 
	filter(fosfomycin != 0) %>% 
	summarise(min(fosfomycin)) %>% 
	pull

max.imipenem <- growth %>% 
	ungroup %>% 
	filter(imipenem != 0) %>% 
	summarise(max(imipenem)) %>% 
	pull

max.fosfomycin <- growth %>% 
	ungroup %>% 
	filter(fosfomycin != 0) %>% 
	summarise(max(fosfomycin)) %>% 
	pull


mic.imipenem <- growth %>% 
	filter(fosfomycin == 0 & growth.prop < growth_means_inhibition) %>% 
	summarise(mic.imipenem = min(imipenem))

mic.fosfomycin <- growth %>% 
	filter(imipenem == 0 & growth.prop < growth_means_inhibition) %>% 
	summarise(mic.fosfomycin = min(fosfomycin))

synergy_tally <-
	growth %>% 
	filter(fic.fosfomycin + fic.imipenem < 0.5 & growth.prop < growth_means_inhibition) %>% 
	tally %>%
	mutate(synergy_count = paste(n, "wells"))

growth <- growth %>%
	left_join(synergy_tally)

growth %>% 
	left_join(mic.imipenem) %>% 
	left_join(mic.fosfomycin) %>%
	mutate(
		`Percent Growth` = (100 * (growth.prop)),
		Interaction = case_when(
			imipenem > 0 & 
				fosfomycin > 0 & 
				fic.fosfomycin + fic.imipenem < 0.5 & 
				growth.prop < growth_means_inhibition ~ "Synergy (No Growth)",
			growth.prop < growth_means_inhibition ~ "No Growth",
			growth.prop >= growth_means_inhibition ~ "Growth",
			TRUE ~ "Additive"),
		Interaction = factor(Interaction, levels = c("Growth", "No Growth", "Synergy (No Growth)")),
		imipenem = case_when(
			imipenem == 0 ~ min.imipenem/2, TRUE ~ imipenem),
		fosfomycin = case_when(
			fosfomycin == 0 ~ min.fosfomycin/2, TRUE ~ fosfomycin)) %>%
	ggplot(aes(x = imipenem, y = fosfomycin)) + 
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
			imipenem >= min.imipenem & fosfomycin >= min.fosfomycin) %>% 
			rbind(growth %>% filter(fic.fosfomycin + fic.imipenem == 0.5)),
		lwd = 1, lty = "twodash", 
		colour = "red") + 
	doc_theme +
	scale_x_continuous(trans = 'log10') +
	scale_y_continuous(trans = 'log10') +
geom_richtext(
	data = . %>% filter(fosfomycin == mic.fosfomycin & imipenem == min.imipenem/2),
	aes(label = str_wrap(paste("MIC:",fosfomycin), 12), y = fosfomycin, hjust = 0),
	fill = "cornsilk",
	size = 3) +
	geom_richtext(
		data = . %>% filter(imipenem == mic.imipenem & fosfomycin == min.fosfomycin/2),
		aes(label = str_wrap(paste("MIC:",imipenem), 12), x = imipenem, hjust = 0), 
		angle = 90,
		fill = "cornsilk",
		size = 3) +
	theme(
		legend.position = 'none', legend.direction = "horizontal",
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
		xlim = c(min.imipenem/2.5, max.imipenem*2),
		ylim = c(min.fosfomycin/2.5, max.fosfomycin*2))


##################################################################

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

growth_means_inhibition <- 0.1
fosfomycin.max <- 2500
meropenem.max <- 16


doc_theme <- theme_ipsum(
	base_family = "Arial", 
	caption_margin = 12,
	axis_title_size = 12,
	axis_col = "black")

growth.wt <- fread(
	"Synergy/amy/growth.fos.mer.med.tsv",
	header = T)

fosfomycin.dose <- data.table(
	row = LETTERS[1:8], 
	fosfomycin = c(fosfomycin.max/2^(0:6), 0))

meropenem.dose <- data.table(
	col = 1:8, 
	meropenem = c(meropenem.max/2^(0:6), 0))

growth <- rbind(
	growth.wt %>% 
		pivot_longer(cols = -`<>`, names_to = "col", values_to = "growth") %>% 
		mutate(gene = "wt")) %>%
	rename(row = `<>`) %>%
	group_by(gene) %>% 
	mutate(
		growth.adj = growth - min(growth), 
		growth.prop = growth.adj/max(growth)) %>%
	mutate(col = as.numeric(col)) %>%
	inner_join(meropenem.dose) %>%
	inner_join(fosfomycin.dose) 


growth$gene <- factor(
	growth$gene,
	levels=c("nt.1","nt.2"),
	labels=c("wildtype~1", "wildtype~2"))

mic <- growth %>% 
	filter(growth.prop < growth_means_inhibition) %>% {
		full_join(
			filter(., meropenem == 0) %>% 
				summarise(mic.fosfomycin = min(fosfomycin)),
			filter(., fosfomycin == 0) %>% 
				summarise(mic.meropenem = min(meropenem)))}

growth <- growth %>% 
	inner_join(mic) %>% 
	mutate(
		fic.fosfomycin = fosfomycin/mic.fosfomycin,
		fic.meropenem = meropenem/mic.meropenem)

synergy <- growth %>% 
	filter(fic.fosfomycin + fic.meropenem == 0.5) %>% 
	nest %>% 
	mutate(fit = map(data, ~ glm(fosfomycin ~ meropenem, data = .))) %>% 
	mutate(predictions = map2(
		fit, 
		data, 
		~augment(
			.x, 
			newdata = expand_grid(
				meropenem = seq(
					min((.y %>% filter(fic.meropenem <= 0.5) %>% select(meropenem))), 
					max((.y %>% filter(fic.meropenem <= 0.5) %>% select(meropenem))), 
					length = 32))))) %>% 
	unnest(predictions) %>% 
	rename(fosfomycin = .fitted) %>% 
	select(gene, meropenem, fosfomycin)

min.meropenem <- growth %>% 
	ungroup %>% 
	filter(meropenem != 0) %>% 
	summarise(min(meropenem)) %>% 
	pull

min.fosfomycin <- growth %>% 
	ungroup %>% 
	filter(fosfomycin != 0) %>% 
	summarise(min(fosfomycin)) %>% 
	pull

max.meropenem <- growth %>% 
	ungroup %>% 
	filter(meropenem != 0) %>% 
	summarise(max(meropenem)) %>% 
	pull

max.fosfomycin <- growth %>% 
	ungroup %>% 
	filter(fosfomycin != 0) %>% 
	summarise(max(fosfomycin)) %>% 
	pull


mic.meropenem <- growth %>% 
	filter(fosfomycin == 0 & growth.prop < growth_means_inhibition) %>% 
	summarise(mic.meropenem = min(meropenem))

mic.fosfomycin <- growth %>% 
	filter(meropenem == 0 & growth.prop < growth_means_inhibition) %>% 
	summarise(mic.fosfomycin = min(fosfomycin))

synergy_tally <-
	growth %>% 
	filter(fic.fosfomycin + fic.meropenem < 0.5 & growth.prop < growth_means_inhibition) %>% 
	tally %>%
	mutate(synergy_count = paste(n, "wells"))

growth <- growth %>%
	left_join(synergy_tally)

growth %>% 
	left_join(mic.meropenem) %>% 
	left_join(mic.fosfomycin) %>%
	mutate(
		`Percent Growth` = (100 * (growth.prop)),
		Interaction = case_when(
			meropenem > 0 & 
				fosfomycin > 0 & 
				fic.fosfomycin + fic.meropenem < 0.5 & 
				growth.prop < growth_means_inhibition ~ "Synergy (No Growth)",
			growth.prop < growth_means_inhibition ~ "No Growth",
			growth.prop >= growth_means_inhibition ~ "Growth",
			TRUE ~ "Additive"),
		Interaction = factor(Interaction, levels = c("Growth", "No Growth", "Synergy (No Growth)")),
		meropenem = case_when(
			meropenem == 0 ~ min.meropenem/2, TRUE ~ meropenem),
		fosfomycin = case_when(
			fosfomycin == 0 ~ min.fosfomycin/2, TRUE ~ fosfomycin)) %>%
	ggplot(aes(x = meropenem, y = fosfomycin)) + 
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
			meropenem >= min.meropenem & fosfomycin >= min.fosfomycin) %>% 
			rbind(growth %>% filter(fic.fosfomycin + fic.meropenem == 0.5)),
		lwd = 1, lty = "twodash", 
		colour = "red") + 
	doc_theme +
	scale_x_continuous(trans = 'log10') +
	scale_y_continuous(trans = 'log10') +
	geom_richtext(
		data = . %>% filter(fosfomycin == mic.fosfomycin & meropenem == min.meropenem/2),
		aes(label = str_wrap(paste("MIC:",fosfomycin), 12), y = fosfomycin, hjust = 0),
		fill = "cornsilk",
		size = 3) +
	geom_richtext(
		data = . %>% filter(meropenem == mic.meropenem & fosfomycin == min.fosfomycin/2),
		aes(label = str_wrap(paste("MIC:",meropenem), 12), x = meropenem, hjust = 0), 
		angle = 90,
		fill = "cornsilk",
		size = 3) +
	theme(
		legend.position = 'none', legend.direction = "horizontal",
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
		xlim = c(min.meropenem/2.5, max.meropenem*2),
		ylim = c(min.fosfomycin/2.5, max.fosfomycin*2))


