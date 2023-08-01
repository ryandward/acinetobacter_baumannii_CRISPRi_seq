conflict_prefer("pheatmap", "pheatmap")

generate_breaks = function(x, n, center = F) {
	
	if (center) {
		m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
		res = seq(-m, m, length.out = n + 1)}
	
	else {
		res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)}
	
	return(res)
}

plot_colors <- c(
	colorRampPalette(
		c("#9a0007", "white"))(144^3) %>% 
		`[`((1:144)^3) %>% 
		`[`(-144),
	"white",
	rev(colorRampPalette(
		c("#005b9f", "white"))(144^3) %>% 
			`[`((1:144)^3) %>% 
			`[`(-144)))

###################################################################
###################################################################

conditions.formatted <- c(
	bquote("Meropenem (Low)"),
	bquote("Meropenem (High)"),
	bquote("Imipenem (Low)"),
	bquote("Imipenem (High)"))

trna <- median_melted_results %>%
	filter(condition %in% c(
		"Meropenem_0.17_T1 - None_0_T1", 
		"Meropenem_0.11_T1 - None_0_T1",
		"Imipenem_0.06_T1 - None_0_T1",
		"Imipenem_0.09_T1 - None_0_T1")) %>% 
	filter(Pathway %like% "tRNA") %>% 
	filter(type == "perfect") %>% 
	pivot_wider(id_cols = c(Pathway, unique_name), names_from = condition, values_from = medLFC) %>%
	data.table

trna.mat <- trna %>% select(-Pathway, -unique_name) %>% data.matrix 

rownames(trna.mat) <- trna %>% select(unique_name) %>% pull

breaks <- generate_breaks(trna.mat, n = 286, center = T)

trna.formatted <- lapply(
	trna[, .I],
	function(x) { 
		bquote(italic(.(trna[x]$unique_name)))})

trna.mat %>% pheatmap(
	cutree_row = 3, 
	breaks = breaks, 
	col = plot_colors,
	clustering_distance_cols = "canberra",
	clustering_distance_rows = "canberra",
	clustering_method = "ward.D2",
	cluster_cols = FALSE,
	angle_col = 315,
	border_color = NA,
	labels_col = conditions.formatted,
	labels_row = as.expression(trna.formatted))

trna.groups <- 
	trna.mat %>% dist(method = "canberra") %>% 
	hclust(method = "ward.D2") %>% 
	cutree(k = 3) %>% data.matrix %>% 
	`colnames<-`(c("group")) %>% 
	data.table(keep.rownames = "unique_name")


####

pg_div <- median_melted_results %>%
	filter(condition %in% c(
		"Meropenem_0.17_T1 - None_0_T1", 
		"Meropenem_0.11_T1 - None_0_T1",
		"Imipenem_0.06_T1 - None_0_T1",
		"Imipenem_0.09_T1 - None_0_T1")) %>% 
	filter(Pathway %like% "Cell Wall") %>% 
	filter(type == "perfect") %>% 
	pivot_wider(id_cols = c(Pathway, unique_name), names_from = condition, values_from = medLFC) %>%
	data.table

pg_div.mat <- pg_div %>% select(-Pathway, -unique_name) %>% data.matrix 

rownames(pg_div.mat) <- pg_div %>% select(unique_name) %>% pull

breaks <- generate_breaks(pg_div.mat, n = 286, center = T)

pg_div.formatted <- lapply(
	pg_div[, .I],
	function(x) { 
		bquote(italic(.(pg_div[x]$unique_name)))})

pg_div.mat %>% pheatmap(
	cutree_row = 3, 
	breaks = breaks, 
	col = plot_colors,
	clustering_distance_cols = "canberra",
	clustering_distance_rows = "canberra",
	clustering_method = "ward.D2",
	cluster_cols = FALSE,
	angle_col = 315,
	border_color = NA,
	labels_col = conditions.formatted,
	labels_row = as.expression(pg_div.formatted))

pg_div.groups <- 
	pg_div.mat %>% dist(method = "canberra") %>% 
	hclust(method = "ward.D2") %>% 
	cutree(k=3) %>% data.matrix %>% 
	`colnames<-`(c("group")) %>% 
	data.table(keep.rownames = "unique_name")



###################################################################
###################################################################

melted_results %>% 
	filter(type == "mismatch") %>%
	filter(condition %in% c("Meropenem_0.17_T1 - None_0_T1", "None_0_T1 - None_0_T0")) %>% 
	mutate(Condition = case_when(
		condition %like% "Meropenem" ~ "Meropenem",
		condition %like% "None_0_T0" ~ "Inducer"
	)) %>%
	pivot_wider(id_cols = c(Pathway, unique_name, spacer), names_from = Condition, values_from = LFC.adj) %>% 
	ggplot(aes(x = Inducer, y = Meropenem)) + 
	geom_point(alpha = 0.15) +
	geom_smooth(colour = "black", fill = "black",
							# data = . %>% filter(Pathway %like% "PG" | Pathway %like% "tRNA" | is.na(Pathway)),
							aes(), 
							method = glm, 
							formula = y ~ x + 0) + 
	ylim(-5.75, 3) +
	doc_theme

#################################


meropenem_subset <- melted_results %>% 
	filter(type == "mismatch") %>%
	mutate(Pathway = case_when(
		Pathway == "Ribosome" ~ "Ribosome",
		# Pathway == "LOS" ~ "LOS",
		# unique_name %like% "nuo" ~ "NDH-1",
		Pathway %like% "Cell Wall" & unique_name %in% (pg_div.groups %>% filter(group == 2) %>% pull(unique_name)) ~ "PG/Division (Responsive)",
		Pathway %like% "Cell Wall" & unique_name %in% (pg_div.groups %>% filter(group != 2) %>% pull(unique_name)) ~ "PG/Division (Other)",
		Pathway %like% "tRNA" & unique_name %in% (trna.groups %>% filter(group == 1) %>% pull(unique_name)) ~  "tRNA Ligase (Responsive)",
		Pathway %like% "tRNA" & unique_name %in% (trna.groups %>% filter(group != 1) %>% pull(unique_name)) ~ "tRNA Ligase (Other)",
		TRUE ~ "Other")) %>%
	filter(condition %in% c("Meropenem_0.17_T1 - None_0_T1", "None_0_T1 - None_0_T0")) %>% 
	mutate(Condition = case_when(
		condition %like% "Meropenem" ~ "Meropenem",
		condition %like% "None_0_T0" ~ "Inducer"
	)) %>%
	pivot_wider(id_cols = c(Pathway, unique_name, spacer), names_from = Condition, values_from = LFC.adj) %>%
	group_by(Pathway)

meropenem_subset %>%
	ggplot(aes(x = Inducer, y = Meropenem)) + 
	geom_point(aes(colour = Pathway), alpha = 0.15) +
	geom_smooth(
		# data = . %>% filter(Pathway %like% "PG" | Pathway %like% "tRNA" | is.na(Pathway)),
		aes(fill = Pathway, colour = Pathway), 
		method = glm, 
		formula = y ~ x + 0)+
	scale_colour_manual(
		na.value = alpha("dark grey", 1),
		values = c(
			"Ribosome" = alpha("#E31A1C", 1),
			# "NDH-1" = alpha("#6A3D9A", 1),
			# "LOS" = alpha("#33a02c", 1),
			"PG/Division (Responsive)" = alpha("#ff7f00", 1),
			"PG/Division (Other)" = alpha("#FDBF6F", 1),
			"tRNA Ligase (Responsive)" = alpha("#1F78B4", 1),
			"tRNA Ligase (Other)" = alpha("#A6CEE3", 1),
			"Other" = "#808080")) +
	scale_fill_manual(
		na.value = alpha("dark grey", 0.5),
		values = c(
			"Ribosome" = "#E31A1C",
			# "NDH-1" = "#6A3D9A",
			# "LOS" = "#33a02c",
			"PG/Division (Responsive)" = alpha("#ff7f00", 0.5),
			"PG/Division (Other)" = alpha("#FDBF6F", 0.5),
			"tRNA Ligase (Responsive)" = alpha("#1F78B4", 0.5),
			"tRNA Ligase (Other)" = alpha("#A6CEE3", 0.5),
			"Other" = "#808080")) +
	ylim(-5.75, 3) +
	theme(legend.position = "none") + 		
	guides(colour = guide_legend(nrow = 3)) + 
	doc_theme +
	theme(legend.position = "none")

###################################################################
###################################################################

meropenem_subset.cors <-meropenem_subset %>%
	summarise(
		p.cor = cor.test(Meropenem, Inducer)$p.value,
		# cor = cor.test(Meropenem, Inducer)$statistic, 
		r = cor.test(Meropenem, Inducer)$estimate,
		r_sq = cor.test(Meropenem, Inducer)$estimate^2,
		# method = cor.test(Meropenem, Inducer)$method,
		df.guides = cor.test(Meropenem, Inducer)$parameter)

meropenem_subset.regressions <- meropenem_subset %>% 
	nest(data = -Pathway) %>%
	mutate(
		fit = map(data, ~ lm(Meropenem ~ Inducer + 0, data = .x)),
		tidied = map(fit, tidy)
	) %>% 
	unnest(tidied)

meropenem.stats <- meropenem_subset.cors %>% 
	inner_join(meropenem_subset.regressions) %>% 
	select(Pathway, r_sq, p.cor, estimate, std.error, p.value, df.guides) %>% 
	rename (slope = estimate, std.err.slope = std.error, p.slope = p.value)

###################################################################

# parameter fits for all genes
meropenem_subset.cors <- meropenem_subset %>%
	mutate(Pathway = "All") %>%
	summarise(
		p.cor = cor.test(Meropenem, Inducer)$p.value,
		# cor = cor.test(Meropenem, Inducer)$statistic, 
		r = cor.test(Meropenem, Inducer)$estimate,
		r_sq = cor.test(Meropenem, Inducer)$estimate^2,
		# method = cor.test(Meropenem, Inducer)$method,
		df.guides = cor.test(Meropenem, Inducer)$parameter) 

meropenem_subset.regressions <- meropenem_subset %>%
	mutate(Pathway = "All") %>%
	nest(data = -Pathway) %>%
	mutate(
		fit = map(data, ~ lm(Meropenem ~ Inducer + 0, data = .x)),
		tidied = map(fit, tidy)
	) %>% 
	unnest(tidied)

meropenem.stats <- meropenem.stats %>%
	rbind(meropenem_subset.cors %>% 
	inner_join(meropenem_subset.regressions) %>% 
	select(Pathway, r_sq, p.cor, estimate, std.error, p.value, df.guides) %>% 
	rename (slope = estimate, std.err.slope = std.error, p.slope = p.value))

meropenem.stats %>% print