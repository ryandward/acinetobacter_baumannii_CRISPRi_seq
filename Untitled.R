

aba %>% inner_join(aba.design) %>% filter(dose != 0.06 & dose != 0.11) %>% filter(!verbose %like% "No drug") %>% mutate(condition = paste(drug, timing, rep)) %>% pivot_wider(id_cols = spacer, names_from = condition, values_from = count, values_fill = 0) -> z

z %>% select(-spacer) %>% data.matrix %>% cor -> plot_matrix

generate_breaks = function(x, n, center = F) {
	
	if (center) {
		m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
		res = seq(-m, m, length.out = n + 1)}
	
	else {
		res = seq(min(0, na.rm = T), max(x, na.rm = T), length.out = n + 1)}
	
	return(res)
}

breaks <- generate_breaks(plot_matrix, n = 286, center = F)

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

pheatmap(
	plot_matrix,
	col = plot_colors,
	breaks = breaks,
	border_color = NA,
	angle_col = 45,
	cutree_rows = 3,
	cutree_cols = 3,
	clustering_distance_rows = "maximum",
	clustering_distance_cols = "maximum")

################################################################################


Cell_Wall <- fread("CL707.tsv")

median_melted_results %>%
	mutate(Pathway = case_when(
		AB030 %in% Cell_Wall$AB030 ~ "Cell Wall",
		TRUE ~ Pathway)) %>%
	filter(type == "perfect") %>% 
	filter(condition %in% c(
		"Rifampicin_0.34_T2 - None_0_T2", 
		"Colistin_0.44_T2 - None_0_T2",
		"Imipenem_0.09_T2 - None_0_T2",
		"Meropenem_0.17_T2 - None_0_T2")) %>% 
	mutate(condition = gsub("_.*", "", condition)) %>%
	pivot_wider(id_cols = c(Pathway, unique_name), 
							names_from = condition, values_from = medLFC) -> z

z %>% filter(Pathway %in% c("LOS", "NADH")) %>% select(-Pathway, -unique_name) -> zz
														
z_title <- z %>% filter(Pathway %in% c("LOS", "NADH")) %>% select(Pathway, unique_name) %>% data.table


genes.formatted <- lapply(
	z_title[, .I],
	function(x) { 
		bquote(bold(.(z_title[x]$Pathway)) ~ italic(.(z_title[x]$unique_name)))})

zz %>% data.matrix -> zz				

rownames(zz) <- z %>% filter(Pathway %in% c("LOS", "NADH", "Cell Wall")) %>% select(Pathway, unique_name) %>% 
	mutate(gene = paste(Pathway, unique_name)) %>% pull(gene)

breaks <- generate_breaks(zz, n = 286, center = T)

pheatmap(
	zz,
	col = plot_colors,
	breaks = breaks,
	border_color = NA,
	angle_col = 45,
	cutree_rows = ,
	cutree_cols = 4,
	labels_row = as.expression(genes.formatted),
	clustering_method = "ward.D2",
	clustering_distance_rows = "correlation",
	clustering_distance_cols = "correlation")


y %>% group_by(condition, spacer) %>% summarise(count = median(count)) %>% inner_join(aba.key) %>%  pivot_wider(id_cols = spacer, names_from = condition, values_from = count, values_fill = 0) %>% 
	inner_join(yy %>% mutate(Pathway = case_when(Pathway %in% c("LOS", "NADH") ~ Pathway, TRUE ~ NA_character_))) %>%
	ggplot(aes(x = `Colistin T2`, y = `Rifampicin T2`)) + 
	geom_point(data = . %>% filter(is.na(Pathway)), alpha = 0.50, colour = "grey") + 
	geom_point(data = . %>% filter(Pathway %in% c("LOS", "NADH")), aes(colour = Pathway)) +
	scale_x_continuous(
	trans = 'log10',
	breaks = c(0, 10^seq(0,7)),
	labels = label_number(scale_cut = cut_short_scale())) + scale_y_continuous(
		trans = 'log10',
		breaks = c(0, 10^seq(0,7)),
		labels = label_number(scale_cut = cut_short_scale())) + doc_theme


median_melted_results %>%
	mutate(Pathway = case_when(
		AB030 %in% Cell_Wall$AB030 ~ "Cell Wall",
		TRUE ~ Pathway)) %>%
	filter(type == "perfect") %>% 
	filter(!condition %like% "0.06" & !condition %like% "0.11") %>%
	filter(condition %in% interest$condition) %>% filter(!condition %like% "None_0_T0") %>%
	mutate(condition = gsub(" - None_0_T[0-9]", "", condition)) %>%
	mutate(condition = gsub("_", " ", condition)) %>%
	mutate(condition = gsub("T", "T", condition)) %>%
	filter(condition %like% "T2") %>%
	mutate(condition = gsub("T2", "", condition)) %>%
	mutate(condition = gsub(" 0\\.[0-9][0-9]", "", condition)) %>%
 	pivot_wider(id_cols = unique_name, names_from = condition, values_from = medLFC, values_fill = 0) %>% select(-unique_name) %>% cor  -> qq

breaks <- generate_breaks(qq, 286, center = T)

pheatmap(
	qq, 
	breaks = breaks, 
	color = plot_colors, 
	clustering_method = "ward.D2",
	clustering_distance_rows = "canberra",
	clustering_distance_cols = "canberra",
	cutree_rows = 3,
	cutree_cols = 3,
	border_color = NA,
	angle_col = 45,
	display_numbers =T, 
	fontsize_number = 16,
	number_color = "black")
	