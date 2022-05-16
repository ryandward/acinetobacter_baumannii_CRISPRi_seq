mixed <- fread("../acinetobacter_baumannii_chemical_genomics/fosf.tsv") %>% 
	as_tibble %>% 
	dplyr::select(AB19606, medLFC, FDR, condition, type) %>% 
	rbind(median_melted_results[condition %like% "Imipenem_0.09_T2 - None_0_T2"] %>% 
					dplyr::select(AB19606, medLFC, FDR, condition, type)) %>%	
	filter(type == "perfect") %>%
	mutate(
		depleted = (medLFC < 1 & FDR < 0.01), 
		resistant = (medLFC > 1 & FDR < 0.01)) %>% 
	select(AB19606, condition, depleted, resistant) %>% 
	pivot_longer(!c(condition, AB19606)) %>% 
	mutate(condition = gsub("_[0-9\\.]+_", " ", condition)) %>% 
	unite("condition", c("condition", "name"), sep = " ") %>% 
	pivot_wider(id_cols = AB19606, names_from = condition, values_from = value) 

mixed.matrix <- mixed %>%
	make_comb_mat(mode = "distinct")

p <- UpSet(
	mixed.matrix, 
	top_annotation = HeatmapAnnotation(
		"Intersection" = anno_barplot(
			comb_size(mixed.matrix),
			border = FALSE, 
			height = unit(8, "cm"),
			add_numbers = T,
			gp = gpar(fill = c("grey", rep("red", 1), rep("grey", 6)), lty = "blank")), 
		show_annotation_name = FALSE),
	right_annotation = rowAnnotation(
		"Genes Significant Beyond Induction at T2" = anno_barplot(
			set_size(mixed.matrix),
			border = FALSE,
			gp = gpar(
				# fill = viridis(set_size(mixed.matrix)/2 %>% length, direction = -1, alpha = 0.5),
				# fill = c(rep("light blue", 2), rep("coral", 2)),
				fill = viridis(2, direction = -1, alpha = 0.35) %>% rep(each = 2), 
				col = c("red", "black", "red", "black"),
				lwd = 3,
				lty = c("dotted", "solid", "dotted", "solid")),
			width = unit(6, "cm"),
			add_numbers = T),
		annotation_name_gp = gpar(fontsize = 8)),
	set_order = set_name(mixed.matrix),
	row_names_gp = grid::gpar(fontsize = 10),
	comb_col = c("black", rep("red", 1), rep("black", 6)),
	row_labels = c(
		as.expression(bquote(Fosfomycin[" Vulnerable "])), 
		as.expression(bquote(Fosfomycin[" Resistant "])),
		as.expression(bquote(Imipenem[" Vulnerable "])),
		as.expression(bquote(Imipenem[" Resistant "]))))

print(p)