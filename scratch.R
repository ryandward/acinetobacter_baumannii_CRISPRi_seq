p_load(tidyverse, ComplexHeatmap, colorspace, viridis)

conditions.full <- tibble(condition = c(
	"None_0_T1 - None_0_T0",
	"None_0_T2 - None_0_T0",
	"Colistin_0.44_T1 - None_0_T0",
	"Colistin_0.44_T2 - None_0_T0",
	"Rifampicin_0.34_T1 - None_0_T0",
	"Rifampicin_0.34_T2 - None_0_T0",
	"Imipenem_0.06_T1 - None_0_T0",
	"Imipenem_0.06_T2 - None_0_T0",
	"Imipenem_0.09_T1 - None_0_T0",
	"Meropenem_0.11_T1 - None_0_T0",
	"Meropenem_0.11_T2 - None_0_T0",
	"Meropenem_0.17_T1 - None_0_T0",
	"Meropenem_0.17_T2 - None_0_T0"))

conditions.interest <- fread("interest.tsv", sep = "\t")

conditions.all <- rbind(conditions.full, conditions.interest) %>% unique

median_melted_results <- fread("Results/median_melted_results.tsv.gz")


depleted.long <- 
	conditions.all %>% 
	inner_join(
		median_melted_results %>%
			filter(condition %in% c(
				"None_0_T2 - None_0_T0",
				"Meropenem_0.11_T2 - None_0_T2",
				"Imipenem_0.06_T2 - None_0_T2") &
					type == "perfect")) 

mixed <- depleted.long %>%
	mutate(
		depleted = (medLFC < 0 & FDR < 0.0001), 
		resistant = (medLFC > 0 & FDR < 0.0001)) %>% 
	select(unique_name, condition, depleted, resistant) %>% 
	pivot_longer(!c(condition, unique_name)) %>% 
	mutate(condition = gsub(" - None_0_T2", "", condition)) %>%
	mutate(condition = gsub(" - None_0_T0", "(induced)", condition)) %>% 
	mutate(condition = gsub("_", " ", condition)) %>% 
	mutate(condition = gsub("T2", "", condition)) %>% 
	unite("condition", c("condition", "name"), sep = " ") %>% 
	pivot_wider(id_cols = unique_name, names_from = condition, values_from = value)

## distinct mode
mixed.matrix.distinct <- mixed %>%
	make_comb_mat(mode = "distinct")

p <- UpSet(
	mixed.matrix.distinct, 
	top_annotation = HeatmapAnnotation(
		"Distinct Sets" = anno_barplot(
			comb_size(mixed.matrix.distinct),
			border = FALSE, 
			height = unit(8, "cm"),
			add_numbers = T), 
annotation_name_gp = gpar(fontsize = 10),		
show_annotation_name = TRUE),
	right_annotation = rowAnnotation(
		"Genes Significant Beyond Induction at T2" = anno_barplot(
			set_size(mixed.matrix.distinct),
			border = FALSE,
			gp = gpar(
				fill = viridis(5, direction = -1, alpha = 0.35, option = "plasma") %>% rep(each = 2), 
				col = c("red", "black", "red", "black"),
				lwd = 3,
				lty = c("dotted", "solid", "dotted", "solid")),
			width = unit(6, "cm"),
			add_numbers = T),
		annotation_name_gp = gpar(fontsize = 8)),
	set_order = set_name(mixed.matrix.distinct),
	row_names_gp = grid::gpar(fontsize = 10),
	comb_col = darken(viridis(5, direction = -1)[comb_degree(mixed.matrix.distinct)], 0.35))

print(p)


## intersect mode
mixed.matrix.intersect <- mixed %>%
	make_comb_mat(mode = "intersect")

p <- UpSet(
	mixed.matrix.intersect, 
	top_annotation = HeatmapAnnotation(
		"Intersect Sets" = anno_barplot(
			comb_size(mixed.matrix.intersect),
			border = FALSE, 
			height = unit(8, "cm"),
			add_numbers = T), 
		annotation_name_gp = gpar(fontsize = 10),		
		show_annotation_name = TRUE),
	right_annotation = rowAnnotation(
		"Genes Significant Beyond Induction at T2" = anno_barplot(
			set_size(mixed.matrix.intersect),
			border = FALSE,
			gp = gpar(
				fill = viridis(5, direction = -1, alpha = 0.35, option = "plasma") %>% rep(each = 2), 
				col = c("red", "black", "red", "black"),
				lwd = 3,
				lty = c("dotted", "solid", "dotted", "solid")),
			width = unit(6, "cm"),
			add_numbers = T),
		annotation_name_gp = gpar(fontsize = 8)),
	set_order = set_name(mixed.matrix.intersect),
	row_names_gp = grid::gpar(fontsize = 10),
	comb_col = darken(viridis(5, direction = -1)[comb_degree(mixed.matrix.intersect)], 0.35))

print(p)