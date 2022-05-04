library(pacman)

p_load(pheatmap, data.table, ggplot2)

melted_results <- fread("Results/melted_results.tsv.gz", sep = "\t")

interest <- fread("Results/interest.tsv", sep = "\t")

generate_breaks = function(x, n, center = F){
	if(center){
		m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
		res = seq(-m, m, length.out = n + 1)
	}
	else{
		res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
	}
	
	return(res)
}

selected_LFC <-
	melted_results[
		condition %in% c(
			"None_0_T1 - None_0_T0",
			"Colistin_0.44_T1 - None_0_T0",
			"None_0_T2 - None_0_T0",
			"Colistin_0.44_T2 - None_0_T0") 
		& Pathway == "LPS"]

median_LFC <- dcast(
	selected_LFC [type == "perfect"], 
	AB030 + unique_name + Pathway ~ condition, 
	value.var = "LFC", 
	fun.aggregate = median)

# change order
median_LFC <- median_LFC[, .(
	AB030,
	unique_name,
	Pathway,
	`None_0_T1 - None_0_T0`,
	`Colistin_0.44_T1 - None_0_T0`,
	`None_0_T2 - None_0_T0`,
	`Colistin_0.44_T2 - None_0_T0`)]

LFC_grid <- data.matrix(
	median_LFC[, c(-1:-3)])

rownames(LFC_grid) <- median_LFC[, paste(unique_name)]

# format = Pathway + italic(Gene Name)
pathway_gene <- lapply(
	median_LFC[, .I],
	function(x) bquote(.(median_LFC[x]$Pathway) ~ italic(.(median_LFC[x]$unique_name))))

formatted_conditions <- c(
	bquote(Induction ~ (t[1])),
	bquote(Induction + Colistin ~ (t[1])),
	bquote(Induction ~ (t[2])),
	bquote(Induction + Colistin ~ (t[2])))

plot_matrix <- LFC_grid

breaks <- 
	generate_breaks(plot_matrix, n = 256, center = T)

plot_colors <- 
	viridis(n = 256, option = "viridis")

plot_colors <- colorRampPalette(c("#ba000d", "white", "#007ac1"), interpolate = "linear")(256)

to_plot_title <- 
	paste("LPS Biosynthetic Genes")

to_plot <- pheatmap(plot_matrix,
										col = plot_colors,
										breaks = breaks,
										border_color = NA,
										cellwidth = 48,
										cellheight = 16,
										main = to_plot_title,
										cluster_cols = F,
										angle_col = 45,
										clustering_method = "ward.D2",
										labels_row = as.expression(pathway_gene),
										labels_col = as.expression(formatted_conditions),
										gaps_col = 2)

print(to_plot)



