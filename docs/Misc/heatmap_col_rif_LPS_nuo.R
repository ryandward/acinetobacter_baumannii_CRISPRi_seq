# source('../meta_counter.R')

library(pacman)

p_load(pheatmap, data.table, ggplot2)

melted_results <- fread("Results/melted_results.tsv.gz", sep = "\t")

interest <- fread("Results/interest.tsv", sep = "\t")

scale_rows = function(x){
	m = apply(x, 1, mean, na.rm = T)
	s = apply(x, 1, sd, na.rm = T)
	return((x - m) / s)
}

scale_mat = function(mat, scale){
	if(!(scale %in% c("none", "row", "column"))){
		stop("scale argument shoud take values: 'none', 'row' or 'column'")
	}
	mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
	return(mat)
}

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
		shift %like% "T2"
		& base %like% "T2"
		& !shift %like% "None" 
		& !Pathway %like% "Ribosome" 
		& !Pathway %like% "Other" 
		& (shift %like% "Col" | shift %like% "Rif")
		& (base %like% "None_0")
		& !Pathway  %like% "tRNA"
		& !Pathway %like% "Cell Wall"
		& type == "perfect"]


median_LFC <- dcast(
	selected_LFC, 
	AB030 + unique_name + Pathway ~ condition, 
	value.var = "LFC", 
	fun.aggregate = median)


LFC_grid <- data.matrix(
	median_LFC[, c(-1:-3)])

rownames(LFC_grid) <- median_LFC[, paste(unique_name)]

# format = Pathway + italic(Gene Name)
newnames <- lapply(
	median_LFC[, .I],
	function(x) bquote(.(median_LFC[x]$Pathway) ~ italic(.(median_LFC[x]$unique_name))))

plot_matrix <- LFC_grid


break_halves <- 
	length(unique(as.vector(plot_matrix)))

breaks <- c(
	seq(min(plot_matrix), 0, length.out = break_halves)[-break_halves],
	seq(0, max(plot_matrix), length.out = break_halves))

plot_colors <- c(
	colorRampPalette(c("#ba000d", "white"))(break_halves)[-break_halves],
	colorRampPalette(c("white", "#007ac1"))(break_halves)[-1])



# breaks <- generate_breaks(plot_matrix, n = 256, center = T)
# 
# plot_colors <- viridis(n = 256, option = "viridis")
# 
# plot_colors <- colorRampPalette(c("#ba000d", "white", "#007ac1"), interpolate = "linear")(256)

to_plot_title <- paste("Log Fold Change Fitness by Pathway")

to_plot <- pheatmap(plot_matrix,
										col = plot_colors,
										breaks = breaks,
										border_color = NA,
										cellwidth = 45,
										cellheight = 14,
										cutree_rows = 2,
										cutree_cols = 2,
										main = to_plot_title,
										angle_col = 45,
										fontsize_col = 16,
										fontsize_row = 12,
										clustering_method = "ward.D2",
										clustering_distance_rows = "maximum",
										clustering_distance_cols = "maximum",
										labels_row = as.expression(newnames))

print(to_plot)



